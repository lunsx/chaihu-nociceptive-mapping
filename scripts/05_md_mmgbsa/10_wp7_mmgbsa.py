# -*- coding: utf-8 -*-
"""WP7 MM-GBSA cross-replicate estimates and convergence analysis.

Single-trajectory MM-GBSA (OBC2 implicit solvent, NoCutoff, entropy omitted)
is evaluated with the same subsystem construction as e40_mmgbsa.py and
08_md2_analysis.py for all three WP7 replicates. The three subsystems share an
identical topology across replicates, so the energy contexts are built once and
reused. Frames are sampled at a fixed stride to control cost.

Convergence is reported three ways: cumulative running mean within each
replicate, time-block means with a block-based standard error that accounts for
serial correlation, and the pooled mean as blocks accumulate across
replicates.

Inputs
  - C:\\dockwork\\md\\system.xml, modeller.pkl
  - C:\\dockwork\\md\\wp7\\rep{0,1,2}\\rep{0,1,2}.dcd
Outputs (repo data/processed_small/wp7)
  - wp7_mmgbsa_per_frame.csv, wp7_mmgbsa_summary.csv,
    wp7_mmgbsa_block_convergence.csv, wp7_mmgbsa_convergence.json
"""
import json
import os
import pickle

import numpy as np
import pandas as pd
import mdtraj as md

from openmm import (
    XmlSerializer,
    System,
    Platform,
    Vec3,
    Context,
    VerletIntegrator,
    HarmonicBondForce,
    HarmonicAngleForce,
    PeriodicTorsionForce,
    NonbondedForce,
    GBSAOBCForce,
    unit,
)

CH = unit.elementary_charge
CH2 = unit.elementary_charge ** 2
NM = unit.nanometer
KJ = unit.kilojoule_per_mole

WORK = r"C:\dockwork\md"
WP7 = os.path.join(WORK, "wp7")
REPO = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\chaihu-nociceptive-mapping"
OUT = os.path.join(REPO, "data", "processed_small", "wp7")
os.makedirs(OUT, exist_ok=True)

SKIP = 6
STRIDE = 4              # sample every 0.16 ns -> 125 frames per replicate
BLOCK_NS = 4.0          # 5 blocks per 20 ns replicate
FRAME_DT_NS = 0.04

with open(os.path.join(WORK, "modeller.pkl"), "rb") as fh:
    blob = pickle.load(fh)
TOP = blob["topology"]
N_PROT = int(blob["n_prot"])
LIG = list(blob["lig_idx"])
top = md.Topology.from_openmm(TOP)
ELEMS = [a.element.symbol for a in TOP.atoms()]

with open(os.path.join(WORK, "system.xml")) as fh:
    FULL = XmlSerializer.deserialize(fh.read())

_bondmap = {}
for fc in FULL.getForces():
    if isinstance(fc, HarmonicBondForce):
        for i in range(fc.getNumBonds()):
            a, b, _, _ = fc.getBondParameters(i)
            a, b = int(a), int(b)
            _bondmap.setdefault(a, []).append(b)
            _bondmap.setdefault(b, []).append(a)


def radius(i):
    e = ELEMS[i]
    if e == "C":
        return 0.17
    if e == "O":
        return 0.15
    if e == "N":
        return 0.155
    if e == "S":
        return 0.20
    if e == "H":
        heavy_nbr = [j for j in _bondmap.get(i, []) if ELEMS[j] != "H"]
        if heavy_nbr and ELEMS[heavy_nbr[0]] in ("O", "N"):
            return 0.10
        return 0.12
    return 0.15


def build_subsystem(keep):
    keep = set(keep)
    old2new = {old: k for k, old in enumerate(sorted(keep))}
    nb_src = next(f for f in FULL.getForces()
                  if isinstance(f, NonbondedForce))
    charges = {
        i: nb_src.getParticleParameters(i)[0].value_in_unit(CH)
        for i in range(nb_src.getNumParticles())
    }
    sysn = System()
    for old in sorted(keep):
        sysn.addParticle(FULL.getParticleMass(old))
    bf = HarmonicBondForce()
    bf.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc, HarmonicBondForce):
            for i in range(fc.getNumBonds()):
                a, b, L, k = fc.getBondParameters(i)
                a, b = int(a), int(b)
                if a in keep and b in keep:
                    bf.addBond(old2new[a], old2new[b], L, k)
    af = HarmonicAngleForce()
    af.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc, HarmonicAngleForce):
            for i in range(fc.getNumAngles()):
                a, b, c, ang, k = fc.getAngleParameters(i)
                a, b, c = int(a), int(b), int(c)
                if a in keep and b in keep and c in keep:
                    af.addAngle(old2new[a], old2new[b],
                                old2new[c], ang, k)
    tf = PeriodicTorsionForce()
    tf.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc, PeriodicTorsionForce):
            for i in range(fc.getNumTorsions()):
                a, b, c, d, per, ph, k = fc.getTorsionParameters(i)
                a, b, c, d = int(a), int(b), int(c), int(d)
                if a in keep and b in keep and c in keep and d in keep:
                    tf.addTorsion(old2new[a], old2new[b], old2new[c],
                                  old2new[d], per, ph, k)
    nb = NonbondedForce()
    nb.setNonbondedMethod(NonbondedForce.NoCutoff)
    nb.setForceGroup(2)
    for old in sorted(keep):
        q, sig, eps = nb_src.getParticleParameters(old)
        nb.addParticle(q.value_in_unit(CH) * CH,
                       sig.value_in_unit(NM) * NM,
                       eps.value_in_unit(KJ) * KJ)
    for j in range(nb_src.getNumExceptions()):
        a, b, q, sig, eps = nb_src.getExceptionParameters(j)
        a, b = int(a), int(b)
        if a in keep and b in keep:
            nb.addException(old2new[a], old2new[b],
                            q.value_in_unit(CH2) * CH2,
                            sig.value_in_unit(NM) * NM,
                            eps.value_in_unit(KJ) * KJ)
    gb = GBSAOBCForce()
    gb.setNonbondedMethod(GBSAOBCForce.NoCutoff)
    gb.setSolventDielectric(78.5)
    gb.setSoluteDielectric(1.0)
    gb.setSurfaceAreaEnergy(True)
    gb.setForceGroup(3)
    for old in sorted(keep):
        gb.addParticle(charges[old], radius(old), 0.7)
    for f in (bf, af, tf, nb, gb):
        sysn.addForce(f)
    return sysn, sorted(keep)


species = {
    "complex": list(range(0, LIG[-1] + 1)),
    "protein": list(range(0, N_PROT)),
    "ligand": LIG,
}
subs = {}
for name, keep in species.items():
    sysn, kept = build_subsystem(keep)
    ctx = Context(sysn, VerletIntegrator(0.001),
                  Platform.getPlatformByName("CPU"))
    subs[name] = (ctx, kept)
    print(name, "atoms:", len(kept))

records = []
for rep in range(3):
    dcd = os.path.join(WP7, f"rep{rep}", f"rep{rep}.dcd")
    tr_all = md.load(dcd, top=top)
    frames = list(range(SKIP, tr_all.n_frames, STRIDE))
    for fi in frames:
        rec = {"replicate": rep, "frame": fi,
               "time_ns": round((fi - SKIP) * FRAME_DT_NS, 3)}
        for name, (ctx, kept) in subs.items():
            ctx.setPositions([Vec3(*tr_all.xyz[fi, a]) for a in kept])
            e_bond = ctx.getState(
                getEnergy=True, groups={1}
            ).getPotentialEnergy().value_in_unit(KJ)
            e_mm = ctx.getState(
                getEnergy=True, groups={2}
            ).getPotentialEnergy().value_in_unit(KJ)
            e_gb = ctx.getState(
                getEnergy=True, groups={3}
            ).getPotentialEnergy().value_in_unit(KJ)
            rec[name + "_bonded"] = e_bond
            rec[name + "_MM"] = e_mm
            rec[name + "_GBSA"] = e_gb
            rec[name + "_total"] = e_bond + e_mm + e_gb
        for comp in ["bonded", "MM", "GBSA", "total"]:
            d = rec["complex_" + comp] - rec["protein_" + comp] \
                - rec["ligand_" + comp]
            rec["dG_" + comp + "_kJ"] = d
            rec["dG_" + comp + "_kcal"] = d / 4.184
        records.append(rec)
    print(f"rep{rep}: {len(frames)} sampled frames")
    del tr_all

pf = pd.DataFrame(records)
pf.to_csv(os.path.join(OUT, "wp7_mmgbsa_per_frame.csv"),
          index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------------------
# Per-replicate summary
# ---------------------------------------------------------------------------
summary_rows = []
for rep in range(3):
    sub = pf[pf.replicate == rep]
    for comp in ["bonded", "MM", "GBSA", "total"]:
        x = sub["dG_" + comp + "_kcal"]
        summary_rows.append({
            "replicate": f"rep{rep}",
            "component": comp,
            "n_frames": len(x),
            "mean_kcal_mol": x.mean(),
            "SD": x.std(ddof=1),
            "naive_SEM": x.std(ddof=1) / np.sqrt(len(x)),
            "median": x.median(),
        })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT, "wp7_mmgbsa_summary.csv"),
                  index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------------------
# Time-block convergence (4 ns blocks)
# ---------------------------------------------------------------------------
block_rows = []
for rep in range(3):
    sub = pf[pf.replicate == rep].sort_values("time_ns")
    t = sub["time_ns"].to_numpy()
    g = sub["dG_total_kcal"].to_numpy()
    block_id = np.floor(t / BLOCK_NS).astype(int)
    for b in sorted(set(block_id)):
        m = g[block_id == b]
        block_rows.append({
            "replicate": rep,
            "block": int(b),
            "block_start_ns": b * BLOCK_NS,
            "block_end_ns": (b + 1) * BLOCK_NS,
            "n_samples": len(m),
            "block_mean_dG_kcal_mol": m.mean(),
            "block_SD": m.std(ddof=1),
        })
block_df = pd.DataFrame(block_rows)
block_df.to_csv(os.path.join(OUT, "wp7_mmgbsa_block_convergence.csv"),
                index=False, encoding="utf-8-sig")

# block-based SEM across all 15 blocks (serial correlation handled by blocking)
all_block_means = block_df["block_mean_dG_kcal_mol"].to_numpy()
n_blocks = len(all_block_means)
blocked_sem = all_block_means.std(ddof=1) / np.sqrt(n_blocks)

# cumulative pooled mean as blocks accumulate in chronological order across reps
cum = []
running = []
for _, row in block_df.iterrows():
    running.append(row["block_mean_dG_kcal_mol"])
    cum.append({
        "replicate": int(row["replicate"]),
        "block": int(row["block"]),
        "cumulative_ns": row["block_end_ns"]
        if row["replicate"] == 0
        else row["block_end_ns"] + 20.0 * row["replicate"],
        "pooled_running_mean_kcal_mol": float(np.mean(running)),
    })
cum_df = pd.DataFrame(cum)

# within-replicate cumulative running mean at sampled frames
within = {}
for rep in range(3):
    sub = pf[pf.replicate == rep].sort_values("time_ns")
    g = sub["dG_total_kcal"].to_numpy()
    within[f"rep{rep}"] = {
        "time_ns": sub["time_ns"].tolist(),
        "cumulative_mean_kcal_mol": np.cumsum(g)
        / np.arange(1, len(g) + 1),
    }

total_summary = summary_df[summary_df.component == "total"]
rep_means = total_summary.set_index("replicate")["mean_kcal_mol"].to_dict()
pooled_mean = float(np.mean([rep_means[f"rep{r}"] for r in range(3)]))
between_sd = float(np.std([rep_means[f"rep{r}"] for r in range(3)],
                          ddof=1))

convergence = {
    "method": "single-trajectory MM-GBSA, OBC2, NoCutoff, entropy omitted; "
              "relative not absolute binding free energy",
    "sampling_stride_ns": STRIDE * FRAME_DT_NS,
    "samples_per_replicate": int((pf.replicate == 0).sum()),
    "replicate_total_dG_kcal_mol": {k: float(v) for k, v in rep_means.items()},
    "pooled_mean_kcal_mol": pooled_mean,
    "between_replicate_SD_kcal_mol": between_sd,
    "between_replicate_SEM_kcal_mol": between_sd / np.sqrt(3),
    "block_size_ns": BLOCK_NS,
    "n_blocks": int(n_blocks),
    "block_based_SEM_kcal_mol": float(blocked_sem),
    "naive_pooled_SEM_kcal_mol": float(
        total_summary["naive_SEM"].mean() / np.sqrt(3)),
    "first_half_mean_kcal_mol": float(
        pf[pf.time_ns < 10.0]["dG_total_kcal"].mean()),
    "second_half_mean_kcal_mol": float(
        pf[pf.time_ns >= 10.0]["dG_total_kcal"].mean()),
    "pooled_running_mean_final_kcal_mol": float(
        cum_df["pooled_running_mean_kcal_mol"].iloc[-1]),
    "within_rep_cumulative_final": {
        rep: float(within[rep]["cumulative_mean_kcal_mol"][-1])
        for rep in within},
}
with open(os.path.join(OUT, "wp7_mmgbsa_convergence.json"), "w",
          encoding="utf-8") as fh:
    json.dump(convergence, fh, ensure_ascii=False, indent=2)

print("\n=== per-replicate total dG (kcal/mol) ===")
print(total_summary[["replicate", "mean_kcal_mol", "SD",
                     "naive_SEM"]].to_string(index=False))
print("\npooled mean %.2f, between-rep SD %.2f, blocked SEM %.2f"
      % (pooled_mean, between_sd, blocked_sem))
print("first half %.2f vs second half %.2f"
      % (convergence["first_half_mean_kcal_mol"],
         convergence["second_half_mean_kcal_mol"]))
print("wrote", OUT)
