# -*- coding: utf-8 -*-
"""WP10 multi-pose MM-GBSA across ten 20 ns arms.

Single-trajectory MM-GBSA (OBC2 implicit solvent, NoCutoff, entropy omitted)
uses the same subsystem construction as 10_wp7_mmgbsa.py. All ten arms share
an identical solute topology (protein atoms 0..n_prot-1, contiguous ligand),
so the three energy contexts (complex, protein, ligand) are built once from
arm0's system.xml and reused; per-frame positions are taken from each arm's
DCD. Frames are sampled at a fixed stride to control cost.

Convergence is reported per arm (cumulative running mean, 4 ns time-block
means) and across arms, where the between-arm spread of the binding estimate
directly quantifies pose dependence.

Inputs
  - arm0 system.xml; every arm modeller.pkl and arm{i}.dcd
Outputs (repo data/processed_small/wp10)
  - wp10_mmgbsa_per_frame.csv, wp10_mmgbsa_summary.csv,
    wp10_mmgbsa_block_convergence.csv, wp10_mmgbsa_cross_arm.json
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

ROOT = r"C:\dockwork\md\wp10_multipose"
REPO = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\chaihu-nociceptive-mapping"
OUT = os.path.join(REPO, "data", "processed_small", "wp10")
os.makedirs(OUT, exist_ok=True)

N_ARMS = 10
SKIP = 6
STRIDE = 4              # sample every 0.16 ns -> 125 frames per arm
BLOCK_NS = 4.0          # 5 blocks per 20 ns arm
FRAME_DT_NS = 0.04

# Build energy contexts once from arm0; solute topology is shared.
with open(os.path.join(ROOT, "arm0", "modeller.pkl"), "rb") as fh:
    blob0 = pickle.load(fh)
N_PROT = int(blob0["n_prot"])
LIG0 = list(blob0["lig_idx"])
with open(os.path.join(ROOT, "arm0", "system.xml")) as fh:
    FULL = XmlSerializer.deserialize(fh.read())

_bondmap = {}
for fc in FULL.getForces():
    if isinstance(fc, HarmonicBondForce):
        for i in range(fc.getNumBonds()):
            a, b, _, _ = fc.getBondParameters(i)
            a, b = int(a), int(b)
            _bondmap.setdefault(a, []).append(b)
            _bondmap.setdefault(b, []).append(a)

ELEMS = [a.element.symbol for a in blob0["topology"].atoms()]


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
    "complex": list(range(0, LIG0[-1] + 1)),
    "protein": list(range(0, N_PROT)),
    "ligand": LIG0,
}
subs = {}
for name, keep in species.items():
    sysn, kept = build_subsystem(keep)
    ctx = Context(sysn, VerletIntegrator(0.001),
                  Platform.getPlatformByName("CPU"))
    subs[name] = (ctx, kept)
    print(name, "atoms:", len(kept))

records = []
for arm in range(N_ARMS):
    with open(os.path.join(ROOT, f"arm{arm}", "modeller.pkl"), "rb") as fh:
        blob = pickle.load(fh)
    top = md.Topology.from_openmm(blob["topology"])
    dcd = os.path.join(ROOT, f"arm{arm}", f"arm{arm}.dcd")
    tr_all = md.load(dcd, top=top)
    frames = list(range(SKIP, tr_all.n_frames, STRIDE))
    for fi in frames:
        rec = {"arm": arm, "frame": fi,
               "time_ns": round((fi - SKIP) * FRAME_DT_NS, 3)}
        for name, (ctx, kept) in subs.items():
            ctx.setPositions([Vec3(*tr_all.xyz[fi, a]) for a in kept])
            e_bond = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(KJ)
            e_mm = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(KJ)
            e_gb = ctx.getState(getEnergy=True, groups={3}).getPotentialEnergy().value_in_unit(KJ)
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
    print(f"arm{arm}: {len(frames)} sampled frames")
    del tr_all

pf = pd.DataFrame(records)
pf.to_csv(os.path.join(OUT, "wp10_mmgbsa_per_frame.csv"),
          index=False, encoding="utf-8-sig")

# Per-arm summary
summary_rows = []
for arm in range(N_ARMS):
    sub = pf[pf.arm == arm]
    for comp in ["bonded", "MM", "GBSA", "total"]:
        x = sub["dG_" + comp + "_kcal"]
        summary_rows.append({
            "arm": arm, "component": comp,
            "n_frames": len(x),
            "mean_kcal_mol": x.mean(),
            "SD": x.std(ddof=1),
            "naive_SEM": x.std(ddof=1) / np.sqrt(len(x)),
            "median": x.median(),
        })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT, "wp10_mmgbsa_summary.csv"),
                  index=False, encoding="utf-8-sig")

# Time-block convergence (4 ns blocks)
block_rows = []
for arm in range(N_ARMS):
    sub = pf[pf.arm == arm].sort_values("time_ns")
    t = sub["time_ns"].to_numpy()
    g = sub["dG_total_kcal"].to_numpy()
    block_id = np.floor(t / BLOCK_NS).astype(int)
    for b in sorted(set(block_id)):
        m = g[block_id == b]
        block_rows.append({
            "arm": arm, "block": int(b),
            "block_start_ns": b * BLOCK_NS,
            "block_end_ns": (b + 1) * BLOCK_NS,
            "n_samples": len(m),
            "block_mean_dG_kcal_mol": m.mean(),
            "block_SD": m.std(ddof=1),
        })
block_df = pd.DataFrame(block_rows)
block_df.to_csv(os.path.join(OUT, "wp10_mmgbsa_block_convergence.csv"),
                index=False, encoding="utf-8-sig")

# Cross-arm pose-dependence spread on total binding dG
tot = summary_df[summary_df.component == "total"].reset_index(drop=True)
arm_means = tot["mean_kcal_mol"].to_numpy()
cross_arm = {
    "method": "single-trajectory MM-GBSA, OBC2, NoCutoff, entropy omitted; "
              "relative not absolute binding free energy",
    "sampling_stride_ns": STRIDE * FRAME_DT_NS,
    "samples_per_arm": int((pf.arm == 0).sum()),
    "arm_total_dG_kcal_mol": {
        int(r["arm"]): float(r["mean_kcal_mol"])
        for _, r in tot.iterrows()},
    "mean_of_arm_means_kcal_mol": float(arm_means.mean()),
    "between_arm_SD_kcal_mol": float(arm_means.std(ddof=1)),
    "between_arm_SEM_kcal_mol": float(arm_means.std(ddof=1) / np.sqrt(N_ARMS)),
    "min_arm_kcal_mol": float(arm_means.min()),
    "max_arm_kcal_mol": float(arm_means.max()),
    "range_kcal_mol": float(arm_means.max() - arm_means.min()),
    "block_size_ns": BLOCK_NS,
    "n_blocks_per_arm": 5,
    "first_half_mean_kcal_mol": float(
        pf[pf.time_ns < 10.0]["dG_total_kcal"].mean()),
    "second_half_mean_kcal_mol": float(
        pf[pf.time_ns >= 10.0]["dG_total_kcal"].mean()),
}
with open(os.path.join(OUT, "wp10_mmgbsa_cross_arm.json"), "w",
          encoding="utf-8") as fh:
    json.dump(cross_arm, fh, ensure_ascii=False, indent=2)

print("\n=== per-arm total dG (kcal/mol) ===")
print(tot[["arm", "mean_kcal_mol", "SD", "naive_SEM"]].to_string(index=False))
print("\nmean of arms %.2f, between-arm SD %.2f, range %.2f"
      % (cross_arm["mean_of_arm_means_kcal_mol"],
         cross_arm["between_arm_SD_kcal_mol"],
         cross_arm["range_kcal_mol"]))
print("first half %.2f vs second half %.2f"
      % (cross_arm["first_half_mean_kcal_mol"],
         cross_arm["second_half_mean_kcal_mol"]))
print("wrote", OUT)
