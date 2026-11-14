# -*- coding: utf-8 -*-
"""Independent replicate (trajectory 2) analysis and two-trajectory comparison.

Reuses the exact RMSD / hydrogen-bond / distance and single-trajectory MM-GBSA
(OBC2, NoCutoff, entropy omitted) logic of e39_md_analysis.py and e40_mmgbsa.py.

Inputs
  - trajectory 2: CHAIHU_MD2_DIR (default C:\\dockwork\\md2), traj2.dcd
  - trajectory 1 summary CSVs in the calculation root (MD_RMSD.csv,
    MD_hydrogen_bonds.csv, MD_hydrogen_bond_occupancy.csv, MMGBSA_summary.csv)
Outputs
  - md2_* per-frame/summary CSVs in the calculation root
  - md2_replicate_comparison.json in data/processed_small (repo)
"""
import collections
import json
import os
import pickle
import sys

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

MD1 = r"C:\dockwork\md"
MD2 = os.environ.get("CHAIHU_MD2_DIR", r"C:\dockwork\md2")
ROOT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
REPO_SMALL = os.path.join(
    ROOT, "chaihu-nociceptive-mapping", "data", "processed_small"
)
SKIP = 6
STRIDE = 2

# ---------------------------------------------------------------------------
# Load topology
# ---------------------------------------------------------------------------
m = pickle.load(open(os.path.join(MD2, "modeller.pkl"), "rb"))
TOP = m["topology"]
N_PROT = int(m["n_prot"])
LIG = list(m["lig_idx"])
top = md.Topology.from_openmm(TOP)
backbone = top.select("backbone and chainid 0")
ca = top.select("name CA and chainid 0")
prot_atoms = list(range(0, N_PROT))


def is_heavy(i):
    return top.atom(i).element.symbol != "H"


lig_heavy = [i for i in LIG if is_heavy(i)]

tr_all = md.load(os.path.join(MD2, "traj2.dcd"), top=top)
prod_idx = list(range(SKIP, tr_all.n_frames))
tr = tr_all.slice(prod_idx)
t_ns = (np.arange(tr_all.n_frames) * 0.04)[SKIP:]
frame_labels = np.array(prod_idx)
print("total frames:", tr_all.n_frames, "| production frames:", tr.n_frames)

# ---------------------------------------------------------------------------
# 1) RMSD
# ---------------------------------------------------------------------------
ref = tr.slice(0)
bb_rmsd = md.rmsd(tr, ref, atom_indices=backbone, parallel=False)
tr.superpose(ref, atom_indices=backbone)


def rmsd_noalign(coords, refc):
    d = coords - refc
    return np.sqrt((d * d).sum(axis=(1, 2)) / coords.shape[1])


lig_rmsd = rmsd_noalign(tr.xyz[:, lig_heavy, :], tr.xyz[0, lig_heavy, :])
pd.DataFrame(
    {
        "frame": frame_labels,
        "time_ns": t_ns,
        "protein_backbone_RMSD_nm": bb_rmsd,
        "ligand_heavy_RMSD_nm": lig_rmsd,
    }
).to_csv(os.path.join(ROOT, "md2_MD_RMSD.csv"), index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------------------
# 2) Hydrogen bonds
# ---------------------------------------------------------------------------
solute = tr.atom_slice(prot_atoms + LIG)
lig_set = set(LIG)
per_frame_count = np.zeros(tr.n_frames, dtype=int)
occ_counter = collections.Counter()
for fi in range(tr.n_frames):
    hb = md.baker_hubbard(solute.slice(fi))
    keys_this = set()
    for donor, hh, acc in hb:
        donor, hh, acc = int(donor), int(hh), int(acc)
        d_in = donor in lig_set
        a_in = acc in lig_set
        if d_in == a_in:
            continue
        per_frame_count[fi] += 1
        prot_idx = acc if d_in else donor
        r = top.atom(prot_idx).residue
        keys_this.add((r.resSeq, r.name, "ligDonor" if d_in else "ligAcceptor"))
    for k in keys_this:
        occ_counter[k] += 1
pd.DataFrame(
    {
        "frame": frame_labels,
        "time_ns": t_ns,
        "n_ligand_protein_HB": per_frame_count,
    }
).to_csv(
    os.path.join(ROOT, "md2_MD_hydrogen_bonds.csv"),
    index=False,
    encoding="utf-8-sig",
)
occ_df = pd.DataFrame(
    [
        (rn, nm, kind, c, round(c / tr.n_frames, 3))
        for (rn, nm, kind), c in occ_counter.items()
    ],
    columns=["resSeq", "residue", "ligand_role", "n_frames", "occupancy"],
).sort_values("occupancy", ascending=False)
occ_df.to_csv(
    os.path.join(ROOT, "md2_MD_hydrogen_bond_occupancy.csv"),
    index=False,
    encoding="utf-8-sig",
)

# ---------------------------------------------------------------------------
# 3) Ligand-protein distance
# ---------------------------------------------------------------------------
prot_heavy = [i for i in prot_atoms if is_heavy(i)]
min_dist = np.zeros(tr.n_frames)
com_dist = np.zeros(tr.n_frames)
lc = solute.xyz[:, lig_heavy, :].mean(axis=1)
pc = solute.xyz[:, prot_heavy, :].mean(axis=1)
for fi in range(tr.n_frames):
    a = solute.xyz[fi, lig_heavy, :]
    b = solute.xyz[fi, prot_heavy, :]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    min_dist[fi] = d.min()
    com_dist[fi] = np.linalg.norm(lc[fi] - pc[fi])
pd.DataFrame(
    {
        "frame": frame_labels,
        "time_ns": t_ns,
        "min_lig_prot_distance_nm": min_dist,
        "lig_prot_COM_distance_nm": com_dist,
    }
).to_csv(
    os.path.join(ROOT, "md2_MD_lig_prot_distance.csv"),
    index=False,
    encoding="utf-8-sig",
)

# ---------------------------------------------------------------------------
# 4) MM-GBSA (single-trajectory, OBC2, stride 2)
# ---------------------------------------------------------------------------
with open(os.path.join(MD2, "system.xml")) as f:
    FULL = XmlSerializer.deserialize(f.read())
ELEMS = [a.element.symbol for a in list(TOP.atoms())]

_bondmap = {}
for _fc in FULL.getForces():
    if isinstance(_fc, HarmonicBondForce):
        for _i in range(_fc.getNumBonds()):
            _a, _b, _, _ = _fc.getBondParameters(_i)
            _a, _b = int(_a), int(_b)
            _bondmap.setdefault(_a, []).append(_b)
            _bondmap.setdefault(_b, []).append(_a)


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
    nb_src = next(f for f in FULL.getForces() if isinstance(f, NonbondedForce))
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
                    af.addAngle(old2new[a], old2new[b], old2new[c], ang, k)
    tf = PeriodicTorsionForce()
    tf.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc, PeriodicTorsionForce):
            for i in range(fc.getNumTorsions()):
                a, b, c, d, per, ph, k = fc.getTorsionParameters(i)
                a, b, c, d = int(a), int(b), int(c), int(d)
                if a in keep and b in keep and c in keep and d in keep:
                    tf.addTorsion(old2new[a], old2new[b], old2new[c], old2new[d], per, ph, k)
    nb = NonbondedForce()
    nb.setNonbondedMethod(NonbondedForce.NoCutoff)
    nb.setForceGroup(2)
    for old in sorted(keep):
        q, sig, eps = nb_src.getParticleParameters(old)
        nb.addParticle(
            q.value_in_unit(CH) * CH,
            sig.value_in_unit(NM) * NM,
            eps.value_in_unit(KJ) * KJ,
        )
    for j in range(nb_src.getNumExceptions()):
        a, b, q, sig, eps = nb_src.getExceptionParameters(j)
        a, b = int(a), int(b)
        if a in keep and b in keep:
            nb.addException(
                old2new[a],
                old2new[b],
                q.value_in_unit(CH2) * CH2,
                sig.value_in_unit(NM) * NM,
                eps.value_in_unit(KJ) * KJ,
            )
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
    ctx = Context(sysn, VerletIntegrator(0.001), Platform.getPlatformByName("CPU"))
    subs[name] = (ctx, kept)
    print(name, "atoms:", len(kept))

frames = list(range(SKIP, tr_all.n_frames, STRIDE))
rows = []
for fi in frames:
    rec = {"frame": fi, "time_ns": round(fi * 0.04, 3)}
    for name, (ctx, kept) in subs.items():
        ctx.setPositions([Vec3(*tr_all.xyz[fi, a]) for a in kept])
        e_bond = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(KJ)
        e_nb = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(KJ)
        e_gb = ctx.getState(getEnergy=True, groups={3}).getPotentialEnergy().value_in_unit(KJ)
        rec[name + "_bonded"] = e_bond
        rec[name + "_MM"] = e_nb
        rec[name + "_GBSA"] = e_gb
        rec[name + "_total"] = e_bond + e_nb + e_gb
    for comp in ["bonded", "MM", "GBSA", "total"]:
        d = rec["complex_" + comp] - rec["protein_" + comp] - rec["ligand_" + comp]
        rec["dG_" + comp + "_kJ"] = d
        rec["dG_" + comp + "_kcal"] = d / 4.184
    rows.append(rec)
pf = pd.DataFrame(rows)
pf.to_csv(os.path.join(ROOT, "md2_MMGBSA_per_frame.csv"), index=False, encoding="utf-8-sig")
summ = []
for comp in ["bonded", "MM", "GBSA", "total"]:
    x = pf["dG_" + comp + "_kcal"]
    summ.append(
        {
            "component": comp,
            "mean_kcal_mol": x.mean(),
            "SD": x.std(ddof=1),
            "SEM": x.std(ddof=1) / np.sqrt(len(x)),
            "median": x.median(),
        }
    )
sf = pd.DataFrame(summ)
sf.to_csv(os.path.join(ROOT, "md2_MMGBSA_summary.csv"), index=False, encoding="utf-8-sig")
print(sf.to_string(index=False))

# ---------------------------------------------------------------------------
# 5) Two-trajectory comparison
# ---------------------------------------------------------------------------
md1_rmsd = pd.read_csv(os.path.join(ROOT, "MD_RMSD.csv"), encoding="utf-8-sig")
md1_hb = pd.read_csv(os.path.join(ROOT, "MD_hydrogen_bonds.csv"), encoding="utf-8-sig")
md1_occ = pd.read_csv(
    os.path.join(ROOT, "MD_hydrogen_bond_occupancy.csv"), encoding="utf-8-sig"
)
md1_gbsa = pd.read_csv(os.path.join(ROOT, "MMGBSA_summary.csv"), encoding="utf-8-sig")
md1_total = md1_gbsa.loc[md1_gbsa.component == "total"].iloc[0]
md2_total = sf.loc[sf.component == "total"].iloc[0]

md1_pf = pd.read_csv(os.path.join(ROOT, "MMGBSA_per_frame.csv"), encoding="utf-8-sig")
common_n = min(len(md1_pf), len(pf))
g1 = md1_pf["dG_total_kcal"].to_numpy()[:common_n]
g2 = pf["dG_total_kcal"].to_numpy()[:common_n]

md1_top = md1_occ.set_index(["resSeq", "residue", "ligand_role"])["occupancy"]
md2_top = occ_df.set_index(["resSeq", "residue", "ligand_role"])["occupancy"]
shared_keys = sorted(set(md1_top.index) & set(md2_top.index))
hb_occ_compare = [
    {
        "resSeq": int(k[0]),
        "residue": k[1],
        "ligand_role": k[2],
        "traj1_occupancy": float(md1_top[k]),
        "traj2_occupancy": float(md2_top[k]),
    }
    for k in shared_keys
]

comparison = {
    "n_production_frames": {"traj1": int(len(md1_rmsd)), "traj2": int(tr.n_frames)},
    "mean_backbone_RMSD_nm": {
        "traj1": float(md1_rmsd["protein_backbone_RMSD_nm"].mean()),
        "traj2": float(bb_rmsd.mean()),
    },
    "mean_ligand_RMSD_nm": {
        "traj1": float(md1_rmsd["ligand_heavy_RMSD_nm"].mean()),
        "traj2": float(lig_rmsd.mean()),
    },
    "mean_HB_count": {
        "traj1": float(md1_hb["n_ligand_protein_HB"].mean()),
        "traj2": float(per_frame_count.mean()),
    },
    "MMGBSA_total_kcal_mol": {
        "traj1": {
            "mean": float(md1_total["mean_kcal_mol"]),
            "SD": float(md1_total["SD"]),
            "SEM": float(md1_total["SEM"]),
            "median": float(md1_total["median"]),
        },
        "traj2": {
            "mean": float(md2_total["mean_kcal_mol"]),
            "SD": float(md2_total["SD"]),
            "SEM": float(md2_total["SEM"]),
            "median": float(md2_total["median"]),
        },
        "mean_difference_traj2_minus_traj1": float(
            md2_total["mean_kcal_mol"] - md1_total["mean_kcal_mol"]
        ),
    },
    "frame_dG_pearson": float(np.corrcoef(g1, g2)[0, 1]),
    "shared_HB_occupancy": hb_occ_compare,
    "note": (
        "Two independent single trajectories (velocity seeds 42 and 2027); "
        "MM-GBSA omits entropy and is relative, not absolute affinity."
    ),
}
out_json = os.path.join(REPO_SMALL, "md2_replicate_comparison.json")
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(comparison, f, ensure_ascii=False, indent=2)
print("wrote", out_json)
