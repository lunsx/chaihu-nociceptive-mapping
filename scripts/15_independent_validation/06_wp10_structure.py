# -*- coding: utf-8 -*-
"""WP10 multi-pose structural analysis across ten 20 ns arms.

Each arm has its own modeller.pkl/DCD but an identical solute topology
(protein atoms 0..n_prot-1, contiguous ligand at n_prot..n_prot+30). Only the
solvent count differs. Metric definitions match 09_wp7_structure.py: protein
backbone and ligand heavy-atom RMSD, per-residue CA/backbone RMSF, ligand-
protein hydrogen bonds with occupancy, and ligand-protein minimum/COM
distances. A protein-aligned ligand centre-of-mass displacement is added to
quantify whether the ligand retains its starting site or drifts. Arms are
processed one at a time and released to keep peak memory low.

Inputs
  - C:\\dockwork\\md\\wp10_multipose\\arm{i}\\modeller.pkl, arm{i}.dcd
Outputs (repo data/processed_small/wp10)
  - per-arm rmsd/rmsf/hb/hb_occupancy/distance CSVs
  - wp10_cross_arm_summary.csv, wp10_hb_occupancy_compare.csv,
    wp10_arm_comparison.json
"""
import collections
import json
import os
import pickle

import numpy as np
import pandas as pd
import mdtraj as md

ROOT = r"C:\dockwork\md\wp10_multipose"
REPO = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\chaihu-nociceptive-mapping"
OUT = os.path.join(REPO, "data", "processed_small", "wp10")
os.makedirs(OUT, exist_ok=True)

N_ARMS = 10
SKIP = 6                 # first six DCD frames are the staged equilibration
FRAME_DT_NS = 0.04       # 40 ps per saved frame

per_frame = {}
rmsf_store = {}
occ_store = {}
frame_counts = {}
res_order = None

for arm in range(N_ARMS):
    arm_dir = os.path.join(ROOT, f"arm{arm}")
    with open(os.path.join(arm_dir, "modeller.pkl"), "rb") as fh:
        blob = pickle.load(fh)
    TOP = blob["topology"]
    N_PROT = int(blob["n_prot"])
    LIG = list(blob["lig_idx"])
    top = md.Topology.from_openmm(TOP)
    backbone = top.select("backbone and chainid 0")
    ca = top.select("name CA and chainid 0")
    prot_atoms = list(range(0, N_PROT))

    def is_heavy(i):
        return top.atom(i).element.symbol != "H"

    lig_heavy = [i for i in LIG if is_heavy(i)]
    if res_order is None:
        res_order = [(top.atom(int(a)).residue.resSeq,
                      top.atom(int(a)).residue.name) for a in ca]

    dcd = os.path.join(arm_dir, f"arm{arm}.dcd")
    tr_all = md.load(dcd, top=top)
    prod_idx = list(range(SKIP, tr_all.n_frames))
    tr = tr_all.slice(prod_idx)
    n_frames = tr.n_frames
    frame_counts[arm] = n_frames
    time_ns = (np.array(prod_idx) - SKIP) * FRAME_DT_NS
    print(f"arm{arm}: total {tr_all.n_frames} frames, production {n_frames}")

    # 1) RMSD (protein alignment) -------------------------------------------
    ref = tr.slice(0)
    bb_rmsd = md.rmsd(tr, ref, atom_indices=backbone, parallel=False)
    tr.superpose(ref, atom_indices=backbone)

    def rmsd_noalign(coords, refc):
        d = coords - refc
        return np.sqrt((d * d).sum(axis=(1, 2)) / coords.shape[1])

    lig_rmsd = rmsd_noalign(tr.xyz[:, lig_heavy, :],
                            tr.xyz[0, lig_heavy, :])

    # 1b) ligand COM displacement in the protein-aligned frame (site drift)
    lig_com = tr.xyz[:, lig_heavy, :].mean(axis=1)
    lig_com_drift = np.linalg.norm(lig_com - lig_com[0], axis=1)

    pd.DataFrame({
        "frame": prod_idx,
        "time_ns": time_ns,
        "protein_backbone_RMSD_nm": bb_rmsd,
        "ligand_heavy_RMSD_nm": lig_rmsd,
        "ligand_COM_drift_nm": lig_com_drift,
    }).to_csv(os.path.join(OUT, f"wp10_arm{arm}_rmsd.csv"),
              index=False, encoding="utf-8-sig")

    # 2) RMSF ----------------------------------------------------------------
    xyz = tr.xyz
    ca_rmsf = np.sqrt(
        ((xyz[:, ca, :] - xyz[:, ca, :].mean(axis=0)) ** 2)
        .sum(axis=2)).mean(axis=0)
    res_ids, res_names, bb_rmsf = [], [], []
    for a in ca:
        res = top.atom(int(a)).residue
        bb_in_res = [i for i in backbone
                     if top.atom(int(i)).residue == res]
        d = xyz[:, bb_in_res, :] - xyz[:, bb_in_res, :].mean(axis=0)
        bb_rmsf.append(np.sqrt((d * d).sum(axis=2)).mean())
        res_ids.append(res.resSeq)
        res_names.append(res.name)
    pd.DataFrame({
        "resSeq": res_ids,
        "residue": res_names,
        "CA_RMSF_nm": ca_rmsf,
        "backbone_RMSF_nm": bb_rmsf,
    }).to_csv(os.path.join(OUT, f"wp10_arm{arm}_rmsf.csv"),
              index=False, encoding="utf-8-sig")
    rmsf_store[arm] = ca_rmsf

    # 3) Hydrogen bonds ------------------------------------------------------
    solute = tr.atom_slice(prot_atoms + LIG)
    stop = solute.topology
    lig_set = set(range(N_PROT, N_PROT + len(LIG)))
    hb_count = np.zeros(n_frames, dtype=int)
    occ_counter = collections.Counter()
    for fi in range(n_frames):
        hb = md.baker_hubbard(solute.slice(fi))
        keys_this = set()
        for donor, hh, acc in hb:
            donor, acc = int(donor), int(acc)
            d_in = donor in lig_set
            a_in = acc in lig_set
            if d_in == a_in:
                continue
            hb_count[fi] += 1
            prot_idx = acc if d_in else donor
            r = stop.atom(prot_idx).residue
            keys_this.add(
                (r.resSeq, r.name, "ligDonor" if d_in else "ligAcceptor"))
        for key in keys_this:
            occ_counter[key] += 1
    pd.DataFrame({
        "frame": prod_idx,
        "time_ns": time_ns,
        "n_ligand_protein_HB": hb_count,
    }).to_csv(os.path.join(OUT, f"wp10_arm{arm}_hb.csv"),
              index=False, encoding="utf-8-sig")
    occ_df = pd.DataFrame(
        [(rn, nm, role, c, round(c / n_frames, 3))
         for (rn, nm, role), c in occ_counter.items()],
        columns=["resSeq", "residue", "ligand_role",
                 "n_frames", "occupancy"])
    occ_df.to_csv(os.path.join(OUT, f"wp10_arm{arm}_hb_occupancy.csv"),
                  index=False, encoding="utf-8-sig")
    occ_store[arm] = {k: c / n_frames for k, c in occ_counter.items()}

    # 4) Ligand-protein distance ---------------------------------------------
    prot_heavy = [i for i in prot_atoms if is_heavy(i)]
    min_dist = np.zeros(n_frames)
    com_dist = np.zeros(n_frames)
    lc = solute.xyz[:, lig_heavy, :].mean(axis=1)
    pc = solute.xyz[:, prot_heavy, :].mean(axis=1)
    for fi in range(n_frames):
        a = solute.xyz[fi, lig_heavy, :]
        b = solute.xyz[fi, prot_heavy, :]
        d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        min_dist[fi] = d.min()
        com_dist[fi] = np.linalg.norm(lc[fi] - pc[fi])
    pd.DataFrame({
        "frame": prod_idx,
        "time_ns": time_ns,
        "min_lig_prot_distance_nm": min_dist,
        "lig_prot_COM_distance_nm": com_dist,
    }).to_csv(os.path.join(OUT, f"wp10_arm{arm}_distance.csv"),
              index=False, encoding="utf-8-sig")

    per_frame[arm] = {
        "time_ns": time_ns,
        "bb_rmsd": bb_rmsd,
        "lig_rmsd": lig_rmsd,
        "lig_com_drift": lig_com_drift,
        "hb": hb_count,
        "min_dist": min_dist,
        "com_dist": com_dist,
    }
    del tr_all, tr, solute
    print(f"arm{arm}: done and released")

# ---------------------------------------------------------------------------
# Cross-arm summary
# ---------------------------------------------------------------------------
METRICS = [
    ("bb_rmsd", "backbone_RMSD_nm"),
    ("lig_rmsd", "ligand_RMSD_nm"),
    ("lig_com_drift", "ligand_COM_drift_nm"),
    ("hb", "HB_count"),
    ("min_dist", "min_distance_nm"),
    ("com_dist", "COM_distance_nm"),
]
summary_rows = []
for arm in range(N_ARMS):
    pf = per_frame[arm]
    half = len(pf["bb_rmsd"]) // 2
    label = "canonical" if arm == 0 else f"blind_pose_{arm}"
    row = {"arm": arm, "label": label,
           "n_production_frames": frame_counts[arm]}
    for metric, mlabel in METRICS:
        x = pf[metric]
        row[f"mean_{mlabel}"] = float(x.mean())
        row[f"SD_{mlabel}"] = float(x.std(ddof=1))
        row[f"first_half_mean_{mlabel}"] = float(x[:half].mean())
        row[f"second_half_mean_{mlabel}"] = float(x[half:].mean())
    summary_rows.append(row)
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT, "wp10_cross_arm_summary.csv"),
                  index=False, encoding="utf-8-sig")

# Between-arm dispersion for each headline metric (pose-dependence scale)
dispersion = {}
for metric, mlabel in METRICS:
    means = np.array([per_frame[a][metric].mean() for a in range(N_ARMS)])
    dispersion[f"arm_mean_{mlabel}"] = {
        "mean_of_arms": float(means.mean()),
        "between_arm_SD": float(means.std(ddof=1)),
        "min": float(means.min()),
        "max": float(means.max()),
        "range": float(means.max() - means.min()),
    }

# Hydrogen-bond occupancy comparison on the union of contacts
all_keys = sorted(set().union(*[set(occ_store[a].keys())
                                for a in range(N_ARMS)]))
hb_rows = []
for key in all_keys:
    rn, nm, role = key
    occs = [occ_store[a].get(key, 0.0) for a in range(N_ARMS)]
    row = {"resSeq": rn, "residue": nm, "ligand_role": role}
    for a in range(N_ARMS):
        row[f"arm{a}_occupancy"] = occs[a]
    row["mean_occupancy"] = float(np.mean(occs))
    hb_rows.append(row)
pd.DataFrame(hb_rows).to_csv(
    os.path.join(OUT, "wp10_hb_occupancy_compare.csv"),
    index=False, encoding="utf-8-sig")

comparison = {
    "design": "ten 20 ns arms, one canonical plus nine blind sites; "
              "4 fs timestep, OpenCL; identical solute topology",
    "production_frames_each": int(frame_counts[0]),
    "between_arm_dispersion": dispersion,
    "per_arm_means": {
        f"arm{a}": {
            mlabel: float(per_frame[a][metric].mean())
            for metric, mlabel in METRICS}
        for a in range(N_ARMS)},
    "n_distinct_ligand_protein_HB_contacts": len(all_keys),
    "arm_mean_ligand_COM_drift_nm": {
        f"arm{a}": float(per_frame[a]["lig_com_drift"].mean())
        for a in range(N_ARMS)},
    "max_ligand_COM_drift_nm": {
        f"arm{a}": float(per_frame[a]["lig_com_drift"].max())
        for a in range(N_ARMS)},
}
with open(os.path.join(OUT, "wp10_arm_comparison.json"), "w",
          encoding="utf-8") as fh:
    json.dump(comparison, fh, ensure_ascii=False, indent=2)

print("\n=== between-arm dispersion ===")
for mlabel, d in dispersion.items():
    print(f"{mlabel}: mean {d['mean_of_arms']:.4f}, SD {d['between_arm_SD']:.4f}, "
          f"range {d['range']:.4f}")
print("distinct HB contacts:", len(all_keys))
print("wrote", OUT)
