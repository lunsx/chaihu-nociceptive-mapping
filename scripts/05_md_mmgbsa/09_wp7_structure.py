# -*- coding: utf-8 -*-
"""WP7 structural robustness: cross-replicate trajectory analysis.

Three independent 20 ns replicates (distinct velocity seeds) are analysed with
the same metric definitions used in e39_md_analysis.py and 08_md2_analysis.py:
protein backbone and ligand heavy RMSD, per-residue CA and backbone RMSF,
ligand-protein hydrogen bonds with occupancy, and ligand-protein minimum and
centre-of-mass distances. Replicates are processed one at a time and released
to keep peak memory low on the 16 GB workstation.

Inputs
  - topology: C:\\dockwork\\md\\modeller.pkl
  - trajectories: C:\\dockwork\\md\\wp7\\rep{0,1,2}\\rep{0,1,2}.dcd
Outputs (repo data/processed_small/wp7)
  - per-replicate rmsd / rmsf / hb / hb_occupancy / distance CSVs
  - wp7_cross_replicate_summary.csv, wp7_rmsf_cross_replicate.csv,
    wp7_hb_occupancy_compare.csv and wp7_replicate_comparison.json
"""
import collections
import json
import os
import pickle

import numpy as np
import pandas as pd
import mdtraj as md

WORK = r"C:\dockwork\md"
WP7 = os.path.join(WORK, "wp7")
REPO = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\chaihu-nociceptive-mapping"
OUT = os.path.join(REPO, "data", "processed_small", "wp7")
os.makedirs(OUT, exist_ok=True)

N_REPS = 3
SKIP = 6                 # first six DCD frames are the staged equilibration
FRAME_DT_NS = 0.04       # 40 ps per saved frame

with open(os.path.join(WORK, "modeller.pkl"), "rb") as fh:
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

# residue order follows the CA atoms; used to align RMSF across replicates
res_order = [(top.atom(int(a)).residue.resSeq,
              top.atom(int(a)).residue.name) for a in ca]

per_frame = {}          # replicate -> small per-frame arrays
rmsf_store = {}         # replicate -> CA RMSF aligned to res_order
occ_store = {}          # replicate -> (resSeq,residue,role) -> occupancy
frame_counts = {}

for rep in range(N_REPS):
    dcd = os.path.join(WP7, f"rep{rep}", f"rep{rep}.dcd")
    tr_all = md.load(dcd, top=top)
    prod_idx = list(range(SKIP, tr_all.n_frames))
    tr = tr_all.slice(prod_idx)
    n_frames = tr.n_frames
    frame_counts[rep] = n_frames
    time_ns = (np.array(prod_idx) - SKIP) * FRAME_DT_NS
    print(f"rep{rep}: total {tr_all.n_frames} frames, production {n_frames}")

    # 1) RMSD ---------------------------------------------------------------
    ref = tr.slice(0)
    bb_rmsd = md.rmsd(tr, ref, atom_indices=backbone, parallel=False)
    tr.superpose(ref, atom_indices=backbone)

    def rmsd_noalign(coords, refc):
        d = coords - refc
        return np.sqrt((d * d).sum(axis=(1, 2)) / coords.shape[1])

    lig_rmsd = rmsd_noalign(tr.xyz[:, lig_heavy, :],
                            tr.xyz[0, lig_heavy, :])
    pd.DataFrame({
        "frame": prod_idx,
        "time_ns": time_ns,
        "protein_backbone_RMSD_nm": bb_rmsd,
        "ligand_heavy_RMSD_nm": lig_rmsd,
    }).to_csv(os.path.join(OUT, f"wp7_rep{rep}_rmsd.csv"),
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
    }).to_csv(os.path.join(OUT, f"wp7_rep{rep}_rmsf.csv"),
              index=False, encoding="utf-8-sig")
    rmsf_store[rep] = ca_rmsf

    # 3) Hydrogen bonds ------------------------------------------------------
    solute = tr.atom_slice(prot_atoms + LIG)
    stop = solute.topology
    lig_set = set(range(N_PROT, N_PROT + len(LIG)))   # sliced topology indices
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
    }).to_csv(os.path.join(OUT, f"wp7_rep{rep}_hb.csv"),
              index=False, encoding="utf-8-sig")
    occ_df = pd.DataFrame(
        [(rn, nm, role, c, round(c / n_frames, 3))
         for (rn, nm, role), c in occ_counter.items()],
        columns=["resSeq", "residue", "ligand_role",
                 "n_frames", "occupancy"])
    occ_df.to_csv(os.path.join(OUT, f"wp7_rep{rep}_hb_occupancy.csv"),
                  index=False, encoding="utf-8-sig")
    occ_store[rep] = {k: c / n_frames for k, c in occ_counter.items()}

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
    }).to_csv(os.path.join(OUT, f"wp7_rep{rep}_distance.csv"),
              index=False, encoding="utf-8-sig")

    per_frame[rep] = {
        "time_ns": time_ns,
        "bb_rmsd": bb_rmsd,
        "lig_rmsd": lig_rmsd,
        "hb": hb_count,
        "min_dist": min_dist,
        "com_dist": com_dist,
    }
    del tr_all, tr, solute
    print(f"rep{rep}: done and released")

# ---------------------------------------------------------------------------
# Cross-replicate summary
# ---------------------------------------------------------------------------
summary_rows = []
for rep in range(N_REPS):
    pf = per_frame[rep]
    half = len(pf["bb_rmsd"]) // 2
    row = {"replicate": f"rep{rep}", "n_production_frames": frame_counts[rep]}
    for metric, label in [
        ("bb_rmsd", "backbone_RMSD_nm"),
        ("lig_rmsd", "ligand_RMSD_nm"),
        ("hb", "HB_count"),
        ("min_dist", "min_distance_nm"),
        ("com_dist", "COM_distance_nm"),
    ]:
        x = pf[metric]
        row[f"mean_{label}"] = float(x.mean())
        row[f"SD_{label}"] = float(x.std(ddof=1))
        row[f"first_half_mean_{label}"] = float(x[:half].mean())
        row[f"second_half_mean_{label}"] = float(x[half:].mean())
    summary_rows.append(row)
summary_df = pd.DataFrame(summary_rows)

# pooled mean and between-replicate SD for each headline metric
pooled = {}
for metric, label in [
    ("bb_rmsd", "backbone_RMSD_nm"),
    ("lig_rmsd", "ligand_RMSD_nm"),
    ("hb", "HB_count"),
    ("min_dist", "min_distance_nm"),
    ("com_dist", "COM_distance_nm"),
]:
    means = np.array([per_frame[r][metric].mean() for r in range(N_REPS)])
    pooled[f"mean_{label}"] = float(means.mean())
    pooled[f"between_rep_SD_{label}"] = float(means.std(ddof=1))
    pooled[f"between_rep_CV_{label}"] = (
        float(means.std(ddof=1) / means.mean()) if means.mean() != 0 else None)
summary_df.to_csv(os.path.join(OUT, "wp7_cross_replicate_summary.csv"),
                  index=False, encoding="utf-8-sig")

# RMSF agreement across replicates
rmsf_mat = np.vstack([rmsf_store[r] for r in range(N_REPS)])
rmsf_cross = pd.DataFrame({
    "resSeq": [k[0] for k in res_order],
    "residue": [k[1] for k in res_order],
    "rep0_CA_RMSF_nm": rmsf_mat[0],
    "rep1_CA_RMSF_nm": rmsf_mat[1],
    "rep2_CA_RMSF_nm": rmsf_mat[2],
    "mean_CA_RMSF_nm": rmsf_mat.mean(axis=0),
    "between_rep_SD_CA_RMSF_nm": rmsf_mat.std(axis=0, ddof=1),
})
rmsf_cross.to_csv(os.path.join(OUT, "wp7_rmsf_cross_replicate.csv"),
                  index=False, encoding="utf-8-sig")
corr = np.corrcoef(rmsf_mat)
pairwise_corr = {
    "rep0_rep1": float(corr[0, 1]),
    "rep0_rep2": float(corr[0, 2]),
    "rep1_rep2": float(corr[1, 2]),
}

# hydrogen bond occupancy comparison on the union of contacts
all_keys = sorted(set().union(*[set(occ_store[r].keys())
                                for r in range(N_REPS)]))
hb_rows = []
for key in all_keys:
    rn, nm, role = key
    hb_rows.append({
        "resSeq": rn,
        "residue": nm,
        "ligand_role": role,
        "rep0_occupancy": occ_store[0].get(key, 0.0),
        "rep1_occupancy": occ_store[1].get(key, 0.0),
        "rep2_occupancy": occ_store[2].get(key, 0.0),
        "mean_occupancy": float(np.mean(
            [occ_store[r].get(key, 0.0) for r in range(N_REPS)])),
    })
pd.DataFrame(hb_rows).to_csv(
    os.path.join(OUT, "wp7_hb_occupancy_compare.csv"),
    index=False, encoding="utf-8-sig")

comparison = {
    "design": "three independent 20 ns replicates, velocity seeds "
              "20262117, 20263117, 20264117; 4 fs timestep, OpenCL",
    "production_frames_each": int(frame_counts[0]),
    "pooled": pooled,
    "per_replicate_means": {
        f"rep{r}": {
            "backbone_RMSD_nm": float(per_frame[r]["bb_rmsd"].mean()),
            "ligand_RMSD_nm": float(per_frame[r]["lig_rmsd"].mean()),
            "HB_count": float(per_frame[r]["hb"].mean()),
            "min_distance_nm": float(per_frame[r]["min_dist"].mean()),
            "COM_distance_nm": float(per_frame[r]["com_dist"].mean()),
        } for r in range(N_REPS)
    },
    "RMSF_pairwise_pearson": pairwise_corr,
    "n_distinct_ligand_protein_HB_contacts": len(all_keys),
    "consistently_occupied_HB_occupancy_ge_0_2": [
        {"resSeq": int(r["resSeq"]), "residue": r["residue"],
         "ligand_role": r["ligand_role"],
         "mean_occupancy": r["mean_occupancy"]}
        for r in hb_rows if r["mean_occupancy"] >= 0.2],
}
with open(os.path.join(OUT, "wp7_replicate_comparison.json"), "w",
          encoding="utf-8") as fh:
    json.dump(comparison, fh, ensure_ascii=False, indent=2)

print("\n=== pooled cross-replicate metrics ===")
for k, v in pooled.items():
    print(f"{k}: {v:.4f}")
print("RMSF pairwise Pearson:", pairwise_corr)
print("distinct HB contacts:", len(all_keys))
print("wrote", OUT)
