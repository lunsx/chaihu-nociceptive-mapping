# -*- coding: utf-8 -*-
"""
07  Conformational / pose clustering for the 5KIR-luteolin trajectory.

Fills the round-6 open item "pocket volume / conformational clustering":
pocket volume is covered by 16_pocket_volume.py; this script does the clustering.

Two independent clusterings on the production frames:
  1. Ligand pose: pairwise heavy-atom RMSD after backbone superposition,
     agglomerative clustering with a 0.15 nm distance threshold.
  2. Protein backbone: pairwise backbone RMSD, agglomerative clustering
     with a 0.10 nm threshold.

Outputs (data/processed_small):
  md_pose_clusters.csv      per-frame cluster labels (ligand and protein)
  md_pose_cluster_summary.json  cluster sizes, medoid frames, within-cluster RMSD
Figures:
  figures/md_pose_clustering.png
"""
import json
import os
import pickle

import numpy as np
import pandas as pd
import mdtraj as md
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

MD = os.environ.get("CHAIHU_MD_DIR", r"C:\dockwork\md")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "data", "processed_small")
FIGDIR = os.path.join(REPO, "figures")

m = pickle.load(open(os.path.join(MD, "modeller.pkl"), "rb"))
top = md.Topology.from_openmm(m["topology"])
lig_atoms = list(m["lig_idx"])
n_prot = int(m["n_prot"])
backbone = top.select("backbone and chainid 0")


def is_heavy(i):
    return top.atom(i).element.symbol != "H"


lig_heavy = [i for i in lig_atoms if is_heavy(i)]

tr_all = md.load(os.path.join(MD, "traj.dcd"), top=top)
SKIP = 6
prod_idx = list(range(SKIP, tr_all.n_frames))
tr = tr_all.slice(prod_idx)
t_ns = (np.arange(tr_all.n_frames) * 0.04)[SKIP:]
print("production frames:", tr.n_frames, "| ligand heavy atoms:", len(lig_heavy),
      "| backbone atoms:", len(backbone))

# Superpose on backbone using the first production frame as reference.
ref = tr.slice(0)
tr.superpose(ref, atom_indices=backbone)


def pairwise_rmsd(atom_idx):
    """Full pairwise RMSD matrix (nm) over the given atom indices."""
    x = tr.xyz[:, atom_idx, :]
    n = tr.n_frames
    mat = np.zeros((n, n))
    for i in range(n):
        d = x[i + 1:] - x[i]
        rms = np.sqrt((d * d).sum(axis=(1, 2)) / len(atom_idx))
        mat[i, i + 1:] = rms
        mat[i + 1:, i] = rms
    return mat


def cluster_from_matrix(mat, threshold_nm):
    """Agglomerative (single-linkage) clustering at a distance threshold."""
    condensed = squareform(mat, checks=False)
    z = linkage(condensed, method="average")
    labels = fcluster(z, t=threshold_nm, criterion="distance")
    return labels


def summarize(mat, labels):
    summary = []
    for cid in sorted(set(labels)):
        members = np.where(labels == cid)[0]
        sub = mat[np.ix_(members, members)]
        # medoid: member with minimum mean distance to the rest of the cluster
        medoid_local = members[np.argmin(sub.mean(axis=1))]
        summary.append({
            "cluster": int(cid),
            "n_frames": int(len(members)),
            "frame_fraction": round(len(members) / len(labels), 3),
            "medoid_frame": int(prod_idx[medoid_local]),
            "medoid_time_ns": round(float(t_ns[medoid_local]), 3),
            "mean_within_cluster_RMSD_nm": round(float(sub.mean()), 4),
            "max_within_cluster_RMSD_nm": round(float(sub.max()), 4),
        })
    return summary


lig_mat = pairwise_rmsd(lig_heavy)
bb_mat = pairwise_rmsd(list(backbone))

LIG_THR = 0.15
# Backbone RMSD is tightly distributed (mean 0.148 nm, p75 0.160, max 0.200);
# a 0.16 nm threshold sits near the upper quartile so the dominant basin is
# retained rather than fragmented by a threshold below the normal fluctuation.
BB_THR = 0.16
lig_labels = cluster_from_matrix(lig_mat, LIG_THR)
bb_labels = cluster_from_matrix(bb_mat, BB_THR)

# Relabel clusters in descending size so cluster 1 is the largest.
def relabel_by_size(labels):
    sizes = pd.Series(labels).value_counts()
    order = {old: new for new, old in enumerate(sizes.index, start=1)}
    return np.array([order[x] for x in labels])


lig_labels = relabel_by_size(lig_labels)
bb_labels = relabel_by_size(bb_labels)

lig_summary = summarize(lig_mat, lig_labels)
bb_summary = summarize(bb_mat, bb_labels)

per_frame = pd.DataFrame({
    "frame": prod_idx,
    "time_ns": np.round(t_ns, 3),
    "ligand_pose_cluster": lig_labels,
    "backbone_cluster": bb_labels,
})
per_frame.to_csv(os.path.join(OUT, "md_pose_clusters.csv"),
                 index=False, encoding="utf-8-sig")

summary = {
    "n_production_frames": int(tr.n_frames),
    "ligand_pose": {
        "atom_indices": "ligand heavy atoms",
        "n_atoms": int(len(lig_heavy)),
        "distance_threshold_nm": LIG_THR,
        "linkage": "average",
        "n_clusters": int(len(lig_summary)),
        "clusters": lig_summary,
    },
    "protein_backbone": {
        "atom_indices": "backbone and chainid 0",
        "n_atoms": int(len(backbone)),
        "distance_threshold_nm": BB_THR,
        "linkage": "average",
        "n_clusters": int(len(bb_summary)),
        "clusters": bb_summary,
    },
}
with open(os.path.join(OUT, "md_pose_cluster_summary.json"), "w",
          encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=1)

print("ligand clusters:", [(c["cluster"], c["n_frames"]) for c in lig_summary])
print("backbone clusters:", [(c["cluster"], c["n_frames"]) for c in bb_summary])

# ---------- Figure ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
axs[0].scatter(t_ns, lig_labels, c=lig_labels, cmap="viridis", s=28, zorder=3)
axs[0].set_ylabel("Ligand pose cluster")
axs[0].set_title("A. Ligand pose cluster over time (threshold 0.15 nm)",
                  loc="left", fontweight="bold")
axs[0].set_yticks(range(1, len(lig_summary) + 1))
axs[1].scatter(t_ns, bb_labels, c=bb_labels, cmap="plasma", s=28, zorder=3)
axs[1].set_ylabel("Backbone cluster")
axs[1].set_xlabel("Time (ns)")
axs[1].set_title("B. Protein backbone cluster over time (threshold 0.16 nm)",
                  loc="left", fontweight="bold")
axs[1].set_yticks(range(1, len(bb_summary) + 1))
fig.tight_layout()
figpath = os.path.join(FIGDIR, "md_pose_clustering.png")
fig.savefig(figpath, dpi=200, bbox_inches="tight", facecolor="white")
print("saved", figpath)
