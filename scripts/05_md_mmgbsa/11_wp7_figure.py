# -*- coding: utf-8 -*-
"""Cross-replicate figure for WP7 (three independent 20 ns MD runs).

Builds a six-panel summary from the small committed WP7 data files:

  A. protein backbone RMSD over time, three replicates;
  B. ligand heavy-atom RMSD over time;
  C. minimum ligand-protein distance over time (pocket retention);
  D. per-residue CA RMSF, cross-replicate mean with between-rep SD band;
  E. single-trajectory MM-GBSA cumulative running mean;
  F. MM-GBSA 4 ns block means with the pooled mean reference.

The figure is descriptive and reads the same metric definitions used by
09_wp7_structure.py and 10_wp7_mmgbsa.py. Only numpy and matplotlib are
required; CSV parsing uses the standard library.

Run:  python scripts/05_md_mmgbsa/11_wp7_figure.py
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
WP7 = ROOT / "data" / "processed_small" / "wp7"
FIGDIR = ROOT / "figures"

REPS = [0, 1, 2]
COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}
LABELS = {0: "rep0", 1: "rep1", 2: "rep2"}

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["font.size"] = 10
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["lines.linewidth"] = 1.3


def read_table(name):
    with open(WP7 / name, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fcol(rows, key):
    return np.array([float(r[key]) for r in rows])


def icol(rows, key):
    return np.array([int(r[key]) for r in rows])


def panel_label(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")


rmsd = {r: read_table(f"wp7_rep{r}_rmsd.csv") for r in REPS}
dist = {r: read_table(f"wp7_rep{r}_distance.csv") for r in REPS}
rmsf_x = read_table("wp7_rmsf_cross_replicate.csv")
per_frame = read_table("wp7_mmgbsa_per_frame.csv")
blocks = read_table("wp7_mmgbsa_block_convergence.csv")

fig, axes = plt.subplots(2, 3, figsize=(15, 8.2))

# A. backbone RMSD
ax = axes[0, 0]
for r in REPS:
    ax.plot(fcol(rmsd[r], "time_ns"), fcol(rmsd[r], "protein_backbone_RMSD_nm"),
            color=COLORS[r], label=LABELS[r])
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Backbone RMSD (nm)")
ax.set_title("Protein backbone stability")
ax.legend(frameon=False, fontsize=9)
panel_label(ax, "A")

# B. ligand heavy-atom RMSD
ax = axes[0, 1]
for r in REPS:
    ax.plot(fcol(rmsd[r], "time_ns"), fcol(rmsd[r], "ligand_heavy_RMSD_nm"),
            color=COLORS[r], label=LABELS[r])
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Ligand heavy-atom RMSD (nm)")
ax.set_title("Ligand pose adjustment")
ax.legend(frameon=False, fontsize=9)
panel_label(ax, "B")

# C. minimum ligand-protein distance
ax = axes[0, 2]
for r in REPS:
    ax.plot(fcol(dist[r], "time_ns"),
            fcol(dist[r], "min_lig_prot_distance_nm"),
            color=COLORS[r], label=LABELS[r])
ax.axhline(0.4, color="gray", linestyle="--", linewidth=1.0)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Minimum distance (nm)")
ax.set_title("Ligand pocket retention")
ax.set_ylim(0.2, 0.45)
ax.legend(frameon=False, fontsize=9)
panel_label(ax, "C")

# D. per-residue RMSF with between-rep SD band
ax = axes[1, 0]
x = icol(rmsf_x, "resSeq")
mean = fcol(rmsf_x, "mean_CA_RMSF_nm")
sd = fcol(rmsf_x, "between_rep_SD_CA_RMSF_nm")
ax.fill_between(x, mean - sd, mean + sd, color="#555555", alpha=0.25,
                label="between-rep SD")
ax.plot(x, mean, color="#222222", label="cross-rep mean")
ax.set_xlabel("Residue sequence number")
ax.set_ylabel("CA RMSF (nm)")
ax.set_title("Residue flexibility")
ax.legend(frameon=False, fontsize=9)
panel_label(ax, "D")

# E. cumulative running mean of MM-GBSA
ax = axes[1, 1]
for r in REPS:
    sub = sorted((row for row in per_frame if int(row["replicate"]) == r),
                 key=lambda row: float(row["time_ns"]))
    vals = fcol(sub, "dG_total_kcal")
    cum = np.cumsum(vals) / np.arange(1, len(vals) + 1)
    ax.plot(fcol(sub, "time_ns"), cum, color=COLORS[r], label=LABELS[r])
pooled = fcol(per_frame, "dG_total_kcal").mean()
ax.axhline(pooled, color="black", linestyle="--", linewidth=1.0,
           label="pooled mean")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Cumulative mean dG (kcal/mol)")
ax.set_title("MM-GBSA running convergence")
ax.legend(frameon=True, framealpha=0.9, edgecolor="none",
          fontsize=8, loc="lower right")
panel_label(ax, "E")

# F. block means
ax = axes[1, 2]
for r in REPS:
    sub = [row for row in blocks if int(row["replicate"]) == r]
    m = (fcol(sub, "block_start_ns") + fcol(sub, "block_end_ns")) / 2.0
    ax.scatter(m, fcol(sub, "block_mean_dG_kcal_mol"), color=COLORS[r],
               s=28, label=LABELS[r], zorder=3)
ax.axhline(pooled, color="black", linestyle="--", linewidth=1.0,
           label="pooled mean")
ax.set_xlabel("Block midpoint (ns)")
ax.set_ylabel("4 ns block mean dG (kcal/mol)")
ax.set_title("Block-level free-energy scatter")
ax.legend(frameon=False, fontsize=9)
panel_label(ax, "F")

for ax in axes.flat:
    ax.grid(True, alpha=0.25, linewidth=0.6)

fig.tight_layout()
FIGDIR.mkdir(parents=True, exist_ok=True)
png = FIGDIR / "wp7_cross_replicate.png"
pdf = FIGDIR / "wp7_cross_replicate.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print("saved", png)
print("saved", pdf)
