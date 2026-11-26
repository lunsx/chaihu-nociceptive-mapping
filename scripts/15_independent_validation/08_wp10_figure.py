# -*- coding: utf-8 -*-
"""Multi-pose validation figure for WP10 (ten 20 ns arms).

Builds a six-panel summary from the small committed WP10 data files:

  A. protein backbone RMSD over time (all arms overlap; protein is stable);
  B. protein-aligned ligand centre-of-mass drift (site retention split);
  C. minimum ligand-protein distance (detached arms separate into solvent);
  D. per-arm MM-GBSA total dG, colour-coded by ligand outcome;
  E. ligand heavy-atom RMSD over time (pose adjustment);
  F. mean ligand-protein hydrogen-bond count per arm.

The figure is descriptive and reads the same metric definitions used by
06_wp10_structure.py and 07_wp10_mmgbsa.py. Only numpy and matplotlib are
required; CSV parsing uses the standard library.

Run:  python scripts/15_independent_validation/08_wp10_figure.py
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
WP10 = ROOT / "data" / "processed_small" / "wp10"
FIGDIR = ROOT / "figures"

ARMS = list(range(10))
# tab10 palette; arm0 (canonical) is restyled as a thick black reference
LINE_COLORS = plt.get_cmap("tab10").colors
LINE_STYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--"]

# Ligand outcome classification from the structural analysis
STATUS = {
    0: "retained", 1: "retained", 2: "retained", 3: "retained",
    4: "retained", 5: "re-associated", 6: "detached", 7: "detached",
    8: "relocated", 9: "retained",
}
STATUS_COLOR = {
    "retained": "#2ca02c",
    "re-associated": "#9467bd",
    "relocated": "#ff7f0e",
    "detached": "#d62728",
}

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["font.size"] = 10
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["lines.linewidth"] = 1.2


def read_table(name):
    with open(WP10 / name, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fcol(rows, key):
    return np.array([float(r[key]) for r in rows])


def panel_label(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")


rmsd = {a: read_table(f"wp10_arm{a}_rmsd.csv") for a in ARMS}
dist = {a: read_table(f"wp10_arm{a}_distance.csv") for a in ARMS}
summary = read_table("wp10_cross_arm_summary.csv")
mmgbsa_summary = read_table("wp10_mmgbsa_summary.csv")

fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))


def plot_arms(ax, source, key, ylabel, title):
    for a in ARMS:
        rows = source[a]
        if a == 0:
            ax.plot(fcol(rows, "time_ns"), fcol(rows, key), color="black",
                    linewidth=2.0, linestyle="-", label="arm0 canonical")
        else:
            ax.plot(fcol(rows, "time_ns"), fcol(rows, key),
                    color=LINE_COLORS[a], linestyle=LINE_STYLES[a],
                    label=f"arm{a}")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)


# A. backbone RMSD
ax = axes[0, 0]
plot_arms(ax, rmsd, "protein_backbone_RMSD_nm",
          "Backbone RMSD (nm)", "Protein backbone stability")
ax.legend(frameon=False, fontsize=7, ncol=2)
panel_label(ax, "A")

# B. ligand COM drift
ax = axes[0, 1]
plot_arms(ax, rmsd, "ligand_COM_drift_nm",
          "Ligand COM drift (nm)", "Ligand site retention")
panel_label(ax, "B")

# C. minimum distance
ax = axes[0, 2]
plot_arms(ax, dist, "min_lig_prot_distance_nm",
          "Minimum distance (nm)", "Ligand contact / detachment")
ax.axhline(0.4, color="gray", linestyle="--", linewidth=1.0)
panel_label(ax, "C")

# D. per-arm MM-GBSA dG bar, colour-coded by outcome
ax = axes[1, 0]
tot = {int(r["arm"]): float(r["mean_kcal_mol"])
       for r in mmgbsa_summary if r["component"] == "total"}
bar_colors = [STATUS_COLOR[STATUS[a]] for a in ARMS]
bars = ax.bar([f"arm{a}" for a in ARMS], [tot[a] for a in ARMS],
              color=bar_colors, edgecolor="black", linewidth=0.6)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_ylabel("MM-GBSA dG (kcal/mol)")
ax.set_title("Binding free energy by starting pose")
ax.tick_params(axis="x", labelrotation=45, labelsize=8)
handles = [plt.Line2D([0], [0], marker="s", color="w",
                      markerfacecolor=c, markersize=9, label=k)
           for k, c in STATUS_COLOR.items()]
ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower left")
panel_label(ax, "D")

# E. ligand heavy RMSD
ax = axes[1, 1]
plot_arms(ax, rmsd, "ligand_heavy_RMSD_nm",
          "Ligand heavy-atom RMSD (nm)", "Ligand pose adjustment")
panel_label(ax, "E")

# F. mean HB count per arm
ax = axes[1, 2]
hb_means = {int(r["arm"]): float(r["mean_HB_count"]) for r in summary}
ax.bar([f"arm{a}" for a in ARMS], [hb_means[a] for a in ARMS],
       color=bar_colors, edgecolor="black", linewidth=0.6)
ax.set_ylabel("Mean hydrogen-bond count")
ax.set_title("Ligand-protein hydrogen bonding")
ax.tick_params(axis="x", labelrotation=45, labelsize=8)
ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper right")
panel_label(ax, "F")

for ax in axes.flat:
    ax.grid(True, alpha=0.25, linewidth=0.6)

fig.tight_layout()
FIGDIR.mkdir(parents=True, exist_ok=True)
png = FIGDIR / "wp10_multipose.png"
pdf = FIGDIR / "wp10_multipose.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print("saved", png)
print("saved", pdf)
