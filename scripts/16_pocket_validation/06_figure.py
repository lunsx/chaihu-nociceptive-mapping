# -*- coding: utf-8 -*-
"""Summary figure for the independent pocket-validation analysis.

Panel A: P2Rank pocket probabilities by rank, with the canonical pocket
(pocket6) highlighted. Panel B: best-pose centroid distance to the canonical
site for targeted (pocket-box) versus whole-protein blind docking. Distinct
markers and a colourblind-safe palette keep the panels readable in grayscale.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
FIG = PATHS["figures"] / "pocket_validation.png"
CANON_POCKET = 6

plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
})


def main():
    pockets = pd.read_csv(OUT_DIR / "pocket_match_summary.csv")
    poses = pd.read_csv(OUT_DIR / "pose_rmsd_summary.csv")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))

    colors = ["#d62728" if r == CANON_POCKET else "#1f77b4"
              for r in pockets["rank"]]
    axA.bar(pockets["rank"], pockets["probability"], color=colors,
            edgecolor="black", linewidth=0.6)
    axA.set_xlabel("P2Rank pocket rank")
    axA.set_ylabel("Pocket probability")
    axA.set_title("A. Independent pocket detection")
    axA.annotate("canonical site\n(pocket6)", xy=(6, 0.323),
                 xytext=(8.5, 0.55),
                 arrowprops=dict(arrowstyle="->", color="#d62728"),
                 color="#d62728", fontsize=9)

    groups = [
        ("Targeted\nVina", "targeted_vina", "o"),
        ("Targeted\nVinardo", "targeted_vinardo", "s"),
        ("Blind\nVina", "blind:blind_vina_2026011901", "^"),
        ("Blind\nVinardo", "blind:blind_vinardo_2026011901", "D"),
    ]
    xs, labels = [], []
    for i, (label, src, marker) in enumerate(groups, 1):
        sub = poses[poses["source"] == src]
        best = sub.loc[sub["affinity"].idxmin()]
        col = "#2ca02c" if i <= 2 else "#7f7f7f"
        axB.scatter(i, best["centroid_distance"], marker=marker, s=110,
                    color=col, edgecolor="black", zorder=3)
        xs.append(i); labels.append(label)
    axB.axhline(2.0, color="#d62728", linestyle=":", linewidth=1.2,
                label="2 Å success threshold")
    axB.set_xticks(xs); axB.set_xticklabels(labels)
    axB.set_ylabel("Best-pose centroid distance (Å)")
    axB.set_title("B. Targeted vs whole-protein blind docking")
    axB.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=300, bbox_inches="tight")
    print("saved", FIG)


if __name__ == "__main__":
    main()
