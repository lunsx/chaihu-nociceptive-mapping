# -*- coding: utf-8 -*-
"""Supplementary figure: suspension-type sensitivity of the healthy reference.

Shows that the low healthy macrophage expression of IL1B, TNF, and PTGS2 holds
within both whole-cell and nuclear suspensions, so the result is not an artifact
of nuclear preparation. Reads the table produced by external_healthy_reference.py
and writes an English figure to the repository figures directory.

Run:
  python scripts/07_singlecell_spatial/fig_suspension_sensitivity.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

IN = PATHS["processed_small"] / "external_healthy_suspension_sensitivity.csv"
OUT = PATHS["figures"] / "fig_suspension_sensitivity.png"

GENES = ["IL1B", "TNF", "PTGS2"]
MASLD_MEDIAN = {"IL1B": 80.0, "TNF": 52.58, "PTGS2": 44.01}


def main():
    df = pd.read_csv(IN)
    combined = df[df["lineage"] == "InfMac+NonInfMac"].set_index("suspension")

    x = np.arange(len(GENES))
    width = 0.25
    colors = {"cell": "#4C78A8", "nucleus": "#F58518", "ALL": "#54A24B"}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: healthy percentages by suspension, small scale so differences show.
    ax = axes[0]
    for i, sus in enumerate(["cell", "nucleus", "ALL"]):
        vals = [combined.loc[sus, g] for g in GENES]
        bars = ax.bar(x + (i - 1) * width, vals, width, label=sus, color=colors[sus])
        ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES)
    ax.set_ylabel("Healthy cells expressing the gene (%)")
    ax.set_ylim(0, 9)
    ax.set_title("A. Healthy macrophages by suspension")
    ax.legend(frameon=False, fontsize=9)

    # Panel B: healthy ALL vs MASLD median on the full scale.
    ax = axes[1]
    healthy = [combined.loc["ALL", g] for g in GENES]
    disease = [MASLD_MEDIAN[g] for g in GENES]
    b1 = ax.bar(x - width / 2, healthy, width, label="Healthy (all)", color="#54A24B")
    b2 = ax.bar(x + width / 2, disease, width, label="MASLD median", color="#B279A2")
    ax.bar_label(b1, fmt="%.2f", padding=2, fontsize=8)
    ax.bar_label(b2, fmt="%.1f", padding=2, fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES)
    ax.set_ylabel("Cells expressing the gene (%)")
    ax.set_ylim(0, 90)
    ax.set_title("B. Healthy vs MASLD")
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
