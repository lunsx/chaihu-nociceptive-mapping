# -*- coding: utf-8 -*-
"""Figure 10: external healthy reference vs MASLD single-cell percentages.

Reads the table produced by c2_healthy_comparison.py and writes an English
figure to the repository figures directory. Paths are resolved from config.

Run:
  python scripts/07_singlecell_spatial/fig10_external_healthy.py
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

CSV = PATHS["processed_small"] / "external_healthy_reference_comparison.csv"
OUT = PATHS["figures"] / "fig_external_healthy_reference.png"

KEY_GENES = ["IL1B", "TNF", "PTGS2", "IL6", "FOS", "ICAM1"]


def safe_text(v, decimals=0):
    """Format a cell value without emitting a negative zero label."""
    if np.isnan(v):
        return ""
    r = float(np.round(v, decimals))
    if r == 0.0:
        r = 0.0
    return f"{r:.{decimals}f}"


def main():
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 9

    df = pd.read_csv(CSV)
    d = df[df["gene"].isin(KEY_GENES)].copy()

    fig, axes = plt.subplots(
        1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.15, 1]}
    )

    ax = axes[0]
    sub = d[d["cell_type"] == "Macrophage/Kupffer"].set_index("gene").loc[KEY_GENES]
    x = np.arange(len(KEY_GENES))
    w = 0.27
    ax.bar(x - w, sub["healthy_pct"], w,
           label="External healthy (Aizarani, 4 donors)", color="#4C78A8")
    ax.bar(x, sub["MASLD_donor_median_pct"], w,
           label="MASLD donor median (19 donors)", color="#E45756")
    ax.bar(x + w, sub["single_GSE344087_HC_pct"], w,
           label="Single GSE344087 HC", color="#9D9D9D")
    ax.set_xticks(x)
    ax.set_xticklabels(KEY_GENES)
    ax.set_ylabel("Expressing cells (%)")
    ax.set_title("Macrophage / Kupffer")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    piv_m = d.pivot_table(index="gene", columns="cell_type",
                          values="MASLD_donor_median_pct")
    piv_h = d.pivot_table(index="gene", columns="cell_type", values="healthy_pct")
    diff = (piv_m - piv_h).loc[KEY_GENES]
    cols = [c for c in ["Macrophage/Kupffer", "Endothelial", "Hepatic stellate",
                        "T/NK", "B", "Cholangiocyte", "Hepatocyte"]
            if c in diff.columns]
    diff = diff[cols]
    im = ax.imshow(diff.values, aspect="auto", cmap="Reds", vmin=0, vmax=80)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(KEY_GENES)))
    ax.set_yticklabels(KEY_GENES)
    for i in range(diff.shape[0]):
        for j in range(diff.shape[1]):
            v = diff.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, safe_text(v, 0), ha="center", va="center",
                        color="white" if v > 45 else "black", fontsize=8)
    ax.set_title("MASLD minus healthy (percentage points)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="percentage points")

    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print("saved", OUT)


if __name__ == "__main__":
    main()
