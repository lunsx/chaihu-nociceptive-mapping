# -*- coding: utf-8 -*-
"""Figure 10 (composite): in silico CellOracle perturbation + L1000 reversal drugs.

Rebuilds the two-panel submission figure from the source tables and writes both
a 300 DPI PNG and a vector PDF. The annotation helper suppresses negative zero.

Run:
  python scripts/08_perturbation/fig10_composite.py
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

ROOT = PATHS["project_root"]
KO_CSV = ROOT / "knockout" / "虚拟敲除过表达_下游变化.csv"
L1000_CSV = ROOT / "lincs" / "Enrichr_LINCS_L1000_Chem_Pert_Consensus_Sigs.csv"
OUT_PNG = ROOT / "投稿图集" / "Fig9_因果扰动.png"
OUT_PDF = ROOT / "投稿图集" / "Fig9_因果扰动.pdf"

GENE_ORDER = ["NFKBIA", "IL1B", "HMOX1", "TNF", "ICAM1", "CXCL2", "SERPINE1",
              "CCL2", "PTGS2", "IL6", "SPP1", "CALCRL", "PTGER3", "PTGES", "IL1A"]
PERT_ORDER = ["KO_RELA", "KO_NFKB1", "KO_JUN", "KO_FOS", "KO_STAT3", "OE_RELA"]


def safe_text(v, decimals=2):
    if np.isnan(v):
        return ""
    r = float(np.round(v, decimals))
    if r == 0.0:
        r = 0.0
    return f"{r:.{decimals}f}"


def main():
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 10

    ko = pd.read_csv(KO_CSV)
    mat = (ko.pivot_table(index="target_gene", columns="perturbation",
                          values="mean_delta")
             .reindex(index=GENE_ORDER, columns=PERT_ORDER))

    l1000 = pd.read_csv(L1000_CSV)
    down = l1000[l1000.term.str.endswith("Down")].copy()
    top = down.sort_values("adjusted_pvalue").head(15).copy()
    top["drug"] = top.term.str.replace(" Down", "")
    top = top.iloc[::-1]

    fig, axes = plt.subplots(
        1, 2, figsize=(16, 7.4),
        gridspec_kw={"width_ratios": [1.15, 1.0]},
    )
    fig.suptitle("Causal perturbation evidence: in silico knockout and L1000",
                 fontsize=16, y=0.99)

    ax = axes[0]
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r",
                   vmin=-0.8, vmax=0.8)
    ax.set_xticks(np.arange(len(PERT_ORDER)))
    ax.set_xticklabels(PERT_ORDER, fontsize=10)
    ax.set_yticks(np.arange(len(GENE_ORDER)))
    ax.set_yticklabels(GENE_ORDER, fontsize=10)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, safe_text(v, 2), ha="center", va="center",
                        color="white" if abs(v) >= 0.55 else "black",
                        fontsize=10)
    ax.set_title("")
    ax.text(0.0, 1.15, "A. In silico knockout/overexpression",
            transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="bottom")
    ax.text(0.0, 1.05, "CellOracle in silico perturbation: predicted expression shift",
            transform=ax.transAxes, fontsize=11, va="bottom")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("mean delta expression")

    ax = axes[1]
    ax.barh(top.drug, -np.log10(top.adjusted_pvalue), color="#3b6ea5")
    ax.set_xlabel("-log10(adjusted p-value)")
    ax.set_title("")
    ax.text(0.0, 1.15, "B. L1000 reversal drugs",
            transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="bottom")
    ax.text(0.0, 1.05,
            "L1000 consensus: drugs that reverse the inflammation-up signature",
            transform=ax.transAxes, fontsize=11, va="bottom")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print("saved", OUT_PNG)
    print("saved", OUT_PDF)


if __name__ == "__main__":
    main()
