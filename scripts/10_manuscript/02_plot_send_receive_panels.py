# -*- coding: utf-8 -*-
"""Panels A-C of the signal send / receive figure, built from the result tables.

Values in panel C are read from the tables, not typed in: for each ligand the
figure prints the cell type with the highest percentage and that percentage, so
the picture cannot drift away from the data behind it.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

ROOT = Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析")
OUT = ROOT / "补强分析_20261114"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

LABELS = [
    ("CD8-positive, alpha-beta T cell", "CD8 T"),
    ("CD4-positive, alpha-beta T cell", "CD4 T"),
    ("conventional dendritic cell", "cDC"),
    ("plasmacytoid dendritic cell", "pDC"),
    ("endothelial cell of pericentral hepatic sinusoid", "LSEC-pc"),
    ("endothelial cell of periportal hepatic sinusoid", "LSEC-pp"),
    ("endothelial cell of artery", "Artery EC"),
    ("vein endothelial cell", "Vein EC"),
    ("periportal region hepatocyte", "Hep-periportal"),
    ("centrilobular region hepatocyte", "Hep-centrilobular"),
    ("midzonal region hepatocyte", "Hep-midzonal"),
    ("intrahepatic cholangiocyte", "Cholangiocyte"),
    ("hepatic stellate cell", "HSC"),
    ("hepatic pit cell", "Pit cell"),
    ("natural killer cell", "NK"),
    ("Kupffer cell", "Kupffer"),
    ("mast cell", "Mast"),
    ("mature B cell", "B"),
    ("plasma cell", "Plasma"),
]


def shorten(x):
    x = str(x)
    for long, short in LABELS:
        if x == long:
            return short
    return x[:16]


def panel_a(ax):
    d = pd.read_csv(ROOT / "human_sc" / "人源单细胞_细胞类型定位.csv")
    genes = ["IL1B", "TNF", "PTGS2", "IL6", "PTGES", "CALCRL", "IL1A"]
    top = (d.groupby("cell_type")["n_cells"].max().sort_values(ascending=False)
           .head(12).index.tolist())
    m = (d[d.cell_type.isin(top) & d.gene.isin(genes)]
         .pivot_table(index="cell_type", columns="gene", values="pct_expressing"))
    m = m.reindex(index=top, columns=genes)
    vmax = float(np.nanmax(m.values))
    im = ax.imshow(m.values, cmap="Reds", vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=6.5, style="italic")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([shorten(t) for t in top], fontsize=6.2)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if np.isfinite(v) and v >= 3:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5,
                        color="white" if v > vmax * 0.62 else "black")
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("% cells", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("A  Human liver (PSC/PBC): ligand expression", fontsize=7.5,
                 loc="left")
    return m


def panel_b(ax):
    d = pd.read_csv(Path(r"C:\dockwork\drg") / "drg_cluster_expression.csv")
    rec = ["PTGER3", "PTGER1", "IL6ST", "IL6R", "CALCRL", "TACR1", "TRPV1", "NTRK1"]
    m = (d[d.gene.isin(rec)].pivot_table(index="cluster", columns="gene",
                                         values="pct_expressing").reindex(columns=rec))
    clusters = sorted(m.index, key=lambda c: int(str(c).lstrip("H") or 0))
    m = m.reindex(index=clusters)
    vmax = float(np.nanmax(m.values))
    im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(rec)))
    ax.set_xticklabels(rec, rotation=45, ha="right", fontsize=6.5, style="italic")
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=6)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if np.isfinite(v) and v >= 5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5,
                        color="white" if v > vmax * 0.62 else "black")
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("% nuclei", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("B  Human DRG: receptor expression", fontsize=7.5, loc="left")
    return m


def panel_c(ax):
    d = pd.read_csv(ROOT / "human_sc" / "人源单细胞_细胞类型定位.csv")
    drg = pd.read_csv(Path(r"C:\dockwork\drg") / "drg_cluster_expression.csv")

    def top_of(gene):
        s = d[d.gene == gene].sort_values("pct_expressing", ascending=False).iloc[0]
        return shorten(s.cell_type), float(s.pct_expressing)

    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.4)
    ax.add_patch(plt.Rectangle((0.05, 0.45), 3.75, 3.45, fc="#FDF1F0",
                               ec="#B4544C", lw=0.8))
    ax.add_patch(plt.Rectangle((6.2, 0.45), 3.75, 3.45, fc="#EFF3FA",
                               ec="#3F5F93", lw=0.8))
    ax.text(0.18, 3.62, "Injured liver", fontsize=7.5, weight="bold", color="#8E2E28")
    ax.text(6.33, 3.62, "Human DRG nociceptor", fontsize=7.5, weight="bold",
            color="#2C4470")

    ligands = ["IL1B", "TNF", "PTGS2", "IL6"]
    ys = [2.92, 2.22, 1.52, 0.82]
    for g, y in zip(ligands, ys):
        ct, pct = top_of(g)
        ax.text(0.18, y + 0.09, f"{g}  ({ct})", fontsize=6.4, color="#5A1E1A")
        ax.text(0.18, y - 0.17, f"{pct:.1f}% of cells", fontsize=5.5,
                color="#7A6A68")

    ax.annotate("", xy=(6.1, 2.1), xytext=(3.9, 2.1),
                arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#B4544C"))
    ax.text(5.0, 2.32, "mediators", fontsize=6.2, ha="center", color="#B4544C")
    ax.text(5.0, 1.86, "blood / perineural", fontsize=5.4, ha="center",
            color="#8A6A66")

    def drg_mean(gene):
        return float(drg[drg.gene == gene].pct_expressing.mean())

    def drg_max(gene):
        s = drg[drg.gene == gene].sort_values("pct_expressing", ascending=False).iloc[0]
        return str(s.cluster), float(s.pct_expressing)

    r1c, r1v = drg_max("PTGER3")
    recs = [("PTGER3 (EP3)", f"{drg_mean('PTGER3'):.0f}% mean, {r1v:.0f}% in {r1c}"),
            ("IL6ST", f"{drg_mean('IL6ST'):.0f}% mean, IL6R only {drg_mean('IL6R'):.0f}%"),
            ("CALCRL", f"{drg_mean('CALCRL'):.0f}% mean"),
            ("TACR1 / TRPV1", f"{drg_mean('TACR1'):.0f}% / {drg_mean('TRPV1'):.0f}% mean")]
    for (name, txt), y in zip(recs, ys):
        ax.text(6.33, y + 0.09, name, fontsize=6.4, color="#1E2E4F")
        ax.text(6.33, y - 0.17, txt, fontsize=5.5, color="#5C6478")

    ax.text(0.05, 0.12, "% of cells expressing the gene; localisation, not flux",
            fontsize=5.2, color="#555555")
    ax.set_title("C  Sending and receiving compartments", fontsize=7.5, loc="left")


def main():
    fig = plt.figure(figsize=(7.2, 6.5), dpi=300)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.8],
                  width_ratios=[1, 1], hspace=0.30, wspace=0.46)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, :]))
    fig.subplots_adjust(left=0.13, right=0.955, top=0.955, bottom=0.055)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig_send_receive_panels.{ext}", dpi=300,
                    bbox_inches="tight")
    print("written", OUT / "fig_send_receive_panels.png")


if __name__ == "__main__":
    main()
