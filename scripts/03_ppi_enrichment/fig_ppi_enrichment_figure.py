# -*- coding: utf-8 -*-
"""Build the PPI network and functional enrichment figure (main Figure 2).

Reads the committed topology and enrichment tables and renders three panels:
  A. top-60 STRING PPI network;
  B. formula-wide enrichment (GO BP/MF, KEGG, Reactome);
  C. GO cellular component.

Pain/inflammatory nodes were previously distinguished only by red fill. They
now additionally carry a black outline and a larger size, so the highlight
survives red-green color vision deficiency and grayscale printing. Only
numpy, pandas, networkx and matplotlib are required.
"""
import os
import numpy as np, pandas as pd
import networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Times New Roman"]

ROOT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
PP = os.path.join(ROOT, "ppi")
EN = os.path.join(ROOT, "enrichment", "formula237")
OUT = os.path.join(ROOT, "投稿图集", "Fig4_PPI与富集.png")

PAIN = {"TRPV1", "PTGS2", "TNF", "IL1B", "IL6", "IL1A", "PTGES", "PTGER3",
        "PTGS1", "TRPA1", "TACR1", "CALCRL", "NGF", "SCN9A", "SCN5A",
        "SCN10A", "TRPV3", "TRPC3", "IL1R1"}
BLUE = "#4A90D9"; RED = "#C0392B"

nodes = pd.read_csv(os.path.join(PP, "PPI_节点拓扑.csv"))
edges = pd.read_csv(os.path.join(PP, "PPI_边表.csv"))
top = nodes.head(60)
top_genes = list(top.gene); top_set = set(top_genes)
deg = dict(zip(nodes.gene, nodes.degree))
sub = nx.Graph(); sub.add_nodes_from(top_genes)
for a, b in zip(edges.A, edges.B):
    if a in top_set and b in top_set:
        sub.add_edge(a, b)
pos = nx.spring_layout(sub, seed=7, k=0.55)

fig = plt.figure(figsize=(22, 9))
outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.62], wspace=0.16)

# ---------- Panel A: PPI network ----------
axA = fig.add_subplot(outer[0]); axA.axis("off")
dmax = max(deg[g] for g in top_genes)
base_sizes = {g: 300 + 1700 * deg[g] / dmax for g in top_genes}
pain_list = [g for g in top_genes if g in PAIN]
other_list = [g for g in top_genes if g not in PAIN]
nx.draw_networkx_edges(sub, pos, ax=axA, alpha=0.25, edge_color="#7F8C8D")
nx.draw_networkx_nodes(sub, pos, ax=axA, nodelist=other_list,
                       node_size=[base_sizes[g] for g in other_list],
                       node_color=BLUE, alpha=0.92)
nx.draw_networkx_nodes(sub, pos, ax=axA, nodelist=pain_list,
                       node_size=[base_sizes[g] * 1.28 for g in pain_list],
                       node_color=RED, alpha=0.95,
                       edgecolors="black", linewidths=1.8)
nx.draw_networkx_labels(sub, pos, labels={g: g for g in pain_list}, ax=axA,
                        font_size=10, font_weight="bold", font_color="white")
nx.draw_networkx_labels(sub, pos, labels={g: g for g in other_list}, ax=axA,
                        font_size=9)
axA.legend(handles=[
    Line2D([0], [0], marker="o", color="w", markerfacecolor=RED,
           markeredgecolor="black", markersize=12, label="Pain/inflammatory gene"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE,
           markersize=12, label="Other compound target")],
    loc="lower left", fontsize=9, frameon=True)
axA.text(0.0, 1.06, "A. PPI network", transform=axA.transAxes,
         fontweight="bold", fontsize=15, va="bottom")
axA.text(0.0, 1.015,
         "Compound-target PPI network (top 60 core nodes; red = pain/inflammatory genes)",
         transform=axA.transAxes, fontsize=9.5, va="bottom")

# ---------- Right column: B and C ----------
right = outer[1].subgridspec(2, 1, height_ratios=[1.3, 1.0], hspace=0.42)
gridB = right[0].subgridspec(2, 2, wspace=0.34, hspace=0.62)
gridC = right[1].subgridspec(1, 1)


def plot_terms(ax, filename, short, top_n):
    r = pd.read_csv(os.path.join(EN, filename))
    r = r[r["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").head(top_n)
    if not len(r):
        ax.set_visible(False); return
    y = np.arange(len(r))[::-1]
    ax.barh(y, -np.log10(r["Adjusted P-value"]), color=BLUE)
    ax.set_yticks(y)
    ax.set_yticklabels([t if len(t) <= 58 else t[:56] + "..." for t in r.Term],
                       fontsize=7)
    ax.set_xlabel("-log10(adj P)", fontsize=8)
    ax.set_title(f"{short} Top terms", fontsize=9.5)
    ax.tick_params(axis="x", labelsize=7)


b_files = [
    ("GO_BP", "GO_BP_GO_Biological_Process_2023.csv"),
    ("GO_MF", "GO_MF_GO_Molecular_Function_2023.csv"),
    ("KEGG", "KEGG_KEGG_2021_Human.csv"),
    ("Reactome", "Reactome_Reactome_2022.csv"),
]
for i, (short, fname) in enumerate(b_files):
    ax = fig.add_subplot(gridB[i // 2, i % 2])
    plot_terms(ax, fname, short, 10)
    if i == 0:
        ax.text(0.0, 1.42, "B. Formula-wide enrichment", transform=ax.transAxes,
                fontweight="bold", fontsize=15, va="bottom")
        ax.text(0.0, 1.18,
                "Functional enrichment of compound targets (237 human genes from pass-filter compounds)",
                transform=ax.transAxes, fontsize=9.5, va="bottom")

axC = fig.add_subplot(gridC[0])
plot_terms(axC, "GO_CC_GO_Cellular_Component_2023.csv", "GO_CC", 12)
axC.text(0.0, 1.28, "C. GO cellular component", transform=axC.transAxes,
         fontweight="bold", fontsize=15, va="bottom")

fig.suptitle("PPI network and functional enrichment", fontsize=18, y=1.0)
plt.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print("saved", OUT)
