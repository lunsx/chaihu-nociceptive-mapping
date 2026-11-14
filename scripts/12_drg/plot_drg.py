"""Publication-style dot plot: mediator receptors across human DRG clusters."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
_DEFAULT_FIG = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "figures"))
FIGDIR = os.environ.get("CHAIHU_FIGDIR", _DEFAULT_FIG)
os.makedirs(FIGDIR, exist_ok=True)

expr = pd.read_csv(os.path.join(WORKDIR, "drg_cluster_expression.csv"))
pct = expr.pivot(index="gene", columns="cluster", values="pct_expressing")
mean = expr.pivot(index="gene", columns="cluster", values="mean_log_expr")
order = sorted(pct.columns, key=lambda h: int(h[1:]))

genes = ["PTGER3", "PTGER4", "PTGER2", "PTGER1", "NTRK1", "TNFRSF1A",
         "IL6ST", "IL6R", "IL1R1", "IL1RAP", "OSMR", "TACR1", "CALCRL",
         "TRPV1", "TRPA1", "SCN9A", "SCN10A"]
genes = [g for g in genes if g in pct.index]

fig, ax = plt.subplots(figsize=(11, 6.2))
P = pct.reindex(index=genes, columns=order)
M = mean.reindex(index=genes, columns=order)
x = np.arange(len(order))
y = np.arange(len(genes))
for i, g in enumerate(genes):
    for j, c in enumerate(order):
        p = P.iloc[i, j]
        m = M.iloc[i, j]
        size = 18 + (p / 100.0) * 260
        color = plt.cm.viridis(np.clip(m / 2.2, 0, 1))
        ax.scatter(j, i, s=size, color=color, edgecolor="0.3", linewidth=0.5, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(order)
ax.set_yticks(y)
ax.set_yticklabels(genes)
ax.set_xlabel("Human DRG cluster (Nguyen et al. 2021)")
ax.set_ylabel("Receptor / channel gene")
ax.set_title("Inflammatory mediator receptors and excitatory channels in human DRG neurons\n"
             "dot size = % expressing; color = mean log expression", fontsize=11)
ax.grid(axis="x", color="0.9", lw=0.7, zorder=0)
ax.set_axisbelow(True)
for spine in ax.spines.values():
    spine.set_color("0.7")

# size legend (outside axes, marker area proportional to pct)
from matplotlib.lines import Line2D
handles = []
for pval in [25, 50, 75, 100]:
    handles.append(Line2D([0], [0], marker="o", linestyle="none",
                          markersize=np.sqrt(18 + pval / 100 * 260) / 2.1,
                          markerfacecolor="white", markeredgecolor="0.4",
                          label=f"{pval}%"))
ax.legend(handles=handles, title="% expressing", loc="upper left",
          bbox_to_anchor=(1.01, 1.0), frameon=False, handletextpad=1.2)
fig.tight_layout()
out = os.path.join(FIGDIR, "drg_receptor_dotplot.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print("saved drg_receptor_dotplot.png")
