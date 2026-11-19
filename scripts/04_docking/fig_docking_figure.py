# -*- coding: utf-8 -*-
"""Build the two-panel molecular docking figure (main Figure 8).

Reads the version-controlled docking table and renders:
  A. grouped binding-energy bar chart
  B. multi-receptor heatmap

The heatmap previously used RdYlGn_r, which is not distinguishable for
red-green color vision deficiency and loses ordering in grayscale. It is
replaced by viridis, which is colorblind-safe and luminance-ordered. Numeric
labels keep every value readable directly, so the message never relies on
color alone.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Times New Roman"]

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(REPO, "data", "processed_small", "Table4_对接结合能.csv")
OUTDIR = os.path.abspath(os.path.join(REPO, "..", "投稿图集"))
OUT = os.path.join(OUTDIR, "Fig7_分子对接.png")

rows = list(csv.DictReader(open(DATA, encoding="utf-8-sig")))
recs = sorted(set(r["receptor"] for r in rows))
lig_order = ["cyperotundone", "kaempferol", "luteolin", "paeoniflorin", "quercetin"]
ligs = [l for l in lig_order if l in set(r["ligand"] for r in rows)]
M = np.full((len(recs), len(ligs)), np.nan)
for r in rows:
    i = recs.index(r["receptor"]); j = ligs.index(r["ligand"])
    M[i, j] = float(r["best_affinity_kcal_mol"])

receptor_sub = {"2AZ5": "TNF", "5KIR": "COX-2"}

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.05, 1]})
fig.suptitle("Molecular docking binding energies", fontsize=16)

# Panel A: grouped bars.
bar_colors = ["#1f77b4", "#ff7f0e"]
x = np.arange(len(ligs)); w = 0.38
for k, rc in enumerate(recs):
    ax1.bar(x + (k - 0.5) * w, M[k], w, label=rc, color=bar_colors[k])
    for xi, v in zip(x + (k - 0.5) * w, M[k]):
        ax1.text(xi, v - 0.12, f"{v:.2f}", ha="center", va="top", fontsize=9)
ax1.axhline(-7, ls="--", c="gray", lw=1)
ax1.text(0.15, -7.05, "-7.0 reference", ha="left", va="bottom",
         fontsize=9, color="gray")
ax1.set_xticks(x); ax1.set_xticklabels(ligs, rotation=20, ha="right")
ax1.set_ylabel("Best affinity (kcal/mol)")
ax1.set_ylim(-8.0, -4.5)
ax1.set_title("A. Binding-energy bar chart", loc="left", fontweight="bold", fontsize=13)
ax1.legend()

# Panel B: colorblind-safe heatmap.
vmin, vmax = -8.0, -5.0
im = ax2.imshow(M, cmap="viridis", aspect="auto", vmin=vmin, vmax=vmax)
ax2.set_xticks(range(len(ligs))); ax2.set_xticklabels(ligs, rotation=30, ha="right")
ax2.set_yticks(range(len(recs)))
ax2.set_yticklabels([f"{r}\n({receptor_sub.get(r, '')})" for r in recs])
for i in range(len(recs)):
    for j in range(len(ligs)):
        norm = (M[i, j] - vmin) / (vmax - vmin)
        r0, g0, b0, _ = im.cmap(norm)
        lum = 0.299 * r0 + 0.587 * g0 + 0.114 * b0
        txt_color = "black" if lum > 0.5 else "white"
        ax2.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                 fontsize=10, fontweight="bold", color=txt_color)
cb = fig.colorbar(im, ax=ax2, shrink=0.85)
cb.set_label("Best binding affinity (kcal/mol)")
ax2.set_title("B. Multi-receptor heatmap", loc="left", fontweight="bold", fontsize=13)

plt.tight_layout(rect=[0, 0, 1, 0.96])
os.makedirs(OUTDIR, exist_ok=True)
plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("saved", OUT)
