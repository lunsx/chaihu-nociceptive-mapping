# -*- coding: utf-8 -*-
"""GSE344087 MASLD full series: Harmony batch integration and marker heatmap.

The primary analysis (09_analyze_masld_full.py) concatenates 20 samples and
clusters without explicit batch correction. Reviewers can ask whether the
embedding is driven by sample/donor rather than biology. This script:

  1. reads the processed object in backed mode and pulls only the PCA embedding
     and the needed obs columns, so the 12.7 GB file is never materialized;
  2. runs Harmony (pure-Python harmonypy) on the 30 PCs with sample as the
     batch variable, fixed random state;
  3. rebuilds neighbors and UMAP on the corrected embedding;
  4. quantifies batch mixing before and after as the fraction of each cell's
     nearest neighbors that share its sample (lower = better mixing);
  5. computes per-cell-type marker means from the raw lognorm layer and renders
     a marker heatmap and a donor-composition stacked bar.

The single healthy-control donor still rules out inferential disease tests;
this is a technical integration check and descriptive annotation evidence,
not a new statistical comparison.

Path configuration: set CHAIHU_MASLD_DIR to the directory holding
masld_full_processed.h5ad; the local analysis path is the default.
"""
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from harmonypy import run_harmony

W = Path(os.environ.get("CHAIHU_MASLD_DIR", r"C:\dockwork\masld"))
SRC = W / "masld_full_processed.h5ad"
K = 15
SEED = 0

MARKERS = {
    "Hepatocyte": ["ALB", "TTR", "APOA1", "APOB"],
    "Endothelial": ["PECAM1", "KDR", "CLEC4G", "CLEC14A"],
    "Cholangiocyte": ["KRT19", "EPCAM", "KRT7"],
    "Macrophage/Kupffer": ["CD68", "LYZ", "C1QA", "MARCO", "C1QB"],
    "T/NK": ["NKG7", "CD3D", "CD3E", "GNLY", "CD8A"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "Plasma": ["MZB1", "JCHAIN", "IGHG1"],
    "Hepatic stellate": ["COL1A1", "ACTA2", "RGS5", "COL3A1"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
}


def same_sample_fraction(embedding, sample_labels, k=K):
    """Mean fraction of each cell's k neighbors sharing its sample."""
    n = embedding.shape[0]
    sub = sc.AnnData(X=np.zeros((n, 1)))
    sub.obsm["X_emb"] = np.asarray(embedding)
    sc.pp.neighbors(sub, n_neighbors=k + 1, use_rep="X_emb", random_state=SEED)
    conn = sub.obsp["connectivities"].tocsr()
    out = np.empty(n)
    for i in range(n):
        nb = conn[i].indices
        nb = nb[nb != i][:k]
        if len(nb) == 0:
            out[i] = np.nan
            continue
        out[i] = float(np.mean(sample_labels[nb] == sample_labels[i]))
    return out


def main():
    a0 = ad.read_h5ad(SRC, backed="r")
    pca = np.asarray(a0.obsm["X_pca"])
    obs = a0.obs[["sample", "donor", "condition", "cell_type", "leiden"]].copy()
    print("backed shape", a0.shape, "pca", pca.shape)

    # marker means from raw lognorm (only marker columns read into memory)
    raw = a0.raw
    raw_vars = raw.var_names.astype(str)
    gene_to_col = {g: j for j, g in enumerate(raw_vars)}
    marker_genes = []
    for genes in MARKERS.values():
        for g in genes:
            if g in gene_to_col and g not in marker_genes:
                marker_genes.append(g)
    cols = [gene_to_col[g] for g in marker_genes]
    raw_mat = raw.X[:, cols]
    raw_mat = raw_mat.toarray() if hasattr(raw_mat, "toarray") else np.asarray(raw_mat)
    marker_df = pd.DataFrame(raw_mat, columns=marker_genes)
    marker_df["cell_type"] = obs["cell_type"].astype(str).to_numpy()
    marker_means = marker_df.groupby("cell_type", observed=True).mean()
    marker_means.to_csv(W / "masld_marker_heatmap.csv", encoding="utf-8-sig")
    a0.file.close()

    sample = obs["sample"].astype(str).to_numpy()
    condition = obs["condition"].astype(str).to_numpy()

    # Harmony (pure-Python 0.0.6 API: matrix + metadata, returns Z_corr (p, n))
    meta_df = pd.DataFrame({"sample": sample})
    hp = run_harmony(pca, meta_df, ["sample"],
                     max_iter_harmony=20, random_state=SEED, verbose=False)
    z_corr = np.asarray(hp.Z_corr)
    harmony = z_corr.T if z_corr.shape[0] == pca.shape[1] else z_corr

    # embeddings and UMAP
    def build(emb, name):
        d = sc.AnnData(X=np.zeros((pca.shape[0], 1)))
        d.obsm["X_emb"] = np.asarray(emb)
        sc.pp.neighbors(d, n_neighbors=K, use_rep="X_emb", random_state=SEED)
        sc.tl.umap(d, random_state=SEED)
        return d, np.asarray(d.obsm["X_umap"])

    pre, pre_umap = build(pca, "pre")
    post, post_umap = build(harmony, "post")

    mix_pre = same_sample_fraction(pca, sample)
    mix_post = same_sample_fraction(harmony, sample)
    mix = pd.DataFrame({
        "condition": condition,
        "same_sample_fraction_pre": np.round(mix_pre, 4),
        "same_sample_fraction_post": np.round(mix_post, 4),
    })
    mix_summary = mix.groupby("condition", observed=True).mean()
    mix_summary.loc["ALL"] = mix[
        ["same_sample_fraction_pre", "same_sample_fraction_post"]].mean()
    mix_summary.to_csv(W / "masld_harmony_mixing.csv", encoding="utf-8-sig")
    print(mix_summary.round(4))

    np.savez_compressed(
        W / "masld_harmony_embeddings.npz",
        pca=pca, harmony=harmony, umap_pre=pre_umap, umap_post=post_umap)
    obs.to_csv(W / "masld_harmony_obs.csv", index=False, encoding="utf-8-sig")

    # UMAP panels colored by sample and condition, before vs after
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    panels = [
        (axes[0, 0], pre_umap, sample, "Before Harmony, colored by sample"),
        (axes[0, 1], post_umap, sample, "After Harmony, colored by sample"),
        (axes[1, 0], pre_umap, condition, "Before Harmony, colored by condition"),
        (axes[1, 1], post_umap, condition, "After Harmony, colored by condition"),
    ]
    cmap = plt.get_cmap("tab20")
    for ax, xy, color, title in panels:
        uniq = pd.unique(color)
        for ui, uval in enumerate(uniq):
            m = color == uval
            ax.scatter(xy[m, 0], xy[m, 1], s=2, color=cmap(ui % 20),
                       label=str(uval), alpha=0.6, linewidths=0)
        ax.set_title(title, fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(markerscale=6, fontsize=6, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(W / "masld_harmony_umap.png", dpi=200)
    plt.close(fig)

    # marker heatmap (z-score per gene across cell types)
    z = marker_means.sub(marker_means.mean(axis=0), axis=1).div(
        marker_means.std(axis=0).replace(0, np.nan), axis=1)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(z.shape[1])); ax.set_xticklabels(z.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(z.shape[0])); ax.set_yticklabels(z.index, fontsize=9)
    ax.set_title("Marker genes, z-scored mean log expression across cell types")
    fig.colorbar(im, ax=ax, fraction=0.025)
    fig.tight_layout()
    fig.savefig(W / "masld_marker_heatmap.png", dpi=200)
    plt.close(fig)

    # donor composition stacked bar
    comp = pd.read_csv(W / "masld_full_composition.csv")
    piv = comp.pivot_table(index="donor", columns="cell_type",
                           values="fraction", aggfunc="first", fill_value=0)
    cond_map = obs.drop_duplicates("donor").set_index("donor")["condition"].to_dict()
    piv = piv.loc[sorted(piv.index, key=lambda d: (cond_map.get(d, ""), d))]
    fig, ax = plt.subplots(figsize=(13, 6))
    bottom = np.zeros(len(piv))
    for ui, ct in enumerate(piv.columns):
        vals = piv[ct].to_numpy()
        ax.bar(range(len(piv)), vals, bottom=bottom, color=cmap(ui % 20), label=ct)
        bottom += vals
    ax.set_xticks(range(len(piv)))
    ax.set_xticklabels([f"{d}\n{cond_map.get(d,'')}" for d in piv.index],
                       rotation=90, fontsize=7)
    ax.set_ylabel("Cell-type fraction")
    ax.set_title("Donor-level cell-type composition (descriptive, no p-values)")
    ax.legend(fontsize=7, ncol=2, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(W / "masld_composition_stacked.png", dpi=200)
    plt.close(fig)

    summary = {
        "n_cells": int(pca.shape[0]),
        "n_samples": int(pd.unique(sample).shape[0]),
        "harmony_random_state": SEED,
        "same_sample_fraction_pre_all": float(np.nanmean(mix_pre)),
        "same_sample_fraction_post_all": float(np.nanmean(mix_post)),
        "mixing_by_condition": mix_summary.round(4).reset_index().to_dict("records"),
        "note": ("single healthy-control donor; descriptive integration check, "
                 "no inferential disease comparison"),
    }
    (W / "masld_harmony_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    print("DONE")


if __name__ == "__main__":
    main()
