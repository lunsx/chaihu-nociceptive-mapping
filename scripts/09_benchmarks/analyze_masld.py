# -*- coding: utf-8 -*-
"""Human MASLD single-cell analysis: QC, clustering, marker annotation, target localization."""
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

W = Path(__file__).resolve().parent
H5DIR = W / "h5"
H5DIR.mkdir(exist_ok=True)

# --- h5 files were extracted separately from the partial tar ---
h5s = sorted(H5DIR.glob("*.h5"))
print("h5 files:", len(h5s))

# --- load per sample ---
adatas = []
for f in h5s:
    a = sc.read_10x_h5(str(f))
    a.var_names_make_unique()
    parts = f.name.split("_")
    sample = next((p for p in parts if p.startswith("S") and p[1:].isdigit()), f.stem)
    a.obs["sample"] = sample
    a.obs["condition"] = "HC" if "_HC_" in f.name else "MASLD"
    adatas.append(a)
adata = sc.concat(adatas, join="outer", fill_value=0, label="sample2", keys=[a.obs["sample"].iloc[0] for a in adatas])
print("merged", adata.shape)

# --- QC ---
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20].copy()
print("after QC", adata.shape)

# --- normalize + reduce ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat")
adata.raw = adata
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2)
print("clusters:", adata.obs.leiden.nunique())

# --- marker-based annotation ---
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
present = {k: [g for g in v if g in adata.raw.var_names] for k, v in MARKERS.items()}
sc.tl.score_genes(adata, gene_list=[g for v in present.values() for g in v], score_name="tmp")
scoremat = {}
for ctype, genes in present.items():
    if genes:
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{ctype}")
        scoremat[ctype] = adata.obs[f"score_{ctype}"]
scores = pd.DataFrame(scoremat)
cluster_map = {}
for cl in adata.obs.leiden.cat.categories:
    m = adata.obs.leiden == cl
    best = scores.loc[m.index[m]].mean().idxmax()
    cluster_map[cl] = best
adata.obs["cell_type"] = pd.Categorical(adata.obs.leiden.map(cluster_map))
print(json.dumps(cluster_map, indent=1))

# --- target localization ---
TARGETS = ["IL1B", "TNF", "PTGS2", "IL6", "PTGES", "PTGER3", "CALCRL", "TACR1", "TRPV1", "IL1A"]
raw = adata.raw.to_adata()
rows = []
for ctype in adata.obs.cell_type.cat.categories:
    sub = raw[adata.obs.cell_type == ctype]
    n = sub.n_obs
    for g in TARGETS:
        if g in sub.var_names:
            x = sub[:, g].X
            x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
            rows.append({"cell_type": ctype, "n_cells": n, "gene": g,
                         "pct_expressing": round(100 * (x > 0).mean(), 2),
                         "mean_log_expr": round(float(x.mean()), 3)})
out = pd.DataFrame(rows)
out.to_csv(W / "masld_target_localization.csv", index=False, encoding="utf-8-sig")
ct_counts = adata.obs.cell_type.value_counts()
ct_counts.to_csv(W / "masld_celltype_counts.csv", encoding="utf-8-sig")
summary = {"total_cells": int(adata.n_obs), "total_genes": int(adata.n_vars),
           "n_samples": int(adata.obs["sample"].nunique()),
           "conditions": adata.obs.condition.value_counts().to_dict(),
           "n_cell_types": int(adata.obs.cell_type.nunique()),
           "cluster_map": cluster_map}
json.dump(summary, open(W / "masld_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
drop_cols = [c for c in adata.obs.columns if c.startswith("score_") or c == "tmp"]
adata_write = adata.copy()
adata_write.obs.drop(columns=drop_cols, inplace=True)
adata_write.write_h5ad(W / "masld_processed.h5ad")
print("DONE")
