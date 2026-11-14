# -*- coding: utf-8 -*-
"""All-section human Visium analysis: 4 PSC011 sections vs 4 C73 sections.

For each section we compute, on in-tissue spots with self-made lognorm:
  - positive-spot fraction and mean normalized expression for each target
  - Moran's I with permutation test (200 perms)
Then aggregate section-level values within donor group and report the
descriptive contrast. No inferential p-values across sections: each group has
one donor, so sections are technical replicates, not independent subjects.

Path configuration: set the environment variable CHAIHU_SPATIAL_DIR to the
directory holding the section h5ad files; the cloud analysis path is the
default. Output tables are written to the same directory.
"""
import gc
import os
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import squidpy as sq

BASE = Path(os.environ.get(
    "CHAIHU_SPATIAL_DIR", "/home/user/compute_work/spatial"))
SEC = BASE / "sections"
TARGETS = ["PTGS2", "IL6", "IL1B", "IL1A", "PTGES", "PTGER3", "TNF", "TRPV1",
           "TRPA1", "TACR1", "CALCRL", "NGF", "SCN9A", "CCL2", "CXCL10", "CXCL2",
           "FOS", "ICAM1", "IL10", "MMP9", "SERPINE1", "SPP1", "HMOX1",
           "ADRB2", "CREB1", "MYC", "PLAU", "F3", "GJA1"]


def load(fname):
    # backed read so the full raw matrix is never materialized
    a0 = ad.read_h5ad(fname, backed="r")
    if "in_tissue" in a0.obs.columns:
        spot_idx = np.flatnonzero(
            a0.obs["in_tissue"].astype(str).values == "1")
    else:
        spot_idx = np.arange(a0.n_obs)
    raw = a0.raw
    sym = (raw.var["feature_name"].astype(str).values
           if "feature_name" in raw.var.columns
           else raw.var_names.astype(str))
    # resolve duplicate gene symbols, keep first occurrence
    keep_gene = ~pd.Series(sym).duplicated().values
    sym_to_col = {sym[c]: c for c in range(len(sym)) if keep_gene[c]}
    target_cols = [sym_to_col[g] for g in TARGETS if g in sym_to_col]
    present = [sym[c] for c in target_cols]
    # library size accumulated in gene chunks to cap memory
    n = len(spot_idx)
    lib = np.zeros(n, dtype=np.float64)
    total_genes = raw.X.shape[1]
    step = 2000
    for c0 in range(0, total_genes, step):
        block = raw.X[spot_idx, c0:c0 + step]
        block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
        lib += block.sum(1)
    # materialize only the small target-gene submatrix
    counts = raw.X[spot_idx][:, target_cols]
    counts = counts.tocsr() if hasattr(counts, "tocsr") else sparse.csr_matrix(counts)
    lib[lib == 0] = 1.0
    inv = sparse.diags(1.0 / lib)
    norm = (inv @ counts) * 1e4
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(norm.toarray())))
    a.var_names = pd.Index(present)
    a.layers["counts"] = counts
    has_coords = "spatial" in a0.obsm
    if has_coords:
        a.obsm["spatial"] = np.asarray(a0.obsm["spatial"][spot_idx])
        sq.gr.spatial_neighbors(a, spatial_key="spatial", coord_type="generic",
                                delaunay=True, n_rings=1)
        sq.gr.spatial_autocorr(a, mode="moran", genes=present, n_perms=200, n_jobs=1)
        moran = a.uns["moranI"].copy()
        moran["gene"] = moran.index
    else:
        moran = pd.DataFrame({"gene": present, "I": np.nan,
                              "pval_sim": np.nan}).set_index("gene")
    a0.file.close()
    gc.collect()
    return a, present, moran


def main():
    files = sorted(SEC.glob("*.h5ad"))
    spot_rows, moran_rows = [], []
    for f in files:
        if f.name.startswith("NORMAL_pooled"):
            donor, tag = "NOR", "e84dbcac"
        else:
            donor, tag = f.name[:3], f.name.split("_")[1][:8]
        a, present, moran = load(f)
        moran["donor"], moran["section"] = donor, tag
        moran_rows.append(moran.reset_index(drop=True))
        for g in present:
            j = list(a.var_names).index(g)
            c = a.layers["counts"][:, j]
            c = c.toarray().ravel() if hasattr(c, "toarray") else np.asarray(c).ravel()
            x = a.X[:, j]
            x = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
            spot_rows.append({"donor": donor, "section": tag, "gene": g,
                              "n_spots": a.n_obs,
                              "pct_positive": round(float(np.mean(c > 0) * 100), 2),
                              "mean_norm": round(float(x.mean()), 4),
                              "moran_I": round(float(moran.loc[g, "I"]), 3),
                              "moran_p": round(float(moran.loc[g, "pval_sim"]), 4)})
        print(f"done {f.name[:45]} spots={a.n_obs}", flush=True)
        del a, moran
        gc.collect()

    spot_df = pd.DataFrame(spot_rows)
    moran_df = pd.concat(moran_rows, ignore_index=True)
    spot_df.to_csv(BASE / "空间_全切片_逐切片.csv", index=False, encoding="utf-8-sig")
    moran_df.to_csv(BASE / "空间_全切片_Moran.csv", index=False, encoding="utf-8-sig")

    # aggregate within donor
    agg = (spot_df.groupby(["donor", "gene"])
           .agg(sections=("section", "nunique"),
                mean_pct_positive=("pct_positive", "mean"),
                mean_norm=("mean_norm", "mean"),
                mean_moran_I=("moran_I", "mean"))
           .reset_index())
    piv = agg.pivot(index="gene", columns="donor",
                    values=["mean_pct_positive", "mean_norm", "mean_moran_I"])
    piv.columns = [f"{m}_{d}" for m, d in piv.columns]
    piv = piv.reset_index()
    if "mean_norm_C73" in piv and "mean_norm_PSC" in piv:
        piv["delta_norm_PSC_minus_C73"] = (piv["mean_norm_PSC"] - piv["mean_norm_C73"]).round(4)
    piv = piv.round(3).sort_values("delta_norm_PSC_minus_C73", ascending=False, na_position="last")
    piv.to_csv(BASE / "空间_全切片_组间汇总.csv", index=False, encoding="utf-8-sig")
    print("\n=== group summary (top by delta) ===")
    print(piv.head(12).to_string(index=False))
    print("\nNote: one donor per group; sections are technical replicates, descriptive only.")


if __name__ == "__main__":
    main()
