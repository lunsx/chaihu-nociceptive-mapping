"""Explore the human DRG h5ad structure."""
import os

import anndata as ad

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
a = ad.read_h5ad(os.path.join(WORKDIR, "human_drg.h5ad"))
print("shape:", a.shape)
print("\nobs columns:", list(a.obs.columns))
for c in a.obs.columns:
    vals = a.obs[c].astype(str)
    if vals.nunique() <= 25:
        print(f"\n{c} ({vals.nunique()} levels):")
        print(vals.value_counts().to_string())
print("\nvar columns:", list(a.var.columns))
print("var index sample:", list(a.var_names[:8]))
print("\nX type:", type(a.X), "dtype:", a.X.dtype)
print("layers:", list(a.layers.keys()))
