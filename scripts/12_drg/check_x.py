"""Check X normalization and locate key genes, plus marker-based cluster hints."""
import os

import numpy as np
import anndata as ad

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
a = ad.read_h5ad(os.path.join(WORKDIR, "human_drg.h5ad"))
X = a.X
print("X max:", X.max(), "min:", X.min())
print("nCount_RNA vs X row sum (first 3):")
print("nCount:", a.obs["nCount_RNA"].values[:3])
print("X sums:", np.asarray(X[:3].sum(1)).ravel())
# Is it integer-ish?
sample = X.data[:1000]
print("integer fraction:", np.mean(np.isclose(sample, np.round(sample))))

# map hgnc symbol -> index
sym_to_idx = {}
for i, sym in enumerate(a.var["hgnc_symbol"].astype(str)):
    sym_to_idx.setdefault(sym, i)

targets = ["IL1R1", "IL1RAP", "TNFRSF1A", "TNFRSF1B", "IL6R", "IL6ST",
           "PTGER1", "PTGER2", "PTGER3", "PTGER4", "CALCRL", "RAMP1",
           "TACR1", "TRPV1", "TRPA1", "NTRK1", "SCN9A", "SCN10A",
           "PIEZO2", "PENK", "TAC1", "NPPB", "P2RX3", "GFRA1", "GFRA2"]
print("\nkey gene presence:")
for g in targets:
    print(g, "yes" if g in sym_to_idx else "NO", end=" | ")
print()
