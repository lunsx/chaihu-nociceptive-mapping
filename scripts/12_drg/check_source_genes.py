"""Confirm the inflammatory ligands/enzymes are mostly NOT produced by DRG neurons,
while their receptors are present (peripheral-source argument)."""
import os

import numpy as np
import anndata as ad

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
a = ad.read_h5ad(os.path.join(WORKDIR, "human_drg.h5ad"))
X = a.X
sym_to_idx = {}
for i, sym in enumerate(a.var["hgnc_symbol"].astype(str)):
    sym_to_idx.setdefault(sym, i)

genes = ["PTGS2", "PTGES", "PTGES2", "PTGES3", "IL1B", "IL1A", "TNF", "IL6",
         "PTGER3", "TNFRSF1A", "IL6ST", "IL6R", "IL1R1", "NGF", "NTRK1",
         "CALCA", "CALCB", "TAC1", "PENK"]
print(f"{'gene':9s} {'pct_all_%':>9} {'mean_log':>9}")
for g in genes:
    if g not in sym_to_idx:
        print(f"{g:9s} absent")
        continue
    j = sym_to_idx[g]
    col = np.asarray(X[:, j].todense()).ravel()
    print(f"{g:9s} {np.mean(col>0)*100:9.2f} {col.mean():9.4f}")
