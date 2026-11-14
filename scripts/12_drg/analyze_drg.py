"""Per-cluster expression of inflammatory mediator receptors and nociceptor markers
in human DRG (Nguyen eLife 2021; 1837 nuclei, 15 clusters, log-normalized X)."""
import os

import numpy as np
import pandas as pd
import anndata as ad

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
a = ad.read_h5ad(os.path.join(WORKDIR, "human_drg.h5ad"))
X = a.X  # sparse CSR, log-normalized
clusters = a.obs["clusters"].astype(str).values
uniq = sorted(set(clusters), key=lambda h: int(h[1:]))

# symbol -> row index (first occurrence)
sym_to_idx = {}
for i, sym in enumerate(a.var["hgnc_symbol"].astype(str)):
    sym_to_idx.setdefault(sym, i)

# gene panels
MARKERS = {
    "lineage": ["SCN10A", "SCN11A", "NTRK1", "TAC1", "CALCA", "CALCB", "NPPB",
                "P2RX3", "GFRA1", "GFRA2", "PIEZO2", "NEFH", "PVALB", "TRPM8",
                "NTRK2", "NTRK3", "PENK", "TAC1"],
    "inflammatory_receptors": ["IL1R1", "IL1RAP", "TNFRSF1A", "TNFRSF1B",
                                "IL6R", "IL6ST", "PTGER1", "PTGER2", "PTGER3",
                                "PTGER4", "CALCRL", "RAMP1", "TACR1", "OSMR",
                                "IL13RA1", "JAK1"],
    "channels": ["TRPV1", "TRPA1", "SCN9A", "SCN10A", "SCN11A", "KCNQ2", "KCNQ3"],
}
all_genes = sorted(set(g for v in MARKERS.values() for g in v))
rows = []
for g in all_genes:
    if g not in sym_to_idx:
        continue
    j = sym_to_idx[g]
    col = np.asarray(X[:, j].todense()).ravel()
    for c in uniq:
        m = clusters == c
        vals = col[m]
        rows.append({
            "gene": g, "cluster": c, "n_cells": int(m.sum()),
            "mean_log_expr": round(float(vals.mean()), 4),
            "pct_expressing": round(float(np.mean(vals > 0) * 100), 2),
            "mean_expr_when_on": round(float(vals[vals > 0].mean()) if np.any(vals > 0) else 0, 4),
        })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(WORKDIR, "drg_cluster_expression.csv"), index=False, encoding="utf-8-sig")
print("saved drg_cluster_expression.csv, rows:", len(df))

# compact pivot of pct expressing for the inflammatory receptors
pct = df.pivot(index="gene", columns="cluster", values="pct_expressing")
order = [c for c in uniq if c in pct.columns]
print("\n=== % expressing: inflammatory mediator receptors across clusters ===")
with pd.option_context("display.width", 250, "display.max_columns", 20):
    print(pct.reindex(columns=order).loc[MARKERS["inflammatory_receptors"]].round(1).to_string())
