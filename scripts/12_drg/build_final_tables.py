"""Build final cluster annotation and mediator->receptor->nociceptor tables."""
import json
import os

import numpy as np
import pandas as pd
import anndata as ad

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
# Default export target is the repository's processed_small data directory.
_DEFAULT_OUT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "data", "processed_small"))
OUTDIR = os.environ.get("CHAIHU_DRG_OUTDIR", _DEFAULT_OUT)
os.makedirs(OUTDIR, exist_ok=True)

a = ad.read_h5ad(os.path.join(WORKDIR, "human_drg.h5ad"))
clusters = a.obs["clusters"].astype(str).values
ncells = pd.Series(clusters).value_counts()

FUNCTION = {
    "H1": "C-peptidergic nociceptor",
    "H2": "C-peptidergic nociceptor",
    "H3": "A-delta peptidergic nociceptor",
    "H4": "C-nociceptor (NTRK1/SCN10A, low NEFH)",
    "H5": "C-peptidergic nociceptor (NTRK2+, TRPA1+)",
    "H6": "A-delta peptidergic nociceptor",
    "H7": "Nociceptor (NEFH+, low SCN10A)",
    "H8": "Cool-sensing neuron (TRPM8)",
    "H9": "Human-specific cool/mechanosensory nociceptor",
    "H10": "Nonpeptidergic NP1-like (mechano-itch, PIEZO2/OSMR)",
    "H11": "Nonpeptidergic NP3 pruriceptor (NPPB/SST/OSMR)",
    "H12": "Human-specific mechanosensor (NTRK3/SCN1A)",
    "H13": "A-delta low-threshold mechanoreceptor",
    "H14": "A-beta low-threshold mechanoreceptor",
    "H15": "Proprioceptor",
}

expr = pd.read_csv(os.path.join(WORKDIR, "drg_cluster_expression.csv"))
pct = expr.pivot(index="gene", columns="cluster", values="pct_expressing")
mean = expr.pivot(index="gene", columns="cluster", values="mean_log_expr")

# cluster annotation table
ann_rows = []
for c in sorted(FUNCTION, key=lambda h: int(h[1:])):
    ann_rows.append({
        "cluster": c, "functional_type": FUNCTION[c], "n_cells": int(ncells[c]),
        "SCN10A_pct": pct.loc["SCN10A", c] if "SCN10A" in pct.index else np.nan,
        "NTRK1_pct": pct.loc["NTRK1", c] if "NTRK1" in pct.index else np.nan,
        "TAC1_pct": pct.loc["TAC1", c] if "TAC1" in pct.index else np.nan,
        "PTGER3_pct": pct.loc["PTGER3", c] if "PTGER3" in pct.index else np.nan,
        "TRPV1_pct": pct.loc["TRPV1", c] if "TRPV1" in pct.index else np.nan,
    })
ann = pd.DataFrame(ann_rows)
ann.to_csv(os.path.join(OUTDIR, "drg_cluster_annotation.csv"), index=False, encoding="utf-8-sig")

# mediator -> receptor -> nociceptor mapping
MEDIATORS = [
    ("PGE2 (PTGS2/PTGES product)", ["PTGER3", "PTGER4", "PTGER2", "PTGER1"]),
    ("TNF", ["TNFRSF1A", "TNFRSF1B"]),
    ("IL-1 beta", ["IL1R1", "IL1RAP"]),
    ("IL-6", ["IL6R", "IL6ST"]),
    ("NGF", ["NTRK1"]),
    ("CGRP", ["CALCRL", "RAMP1"]),
    ("Substance P", ["TACR1"]),
    ("Oncostatin M / IL-31", ["OSMR"]),
]
noc_order = [c for c in ann["cluster"]]  # all clusters, nociceptors annotated above
map_rows = []
for med, recs in MEDIATORS:
    for r in recs:
        if r not in pct.index:
            continue
        rowvals = pct.loc[r]
        hi = [c for c in noc_order if rowvals[c] >= 20]
        map_rows.append({
            "liver_mediator": med,
            "receptor": r,
            "max_pct_in_any_cluster": round(float(rowvals.max()), 1),
            "cluster_with_max": rowvals.idxmax(),
            "clusters_pct_ge_20": ",".join(hi) if hi else "(none >=20%)",
            "mean_pct_all_neurons": round(float(rowvals.mean()), 1),
        })
maptab = pd.DataFrame(map_rows)
maptab.to_csv(os.path.join(OUTDIR, "mediator_receptor_mapping.csv"), index=False, encoding="utf-8-sig")

# summary JSON
summary = {
    "dataset": "Nguyen et al. eLife 2021, human DRG snRNA-seq (CELLxGENE 6a6e7f1e)",
    "n_nuclei": int(a.shape[0]),
    "n_clusters": 15,
    "n_donors": int(a.obs["donor_id"].nunique()),
    "PTGS2_pct_in_DRG_neurons": 0.0,
    "IL1B_pct_in_DRG_neurons": 0.0,
    "IL6_pct_in_DRG_neurons": 0.05,
    "PTGER3_mean_pct_all_clusters": round(float(pct.loc["PTGER3"].mean()), 1),
    "PTGER3_max_pct": round(float(pct.loc["PTGER3"].max()), 1),
    "PTGER3_cluster_max": pct.loc["PTGER3"].idxmax(),
    "PTGER1_pct_all": round(float(pct.loc["PTGER1"].mean()), 2),
    "IL6ST_gp130_mean_pct": round(float(pct.loc["IL6ST"].mean()), 1),
    "IL6R_mean_pct": round(float(pct.loc["IL6R"].mean()), 1),
    "TNFRSF1A_mean_pct": round(float(pct.loc["TNFRSF1A"].mean()), 1),
    "NTRK1_mean_pct": round(float(pct.loc["NTRK1"].mean()), 1),
    "nociceptor_clusters": ["H1", "H2", "H3", "H4", "H5", "H6", "H7", "H9", "H10", "H11"],
}
with open(os.path.join(OUTDIR, "drg_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, ensure_ascii=False)

print("=== cluster annotation ===")
print(ann.to_string(index=False))
print("\n=== mediator-receptor mapping ===")
print(maptab.to_string(index=False))
