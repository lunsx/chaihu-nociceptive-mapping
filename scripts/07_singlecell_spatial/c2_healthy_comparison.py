# -*- coding: utf-8 -*-
"""C2: external healthy human liver reference comparison.

The Aizarani healthy liver atlas (CELLxGENE 135bdd30, 73,295 cells, four
donors) is used as an independent external normal reference. Per-cell-type
percentages of expressing cells (raw count > 0) are computed for the target
genes and mapped to the comparable MASLD labels. The two collections differ in
donors, platform, processing, and annotation, so this is a descriptive
external-reference comparison, not a matched statistical test, and no
cross-dataset p-values are computed.

Inputs are resolved from config / environment (no hard-coded machine paths).
Writes ``data/processed_small/external_healthy_reference_comparison.csv``.

Run:
  python scripts/07_singlecell_spatial/c2_healthy_comparison.py
"""
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

HEALTHY = Path(os.environ.get(
    "CHAIHU_HEALTHY_H5AD",
    PATHS["healthy_liver"] / "healthy_liver_integrated.h5ad",
))
MASLD = PATHS["masld"] / "masld_full_localization.csv"
OUT = PATHS["processed_small"] / "external_healthy_reference_comparison.csv"

# Canonical cell type -> comparable MASLD label.
TYPE_MAP = {
    "Kupffer cell": "Macrophage/Kupffer",
    "M1 macrophage": "Macrophage/Kupffer",
    "hepatic stellate cell": "Hepatic stellate",
    "endothelial cell of pericentral hepatic sinusoid": "Endothelial",
    "blood vessel endothelial cell": "Endothelial",
    "B cell": "B",
    "natural killer cell": "T/NK",
    "cholangiocyte": "Cholangiocyte",
    "periportal region hepatocyte": "Hepatocyte",
    "midzonal region hepatocyte": "Hepatocyte",
    "centrilobular region hepatocyte": "Hepatocyte",
    "erythroid lineage cell": "Erythroid",
}

TARGETS = [
    "IL1B", "TNF", "PTGS2", "IL6", "IL1A", "PTGES", "PTGER1", "PTGER3",
    "CALCRL", "TACR1", "TRPV1", "TNFRSF1A", "IL6R", "IL6ST", "ICAM1",
    "SPP1", "FOS",
]


def read_str(node):
    if isinstance(node, h5py.Dataset):
        return node.asstr()[:]
    if "categories" in node:
        cats = node["categories"].asstr()[:]
        codes = node["codes"][:]
        return np.array([cats[c] if c >= 0 else "" for c in codes])
    raise TypeError(node.name)


def main():
    with h5py.File(HEALTHY, "r") as h:
        feat = read_str(h["var"]["feature_name"])
        sym2cols = {}
        for g in TARGETS:
            cols = np.where(feat == g)[0]
            if len(cols):
                sym2cols[g] = cols
        print("targets present:", len(sym2cols), "missing:",
              [g for g in TARGETS if g not in sym2cols])

        wanted = np.unique(np.concatenate(list(sym2cols.values())))
        col2local = {c: i for i, c in enumerate(wanted)}

        rX = h["raw"]["X"]
        data, indices, indptr = rX["data"], rX["indices"], rX["indptr"]
        n = indptr.shape[0] - 1
        sub = np.zeros((n, len(wanted)))
        chunk = 4000
        for s in range(0, n, chunk):
            e = min(s + chunk, n)
            a, b = indptr[s], indptr[e]
            d = data[a:b].astype(np.float64)
            ix = indices[a:b]
            rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
            mask = np.isin(ix, wanted)
            if mask.any():
                loc = np.array([col2local[c] for c in ix[mask]])
                np.add.at(sub, (rows[mask], loc), d[mask])
        ctype = read_str(h["obs"]["cell_type"])
        mapped = np.array([TYPE_MAP.get(t, "OTHER") for t in ctype])

    gene_mat = {
        g: sub[:, [col2local[c] for c in cols]].sum(axis=1)
        for g, cols in sym2cols.items()
    }

    healthy_rows = []
    for g, vals in gene_mat.items():
        for t in sorted(set(TYPE_MAP.values())):
            idx = np.where(mapped == t)[0]
            if len(idx) == 0:
                continue
            pct = float((vals[idx] > 0).mean() * 100)
            healthy_rows.append({
                "gene": g, "cell_type": t, "n_healthy_cells": len(idx),
                "healthy_pct": round(pct, 2),
            })
    hdf = pd.DataFrame(healthy_rows)

    m = pd.read_csv(MASLD)
    md = m[m["condition"] == "MASLD"][[
        "gene", "cell_type", "n_donors", "pct_donor_median",
        "pct_donor_min", "pct_donor_max",
    ]]
    hc = m[m["condition"] == "HC"][["gene", "cell_type", "pct_donor_median"]]
    hc = hc.rename(columns={"pct_donor_median": "single_GSE344087_HC_pct"})

    out = hdf.merge(md, on=["gene", "cell_type"], how="left")
    out = out.merge(hc, on=["gene", "cell_type"], how="left")
    out = out.rename(columns={
        "pct_donor_median": "MASLD_donor_median_pct",
        "pct_donor_min": "MASLD_donor_min_pct",
        "pct_donor_max": "MASLD_donor_max_pct",
    })
    out.to_csv(OUT, index=False)
    print("saved", OUT, "rows:", len(out))


if __name__ == "__main__":
    main()
