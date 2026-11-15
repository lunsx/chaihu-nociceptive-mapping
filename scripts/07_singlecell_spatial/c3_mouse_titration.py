# -*- coding: utf-8 -*-
"""C3: third independent pseudobulk titration on the mouse MASH collection.

The mouse MASH LoxP/SIRT2 collection (singlecell/MASH_liver.h5ad, 21,715 cells,
14 cell types) is run through the same controlled titration as the human MASLD
and PSC/PBC collections, giving a cross-species replication of the bulk-averaging
dilution mechanism. Raw counts are in raw.X (CSR). The 60 human target symbols
are mapped to mouse symbols by ortholog capitalization; CXCL11 has no mouse
ortholog and is excluded.

Inputs are resolved from config / environment. Writes
``mouse_pseudobulk_titration.csv`` and ``mouse_pseudobulk_titration_summary.json``
to ``data/processed_small/``.

Run:
  python scripts/07_singlecell_spatial/c3_mouse_titration.py
"""
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

H5AD = Path(os.environ.get(
    "CHAIHU_MOUSE_H5AD",
    PATHS["singlecell"] / "MASH_liver.h5ad",
))
GENE_LIST = PATHS["processed_small"] / "target_gene_list_60.csv"
OUT = PATHS["processed_small"] / "mouse_pseudobulk_titration.csv"
OUT_JSON = PATHS["processed_small"] / "mouse_pseudobulk_titration_summary.json"

MIX_SIZE = 1000
FRACTIONS = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.0]
N_DRAWS = 40
SEED = 20261115
CHUNK = 4000


def read_str(node):
    if isinstance(node, h5py.Dataset):
        return node.asstr()[:]
    if "categories" in node:
        cats = node["categories"].asstr()[:]
        codes = node["codes"][:]
        return np.array([cats[c] if c >= 0 else "" for c in codes])
    raise TypeError(node.name)


def load_counts():
    human = pd.read_csv(GENE_LIST)["gene"].tolist()
    human = list(dict.fromkeys(human))
    with h5py.File(H5AD, "r") as h:
        rX = h["raw"]["X"]
        data, indices, indptr = rX["data"], rX["indices"], rX["indptr"]
        n_cells = indptr.shape[0] - 1
        feat = read_str(h["raw"]["var"]["gene_symbol"])

        sym2cols, excluded = {}, []
        for g in human:
            mouse_sym = g[0] + g[1:].lower()
            cols = np.where(feat == mouse_sym)[0]
            if len(cols):
                sym2cols[g] = cols
            else:
                excluded.append(g)
        print("mapped genes:", len(sym2cols), "excluded:", excluded)

        wanted_cols = np.unique(np.concatenate(list(sym2cols.values())))
        col2local = {c: i for i, c in enumerate(wanted_cols)}
        total = np.zeros(n_cells)
        sub = np.zeros((n_cells, len(wanted_cols)))
        for s in range(0, n_cells, CHUNK):
            e = min(s + CHUNK, n_cells)
            a, b = indptr[s], indptr[e]
            d = data[a:b].astype(np.float64)
            ix = indices[a:b]
            rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
            np.add.at(total, rows, d)
            mask = np.isin(ix, wanted_cols)
            if mask.any():
                loc = np.array([col2local[c] for c in ix[mask]])
                np.add.at(sub, (rows[mask], loc), d[mask])
        ctype = read_str(h["obs"]["cell_type"])

    gene_mat = {
        g: sub[:, [col2local[c] for c in cols]].sum(axis=1)
        for g, cols in sym2cols.items()
    }
    return gene_mat, total, ctype, excluded


def main():
    gene_mat, total, ctype, excluded = load_counts()
    rng = np.random.default_rng(SEED)
    rows, summary = [], {}
    for gene, vals in gene_mat.items():
        means = pd.Series(vals).groupby(ctype).mean()
        producer = means.idxmax()
        if means.max() <= 0:
            continue
        prod_idx = np.where(ctype == producer)[0]
        other_idx = np.where(ctype != producer)[0]
        if len(prod_idx) < 5 or len(other_idx) < MIX_SIZE:
            continue
        rec = {"producer_celltype": producer, "points": []}
        cpm_by_f = {}
        for fraction in FRACTIONS:
            cpms = []
            n_p = int(round(MIX_SIZE * fraction))
            n_o = MIX_SIZE - n_p
            for _ in range(N_DRAWS):
                ps = rng.choice(prod_idx, n_p, replace=len(prod_idx) < n_p)
                oi = rng.choice(other_idx, n_o, replace=len(other_idx) < n_o)
                sel = np.concatenate([ps, oi]).astype(int)
                gs, ts = vals[sel].sum(), total[sel].sum()
                cpms.append(gs / ts * 1e6 if ts else np.nan)
            cpms = np.array(cpms)
            cpm_by_f[fraction] = cpms.mean()
            rows.append({
                "gene": gene, "producer_celltype": producer,
                "producer_fraction": fraction,
                "cpm_mean": round(float(cpms.mean()), 4),
                "cpm_sd": round(float(cpms.std()), 4),
            })
            rec["points"].append({
                "fraction": fraction, "cpm": round(float(cpms.mean()), 3)})
        full, bg = cpm_by_f[1.0], cpm_by_f[0.0]
        net = full - bg
        if net <= 0 or net < 0.2 * full:
            rec["net_signal_retained"] = None
            rec["specificity"] = "not_producer_specific"
            summary[gene] = rec
            continue
        rec["net_signal_retained"] = {
            str(f): round(float(np.clip((cpm_by_f[f] - bg) / net, 0, 1.5)), 3)
            for f in FRACTIONS}
        rec["specificity"] = "producer_specific"
        summary[gene] = rec

    pd.DataFrame(rows).to_csv(OUT, index=False)
    specific = [g for g in summary if summary[g]["net_signal_retained"] is not None]
    agg = {
        str(f): round(float(np.mean(
            [summary[g]["net_signal_retained"][str(f)] for g in specific])), 3)
        for f in FRACTIONS}
    out_json = {
        "collection": "mouse MASH LoxP/SIRT2 (MASH_liver.h5ad)",
        "n_cells": int(len(total)), "n_genes": len(summary),
        "n_producer_specific": len(specific),
        "excluded_no_ortholog": excluded,
        "not_specific": [g for g in summary
                         if summary[g]["net_signal_retained"] is None],
        "mix_size": MIX_SIZE, "n_draws": N_DRAWS,
        "mean_net_signal_retained": agg, "per_gene": summary,
    }
    OUT_JSON.write_text(json.dumps(out_json, ensure_ascii=False, indent=1), encoding="utf-8")
    print("genes analyzed:", len(summary), "producer-specific:", len(specific))
    print("mean net signal retained:", agg)


if __name__ == "__main__":
    main()
