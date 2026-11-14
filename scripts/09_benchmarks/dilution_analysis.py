# -*- coding: utf-8 -*-
"""Bulk dilution analysis.

Question: why are the bulk meta effects (Hedges g) small even though single-cell
data show clear disease changes? Test whether the answer is dilution by the
fraction of cells that express each gene.

All inputs are local public-data re-analyses. The human single-cell disease is
primary sclerosing cholangitis (PSC), an inflammatory liver condition; the
dilution mechanism itself is a generic property of bulk averaging and does not
depend on the disease label.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "补强分析_20261114"

# --- single-cell whole-cell-level disease comparison (cell_type == ALL) ---
sc = pd.read_csv(ROOT / "human_sc" / "人源单细胞_疾病对比.csv")
sc_all = sc[sc["cell_type"] == "ALL"].copy()
# keep one row per gene (the inflammatory disease present in the file)
sc_all = sc_all.sort_values("pct_disease", ascending=False).drop_duplicates("gene")

# --- bulk meta pooled results ---
meta = pd.read_csv(ROOT / "multicohort" / "Meta_随机效应合并结果.csv")

# --- per-cohort effects (large well-powered cohort) ---
per = pd.read_csv(ROOT / "multicohort" / "Meta_各队列效应量.csv")
big = per[per["cohort"] == "GSE281797"][["gene", "hedges_g"]].rename(columns={"hedges_g": "g_big_cohort"})

# --- deterministic within-gene mixing from the localization table ---
loc = pd.read_csv(ROOT / "human_sc" / "人源单细胞_细胞类型定位.csv")
def weighted(gdf, val):
    w = gdf["n_cells"]
    return (gdf[val] * w).sum() / w.sum()
rows = []
for gene, gdf in loc.groupby("gene"):
    p = weighted(gdf, "pct_expressing") / 100.0
    m_when = weighted(gdf, "mean_when_expressed")
    bulk_mean = (gdf["pct_expressing"] / 100.0 * gdf["mean_when_expressed"] * gdf["n_cells"]).sum() / gdf["n_cells"].sum()
    rows.append({"gene": gene, "pct_overall": p * 100, "mean_when_expressed": m_when,
                 "zero_inclusive_bulk_mean": bulk_mean,
                 "within_gene_attenuation": (m_when / bulk_mean) if bulk_mean > 0 else np.nan})
mix = pd.DataFrame(rows)

# --- merge ---
m = sc_all.merge(meta[["gene", "pooled_hedges_g", "p", "I2_percent"]], on="gene", how="inner")
m = m.merge(big, on="gene", how="left")
m = m.merge(mix, on="gene", how="left")

# dilution factor: a signal confined to fraction p of cells is attenuated ~1/p
m["p_fraction"] = m["pct_disease"] / 100.0
m["dilution_factor"] = 1.0 / m["p_fraction"].clip(lower=0.001)
m["sc_bulk_equivalent_delta"] = m["mean_delta"]

# --- correlations across overlapping genes ---
def corr(x, y, method):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    if method == "spearman":
        r, p = spearmanr(d.x, d.y)
    else:
        r, p = pearsonr(d.x, d.y)
    return float(r), float(p), len(d)

stats = {
    "n_overlap_genes": int(m.shape[0]),
    "median_pct_disease_percent": float(m["pct_disease"].median()),
    "median_abs_pooled_g": float(m["pooled_hedges_g"].abs().median()),
    "median_within_gene_attenuation": float(m["within_gene_attenuation"].median()),
    # deterministic mechanism: within the same single-cell dataset, the zero-inclusive
    # mean change should scale with the expressing fraction
    "pct_vs_sc_delta_same_dataset_spearman": corr(m["pct_disease"], m["sc_bulk_equivalent_delta"], "spearman"),
    # honest cross-gene test against the independent bulk meta (may be null because
    # bulk g are noisy and near zero)
    "pct_vs_abs_bulk_g_spearman": corr(m["pct_disease"], m["pooled_hedges_g"].abs(), "spearman"),
    "pct_vs_abs_g_big_cohort_spearman": corr(m["pct_disease"], m["g_big_cohort"].abs(), "spearman"),
}

m_out = m[["gene", "pct_disease", "pct_normal", "dilution_factor", "mean_when_expressed",
           "zero_inclusive_bulk_mean", "within_gene_attenuation", "mean_delta",
           "pooled_hedges_g", "p", "g_big_cohort", "I2_percent"]].sort_values("pct_disease", ascending=False)
m_out.to_csv(OUT / "dilution_gene_table.csv", index=False, encoding="utf-8-sig")
with open(OUT / "dilution_stats.json", "w", encoding="utf-8") as f:
    json.dump(stats, f, ensure_ascii=False, indent=2)

print(json.dumps(stats, ensure_ascii=False, indent=2))
print("\ncore genes:")
core = ["PTGS2", "IL1B", "IL6", "TNF", "TRPV1", "PTGES", "PTGER3", "IL1A"]
print(m_out[m_out.gene.isin(core)].to_string(index=False))
