# -*- coding: utf-8 -*-
"""Symmetric de-circularized benchmark across two formulas.

This is the reproducible version of the cross-formula analysis used in the
manuscript. It runs entirely on the small committed files in
data/processed_small and needs no live database access.

Inputs (data/processed_small):
  external_pain_gene_set.csv      literature-curated pain genes (reference labels)
  gse89632_significant_DEG.csv    significant disease DEG genes (selection rule)
  chaihu_filtered_gene_universe.csv   Chaihu Shugan San filtered-molecule genes
  jinling_filtered_gene_universe.csv  Jinlingzi San filtered-molecule genes

For each formula the disease axis and the pain reference are held fixed; only the
formula-derived gene universe changes. The selected set is the DEG intersection
and the positive labels are the external pain genes detectable in that universe.
A 10,000-draw permutation gives the null, alongside a one-sided Fisher test.

Output:
  cross_formula_benchmark.csv      metrics table (also stored in processed_small)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data" / "processed_small"
N_PERM = 10000
RNG = np.random.default_rng(20261116)

PAIN = set(pd.read_csv(DATA / "external_pain_gene_set.csv")["gene"])
DEG = set(pd.read_csv(DATA / "gse89632_significant_DEG.csv")["gene"])

UNIVERSES = {
    "Chaihu_Shugan_San": set(pd.read_csv(DATA / "chaihu_filtered_gene_universe.csv")["gene"]),
    "Jinlingzi_San": set(pd.read_csv(DATA / "jinling_filtered_gene_universe.csv")["gene"]),
}


def run_formula(name, universe):
    universe = set(universe)
    positive = PAIN & universe
    selected = DEG & universe
    k = len(selected)
    obs_tp = len(selected & positive)
    obs_prec = obs_tp / k if k else np.nan

    pool = np.array(sorted(universe))
    null_prec = []
    for _ in range(N_PERM):
        draw = set(RNG.choice(pool, size=k, replace=False))
        null_prec.append(len(draw & positive) / k)
    null_prec = np.array(null_prec)

    emp_p = (np.sum(null_prec >= obs_prec) + 1) / (N_PERM + 1)
    expected = k * len(positive) / len(universe)

    fp = k - obs_tp
    fn = len(positive) - obs_tp
    tn = len(universe) - obs_tp - fp - fn
    _, fisher_p = fisher_exact([[obs_tp, fp], [fn, tn]], alternative="greater")

    three_way = sorted(DEG & PAIN & universe)
    return {
        "formula": name,
        "universe_size": len(universe),
        "n_external_pain_in_universe": len(positive),
        "selected_DEG_intersect": k,
        "observed_pain_overlap": obs_tp,
        "observed_precision": round(obs_prec, 4),
        "expected_overlap_by_chance": round(expected, 2),
        "enrichment_vs_chance": round(obs_tp / expected, 2) if expected else None,
        "null_precision_mean": round(float(null_prec.mean()), 4),
        "null_precision_p95": round(float(np.percentile(null_prec, 95)), 4),
        "null_precision_max": round(float(null_prec.max()), 4),
        "fisher_p": fisher_p,
        "empirical_p": round(float(emp_p), 5),
        "three_way_genes": ";".join(three_way),
    }


results = pd.DataFrame([run_formula(name, uni) for name, uni in UNIVERSES.items()])
results.to_csv(REPO_ROOT / "scripts" / "13_cross_formula" / "cross_formula_benchmark.csv",
               index=False, encoding="utf-8-sig")
print(results.to_string(index=False))
