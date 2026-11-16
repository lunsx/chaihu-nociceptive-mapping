"""Symmetric de-circularized benchmark across two formulas.

The procedure is identical to scripts/13_cross_formula/01_cross_formula_benchmark.py
in the repository, but it reads the small reference files bundled with the
package so it runs in a clean install without network access. The disease axis
and the external pain reference are held fixed; only the formula-derived gene
universe changes. A 10,000-draw permutation gives the null, alongside a
one-sided Fisher test.
"""
from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

N_PERM = 10000
SEED = 20261116


def _load_genes(name):
    with files("decimap.data").joinpath(name).open("r", encoding="utf-8") as fh:
        return set(pd.read_csv(fh)["gene"])


def load_inputs():
    """Return the fixed pain labels, DEG rule, and the two gene universes."""
    pain = _load_genes("external_pain_gene_set.csv")
    deg = _load_genes("gse89632_significant_DEG.csv")
    universes = {
        "Chaihu_Shugan_San": _load_genes("chaihu_filtered_gene_universe.csv"),
        "Jinlingzi_San": _load_genes("jinling_filtered_gene_universe.csv"),
    }
    return pain, deg, universes


def _run_formula(name, universe, pain, deg, rng):
    universe = set(universe)
    positive = pain & universe
    selected = deg & universe
    k = len(selected)
    obs_tp = len(selected & positive)
    obs_prec = obs_tp / k if k else np.nan

    pool = np.array(sorted(universe))
    null_prec = []
    for _ in range(N_PERM):
        draw = set(rng.choice(pool, size=k, replace=False))
        null_prec.append(len(draw & positive) / k)
    null_prec = np.array(null_prec)

    emp_p = (np.sum(null_prec >= obs_prec) + 1) / (N_PERM + 1)
    expected = k * len(positive) / len(universe)

    fp = k - obs_tp
    fn = len(positive) - obs_tp
    tn = len(universe) - obs_tp - fp - fn
    _, fisher_p = fisher_exact([[obs_tp, fp], [fn, tn]], alternative="greater")

    three_way = sorted(deg & pain & universe)
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


def run_benchmark():
    """Run the symmetric benchmark and return a metrics DataFrame."""
    pain, deg, universes = load_inputs()
    rng = np.random.default_rng(SEED)
    rows = []
    for name, universe in universes.items():
        rows.append(_run_formula(name, universe, pain, deg, rng))
    return pd.DataFrame(rows)
