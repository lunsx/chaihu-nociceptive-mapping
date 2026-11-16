"""Degree-aware de-circularized enrichment inference.

Classical enrichment tests treat every target gene as exchangeable. In a TCMSP
compound-gene graph that assumption is false: a few hub proteins are mapped by
most compounds, and curated reference sets are enriched for those same hubs. A
uniform permutation therefore overstates the evidence.

The primary estimator here replaces the uniform null with a fixed-size
probability-proportional-to-degree sample. A Pareto (order-sampling) design draws
u_g ~ Uniform(0, 1), forms the key k_g = u_g ** (1 / degree_g), and selects the
largest k keys, so high-degree genes enter the null more often. Monte Carlo gives
each gene's marginal inclusion probability and the null overlap distribution.
The bias-corrected enrichment divides the observed overlap by the degree-weighted
expectation rather than the uniform expectation. The module also provides:

  * empirical calibration from degree-matched negative-control reference sets;
  * a conditional portability test that measures shared three-way convergence
    within the universe intersection, beyond what label independence predicts.

Only numpy, scipy and pandas are required.
"""
from __future__ import annotations

from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.stats import norm

N_NULL = 20000
N_BOOT = 500
N_CONTROLS = 500
SEED = 20261116


def _load_genes(name):
    with files("decimap.data").joinpath(name).open("r", encoding="utf-8") as fh:
        return set(pd.read_csv(fh)["gene"])


def _load_degree(name):
    with files("decimap.data").joinpath(name).open("r", encoding="utf-8") as fh:
        frame = pd.read_csv(fh)
    return dict(zip(frame["gene"], frame["degree"]))


def run_degree_analysis():
    """Degree-aware enrichment for both formulas from bundled data."""
    pain = _load_genes("external_pain_gene_set.csv")
    deg = _load_genes("gse89632_significant_DEG.csv")
    configs = {
        "Chaihu_Shugan_San": ("chaihu_filtered_gene_universe.csv",
                              "chaihu_gene_degree.csv"),
        "Jinlingzi_San": ("jinling_filtered_gene_universe.csv",
                          "jinling_gene_degree.csv"),
    }
    rng = np.random.default_rng(SEED)
    rows = []
    for formula, (universe_file, degree_file) in configs.items():
        universe = _load_genes(universe_file)
        degree = _load_degree(degree_file)
        result = degree_aware_enrichment(
            universe, deg & universe, pain, degree, rng=rng)
        result.pop("inclusion_probabilities", None)
        result = {"formula": formula, **result}
        rows.append(result)
    return pd.DataFrame(rows)


def run_calibration():
    """Empirical calibration for the first formula from bundled data."""
    pain = _load_genes("external_pain_gene_set.csv")
    deg = _load_genes("gse89632_significant_DEG.csv")
    universe = _load_genes("chaihu_filtered_gene_universe.csv")
    degree = _load_degree("chaihu_gene_degree.csv")
    rng = np.random.default_rng(SEED)
    return pd.DataFrame([empirical_calibration(
        universe, deg & universe, pain, degree, rng=rng)])


def run_portability():
    """Conditional portability across both formulas from bundled data."""
    deg = _load_genes("gse89632_significant_DEG.csv")
    pain = _load_genes("external_pain_gene_set.csv")
    u1 = _load_genes("chaihu_filtered_gene_universe.csv")
    u2 = _load_genes("jinling_filtered_gene_universe.csv")
    rng = np.random.default_rng(SEED)
    return pd.DataFrame([conditional_portability(
        u1, u2, deg & u1, deg & u2, pain, rng=rng)])


def _as_set(values):
    return set(values)


def pareto_null(universe, weights, k, positive_mask, n_draw=N_NULL, rng=None):
    """Fixed-size probability-proportional-to-degree null.

    Returns marginal inclusion probabilities, the null overlap array and the
    degree-weighted expected overlap.
    """
    rng = rng or np.random.default_rng(SEED)
    n = len(universe)
    weights = np.asarray(weights, dtype=float)
    pi = np.zeros(n)
    null_overlap = np.zeros(n_draw, dtype=np.int64)
    for r in range(n_draw):
        u = rng.random(n)
        key = u ** (1.0 / weights)
        chosen = np.argpartition(key, n - k)[-k:]
        pi[chosen] += 1
        null_overlap[r] = int(positive_mask[chosen].sum())
    pi /= n_draw
    return pi, null_overlap, float(pi[positive_mask].sum())


def degree_aware_enrichment(universe, selected, positive, degree,
                            n_draw=N_NULL, n_boot=N_BOOT, rng=None):
    """Naive and degree-aware enrichment statistics for one formula."""
    rng = rng or np.random.default_rng(SEED)
    universe = sorted(_as_set(universe))
    selected = _as_set(selected) & set(universe)
    positive = _as_set(positive) & set(universe)
    weights = np.array([degree.get(g, 1) for g in universe], dtype=float)
    sel_mask = np.isin(universe, list(selected))
    pos_mask = np.isin(universe, list(positive))
    n, k, m = len(universe), int(sel_mask.sum()), int(pos_mask.sum())
    obs = int((sel_mask & pos_mask).sum())

    naive_expected = k * m / n if n else np.nan
    pi, null_overlap, weighted_expected = pareto_null(
        universe, weights, k, pos_mask, n_draw, rng)
    emp_p = (np.sum(null_overlap >= obs) + 1) / (n_draw + 1)

    # Degree-perturbation sensitivity interval: the committed degrees are exact,
    # so this measures how the corrected ratio moves when the degree signal is
    # perturbed by multiplicative noise, rather than a sampling CI.
    boot_denominator = []
    for _ in range(n_boot):
        perturbed = weights * rng.normal(1.0, 0.1, size=n)
        perturbed = np.clip(perturbed, 1.0, None)
        pi_b = np.zeros(n)
        for r in range(40):
            u = rng.random(n)
            key = u ** (1.0 / perturbed)
            pi_b[np.argpartition(key, n - k)[-k:]] += 1
        pi_b /= 40
        boot_denominator.append(pi_b[pos_mask].sum())
    boot_ratio = obs / np.array(boot_denominator)
    sens_low, sens_high = np.percentile(boot_ratio, [2.5, 97.5])

    inclusion = {g: round(float(pi[i]), 4) for i, g in enumerate(universe)}
    return {
        "universe_size": n,
        "selected_size": k,
        "positive_in_universe": m,
        "observed_overlap": obs,
        "naive_expected": round(naive_expected, 4),
        "naive_enrichment": round(obs / naive_expected, 4) if naive_expected else np.nan,
        "degree_weighted_expected": round(weighted_expected, 4),
        "corrected_enrichment": round(obs / weighted_expected, 4) if weighted_expected else np.nan,
        "sensitivity_low": round(float(sens_low), 4),
        "sensitivity_high": round(float(sens_high), 4),
        "null_p95": round(float(np.percentile(null_overlap, 95)), 4),
        "null_max": int(null_overlap.max()),
        "degree_weighted_empirical_p": round(float(emp_p), 5),
        "degree_inflation_expected": round(weighted_expected - naive_expected, 4),
        "inclusion_probabilities": inclusion,
    }


def empirical_calibration(universe, selected, positive, degree,
                          n_controls=N_CONTROLS, n_null=2000, rng=None):
    """Calibrate the degree-aware statistic with degree-matched control sets."""
    rng = rng or np.random.default_rng(SEED)
    universe = sorted(_as_set(universe))
    selected = _as_set(selected) & set(universe)
    positive = _as_set(positive) & set(universe)
    weights = np.array([degree.get(g, 1) for g in universe], dtype=float)
    sel_mask = np.isin(universe, list(selected))
    k = int(sel_mask.sum())

    def statistic(pos_set):
        pos_m = np.isin(universe, list(pos_set))
        pi, null_t, expected = pareto_null(
            universe, weights, k, pos_m, n_null, rng)
        obs_here = int((sel_mask & pos_m).sum())
        if not expected or not null_t.var():
            return np.nan
        # standardized against the empirical null
        return (obs_here - null_t.mean()) / null_t.std(ddof=1)

    observed_z = statistic(positive)
    control_z = []
    positive_arr = np.isin(universe, list(positive))
    n_positive = int(positive_arr.sum())
    eligible_index = np.flatnonzero(~positive_arr)
    for _ in range(n_controls):
        control_idx = rng.choice(eligible_index,
                                 size=min(n_positive, len(eligible_index)),
                                 replace=False)
        control_set = set(np.array(universe)[control_idx])
        z = statistic(control_set)
        if not np.isnan(z):
            control_z.append(z)
    control_z = np.array(control_z)

    mu0 = float(np.mean(control_z))
    sigma0 = float(np.std(control_z, ddof=1))
    robust_sigma = 1.4826 * float(np.median(np.abs(control_z - np.median(control_z))))
    calibrated_z = (observed_z - mu0) / sigma0 if sigma0 else np.nan
    return {
        "observed_degree_z": round(float(observed_z), 4),
        "n_negative_control_sets": len(control_z),
        "empirical_null_mean": round(mu0, 4),
        "empirical_null_sd": round(sigma0, 4),
        "empirical_null_robust_sd": round(robust_sigma, 4),
        "calibrated_z": round(float(calibrated_z), 4),
        "calibrated_p": round(float(norm.sf(calibrated_z)), 6) if not np.isnan(calibrated_z) else np.nan,
        "empirical_calibrated_p": round(float(np.mean(control_z >= observed_z)), 5),
    }


def conditional_portability(universe1, universe2, selected1, selected2,
                            positive, n_boot=N_BOOT, rng=None):
    """Shared three-way convergence within the universe intersection."""
    rng = rng or np.random.default_rng(SEED)
    u1, u2 = _as_set(universe1), _as_set(universe2)
    s1, s2 = _as_set(selected1) & u1, _as_set(selected2) & u2
    positive = _as_set(positive)
    inter = u1 & u2
    union = u1 | u2
    jaccard = len(inter) / len(union) if union else np.nan

    # DEG labels are one shared external rule, so selected genes in the
    # intersection coincide; report that count for context, then test the
    # three-way (DEG and pain) enrichment within the intersection.
    deg_i = len((s1 | s2) & inter)
    pain_i = len(positive & inter)
    obs_threeway = len(s1 & s2 & positive & inter)
    expected_threeway = deg_i * pain_i / len(inter) if inter else np.nan
    variance = (deg_i * (pain_i / len(inter)) * (1 - pain_i / len(inter))
                * (len(inter) - deg_i) / (len(inter) - 1)) if len(inter) > 1 else 0.0
    z = ((obs_threeway - expected_threeway) / np.sqrt(variance)) if variance else np.nan

    boot_threeway = rng.hypergeometric(pain_i, len(inter) - pain_i, deg_i,
                                        size=n_boot)
    null_excess = boot_threeway - expected_threeway
    null_p2_5, null_p97_5 = np.percentile(null_excess, [2.5, 97.5])
    emp_p = float(np.mean(boot_threeway >= obs_threeway))

    return {
        "jaccard_universe": round(jaccard, 4),
        "intersection_size": len(inter),
        "deg_genes_in_intersection": deg_i,
        "pain_genes_in_intersection": pain_i,
        "observed_shared_threeway": obs_threeway,
        "expected_shared_threeway": round(expected_threeway, 4),
        "excess_shared_threeway": round(obs_threeway - expected_threeway, 4),
        "null_excess_p2.5": round(float(null_p2_5), 4),
        "null_excess_p97.5": round(float(null_p97_5), 4),
        "threeway_z": round(float(z), 4),
        "threeway_empirical_p": round(emp_p, 5),
    }
