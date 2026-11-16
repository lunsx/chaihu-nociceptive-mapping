"""Simulation-based validation of the degree-aware estimator.

The estimator cannot be validated against a wet-lab ground truth here, so it is
validated on synthetic compound-gene graphs with known truth. Each replicate
builds a strongly skewed degree distribution (a few mega-hubs, a long tail) and
applies two tests: a uniform (exchangeable) permutation and the degree-weighted
Pareto permutation implemented in degree_aware.py.

Scenarios:
  S0 exchangeable null  uniform selected and reference, no association.
  S1 hub-confounded null  both sets independently favour hubs; the apparent
    overlap is induced by degree, not by a true association.
  S2 alternative  low-degree reference genes are planted into the selected set,
    representing genuine, non-hub signal.

The output table reports rejection rates at alpha = 0.05, i.e. empirical type I
error for S0/S1 and power for S2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

N_GENES = 250
K_SEL = 50
M_POS = 12
N_REPS = 500
N_NULL = 1000
ALPHA = 0.05
SEED = 20261116


def make_degrees(n_genes=N_GENES, rng=None):
    """Deterministic-shape skewed degree vector with shuffled gene positions."""
    rng = rng or np.random.default_rng(SEED)
    degrees = np.ones(n_genes, dtype=float)
    degrees[0:6] = rng.integers(80, 110, 6)
    degrees[6:26] = rng.integers(12, 35, 20)
    degrees[26:60] = rng.integers(2, 8, 34)
    rng.shuffle(degrees)
    return degrees


def pareto_draw(weights, size, rng):
    u = rng.random(len(weights))
    key = u ** (1.0 / np.maximum(weights, 1.0))
    return np.argpartition(key, len(weights) - size)[-size:]


def uniform_draw(n, size, rng):
    return rng.choice(n, size=size, replace=False)


def permutation_p(degrees, selected, positive, weighted, n_null, rng):
    pos_mask = np.zeros(len(degrees), dtype=bool)
    pos_mask[positive] = True
    observed = int(pos_mask[selected].sum())
    null_hits = np.zeros(n_null, dtype=np.int64)
    for r in range(n_null):
        if weighted:
            draw = pareto_draw(degrees, len(selected), rng)
        else:
            draw = uniform_draw(len(degrees), len(selected), rng)
        null_hits[r] = int(pos_mask[draw].sum())
    return (np.sum(null_hits >= observed) + 1) / (n_null + 1)


def _replicate(degrees, confounded, n_plant, n_null, rng):
    if n_plant:
        positive = pareto_draw(degrees, M_POS, rng)
        selected = pareto_draw(degrees, K_SEL, rng)
        order = np.argsort(degrees[positive])
        low_positive = positive[order][:n_plant]
        selected = np.unique(np.concatenate([selected, low_positive]))[:K_SEL]
    elif confounded:
        selected = pareto_draw(degrees, K_SEL, rng)
        positive = pareto_draw(degrees, M_POS, rng)
    else:
        selected = uniform_draw(len(degrees), K_SEL, rng)
        positive = uniform_draw(len(degrees), M_POS, rng)
    uniform_p = permutation_p(degrees, selected, positive, False, n_null, rng)
    degree_p = permutation_p(degrees, selected, positive, True, n_null, rng)
    return uniform_p < ALPHA, degree_p < ALPHA


def run_simulation(n_reps=N_REPS, n_null=N_NULL, seed=SEED):
    """Return a tidy DataFrame of empirical rejection rates by scenario."""
    rng = np.random.default_rng(seed)
    rows = []

    for scenario, confounded, n_plant in [
            ("S0_exchangeable_null", False, 0),
            ("S1_hub_confounded_null", True, 0)]:
        u, d = [], []
        for _ in range(n_reps):
            degrees = make_degrees(rng=rng)
            a, b = _replicate(degrees, confounded, n_plant, n_null, rng)
            u.append(a); d.append(b)
        rows.append((scenario, "uniform_permutation", round(float(np.mean(u)), 4)))
        rows.append((scenario, "degree_aware_permutation", round(float(np.mean(d)), 4)))

    for n_plant in (2, 4, 6, 8):
        d = []
        for _ in range(n_reps):
            degrees = make_degrees(rng=rng)
            _, b = _replicate(degrees, True, n_plant, n_null, rng)
            d.append(b)
        rows.append((f"S2_alternative_plant{n_plant}", "degree_aware_permutation",
                      round(float(np.mean(d)), 4)))

    return pd.DataFrame(rows, columns=["scenario", "test", "rejection_rate_alpha05"])
