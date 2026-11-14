# -*- coding: utf-8 -*-
"""Quantitative framing of the null random-effects meta-analysis.

Adds, for every pooled gene: 95% confidence interval, 95% prediction interval,
the upper bound of |effect| still compatible with the data, and a two one-sided
equivalence test (TOST) against a pre-specified small-effect margin.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "补强分析_20261114"
Z = norm.ppf(0.975)

df = pd.read_csv(ROOT / "multicohort" / "Meta_随机效应合并结果.csv")

def tost_p(g, se, margin):
    # H01: g <= -margin (test g > -margin); H02: g >= margin (test g < margin)
    p_low = 1 - norm.cdf((g + margin) / se)
    p_high = norm.cdf((g - margin) / se)
    return max(p_low, p_high)

rows = []
for _, r in df.iterrows():
    g, se, tau2 = r["pooled_hedges_g"], r["se"], r["tau2"]
    ci_lo, ci_hi = g - Z * se, g + Z * se
    pi_lo, pi_hi = g - Z * np.sqrt(se ** 2 + tau2), g + Z * np.sqrt(se ** 2 + tau2)
    rows.append({
        "gene": r["gene"], "k": int(r["k"]), "pooled_g": g, "se": se,
        "ci95_low": ci_lo, "ci95_high": ci_hi,
        "pi95_low": pi_lo, "pi95_high": pi_hi,
        "abs_effect_upper_bound": max(abs(ci_lo), abs(ci_hi)),
        "tost_p_margin0.2": tost_p(g, se, 0.2),
        "tost_p_margin0.3": tost_p(g, se, 0.3),
        "p": r["p"], "I2_percent": r["I2_percent"],
    })
out = pd.DataFrame(rows)
out["equivalent_0.2"] = out["tost_p_margin0.2"] < 0.05
out["equivalent_0.3"] = out["tost_p_margin0.3"] < 0.05
out = out.sort_values("p")
out.to_csv(OUT / "meta_ci_table.csv", index=False, encoding="utf-8-sig")

core = ["PTGS2", "TNF", "IL1B", "IL6", "TACR1", "CALCRL", "TRPV1"]
cols = ["gene", "pooled_g", "ci95_low", "ci95_high", "abs_effect_upper_bound",
        "pi95_low", "pi95_high", "tost_p_margin0.2", "tost_p_margin0.3"]
print(out[out.gene.isin(core)][cols].round(3).to_string(index=False))
print("\nTOST equivalence counts (margin 0.2):", int(out["equivalent_0.2"].sum()), "of", len(out))
print("TOST equivalence counts (margin 0.3):", int(out["equivalent_0.3"].sum()), "of", len(out))
print("median upper bound of |g|:", round(out["abs_effect_upper_bound"].median(), 3))
