# -*- coding: utf-8 -*-
"""Independent method comparison and prediction audit (pure computation).

Adds classical exchangeable-null tests (exact hypergeometric, Fisher), an
exploratory rank-based check on meta Hedges g, and a per-cohort confusion
audit with exact/Wilson intervals and a bootstrap. Reads only committed
files and writes new additive tables; no frozen data are modified.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import (beta, fisher_exact, hypergeom, mannwhitneyu)

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT = PATHS["processed_small"]
SEED = 20261119


def wilson_ci(tp_correct, n, conf=0.95):
    if n == 0:
        return float("nan"), float("nan")
    z = 1.959964
    p = tp_correct / n
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    denom = 1 + z * z / n
    lo = (p + z * z / (2 * n) - half) / denom
    hi = (p + z * z / (2 * n) + half) / denom
    return lo, hi


def enrichment_comparison():
    da = pd.read_csv(OUT / "degree_aware_enrichment.csv")
    rows = []
    for _, r in da.iterrows():
        n, k, m, a = int(r.universe_size), int(r.selected_size), int(
            r.positive_in_universe), int(r.observed_overlap)
        b, c, d = k - a, m - a, n - k - m + a
        hyper_p = hypergeom.sf(a - 1, n, m, k)
        odds, fish_p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "formula": r.formula,
            "universe_n": n, "selected_k": k, "pain_m": m, "overlap_a": a,
            "naive_fold": r.naive_enrichment,
            "degree_aware_corrected_fold": r.corrected_enrichment,
            "hypergeom_p": round(float(hyper_p), 6),
            "fisher_odds_greater": round(float(odds), 4),
            "fisher_p": round(float(fish_p), 6),
            "degree_aware_p": r.degree_weighted_empirical_p,
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "method_comparison_enrichment.csv", index=False)
    return frame


def rank_based_check():
    meta = pd.read_csv(OUT / "Table5_Meta效应BH.csv")
    pain = set(pd.read_csv(Path(__file__).resolve().parents[2]
              / "src" / "decimap" / "data" / "external_pain_gene_set.csv")["gene"])
    meta["is_pain"] = meta["gene"].isin(pain)
    g_pain = meta.loc[meta.is_pain, "pooled_hedges_g"].to_numpy()
    g_other = meta.loc[~meta.is_pain, "pooled_hedges_g"].to_numpy()
    u, p = mannwhitneyu(g_pain, g_other, alternative="two-sided")
    res = {
        "n_genes": int(len(meta)),
        "n_pain_in_list": int(meta.is_pain.sum()),
        "pain_genes": sorted(meta.loc[meta.is_pain, "gene"]),
        "median_g_pain": round(float(np.median(g_pain)), 4),
        "median_g_other": round(float(np.median(g_other)), 4),
        "mean_g_pain": round(float(np.mean(g_pain)), 4),
        "mean_g_other": round(float(np.mean(g_other)), 4),
        "mannwhitney_p_two_sided": round(float(p), 6),
        "note": "Exploratory rank-based check on pooled Hedges g; no "
                "directional enrichment across the 54-gene list is assumed.",
    }
    with open(OUT / "rank_based_check.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    return res


def confusion_for(detail, cohort, method_col):
    sub = detail[detail.cohort == cohort]
    tp = tn = fp = fn = excluded = 0
    for _, r in sub.iterrows():
        obs, pred = r.observed, r[method_col]
        if obs == "zero" or pd.isna(pred) or str(pred).strip() == "":
            excluded += 1
            continue
        if obs == "up" and pred == "up":
            tp += 1
        elif obs == "down" and pred == "down":
            tn += 1
        elif obs == "down" and pred == "up":
            fp += 1
        elif obs == "up" and pred == "down":
            fn += 1
        else:
            excluded += 1
    n = tp + tn + fp + fn
    correct = tp + tn
    acc = correct / n if n else float("nan")
    if n:
        if correct == 0:
            cp_lo = 0.0
        else:
            cp_lo = beta.ppf(0.025, correct, n - correct + 1)
        if correct == n:
            cp_hi = 1.0
        else:
            cp_hi = beta.ppf(0.975, correct + 1, n - correct)
    else:
        cp_lo, cp_hi = float("nan"), float("nan")
    w_lo, w_hi = wilson_ci(correct, n)
    up_baseline = int((sub.observed == "up").sum())
    return {
        "cohort": cohort, "method": method_col, "n_scored": n,
        "excluded": excluded, "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "accuracy": round(acc, 4),
        "clopper_pearson_lo": round(float(cp_lo), 4),
        "clopper_pearson_hi": round(float(cp_hi), 4),
        "wilson_lo": round(float(w_lo), 4),
        "wilson_hi": round(float(w_hi), 4),
        "always_up_accuracy": round(up_baseline / len(sub), 4),
    }


def prediction_audit():
    detail = pd.read_csv(OUT / "heldout_prediction_detail.csv")
    cohorts = sorted(detail.cohort.unique())
    rows = []
    for cohort in cohorts:
        rows.append(confusion_for(detail, cohort, "prior_pred"))
        rows.append(confusion_for(detail, cohort, "loco_pred"))
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "prediction_audit.csv", index=False)

    # Bootstrap CI for the fully held-out cohort prior accuracy.
    held = detail[detail.cohort == "GSE147304"]
    flags = []
    for _, r in held.iterrows():
        obs, pred = r.observed, r.prior_pred
        if obs == "zero" or pd.isna(pred) or str(pred).strip() == "":
            continue
        flags.append(obs == pred)
    flags = np.array(flags, dtype=bool)
    rng = np.random.default_rng(SEED)
    boot = [rng.choice(flags, size=len(flags), replace=True).mean()
            for _ in range(2000)]
    ci = {
        "heldout_cohort": "GSE147304",
        "n_scored": int(len(flags)),
        "point_accuracy": round(float(flags.mean()), 4),
        "bootstrap_95ci_lo": round(float(np.percentile(boot, 2.5)), 4),
        "bootstrap_95ci_hi": round(float(np.percentile(boot, 97.5)), 4),
    }
    with open(OUT / "prediction_audit_ci.json", "w", encoding="utf-8") as fh:
        json.dump(ci, fh, indent=2)
    return frame, ci


def main():
    print("== enrichment comparison ==")
    print(enrichment_comparison().to_string(index=False))
    print("\n== rank-based check ==")
    print(json.dumps(rank_based_check(), indent=2))
    print("\n== prediction audit ==")
    frame, ci = prediction_audit()
    print(frame.to_string(index=False))
    print(json.dumps(ci, indent=2))


if __name__ == "__main__":
    main()
