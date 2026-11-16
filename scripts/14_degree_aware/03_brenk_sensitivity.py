# -*- coding: utf-8 -*-
"""Brenk named-alert exclusion sensitivity for the degree-aware enrichment.

Each scenario excludes the compounds carrying one named Brenk alert class (or a
set of reactive toxicophores, or every flagged compound), rebuilds the
compound-gene graph from the retained compounds, and reruns the degree-aware
enrichment. Writes data/processed_small/brenk_sensitivity.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

sys.path.insert(0, str(PATHS["repo_root"] / "src"))
from decimap.degree_aware import degree_aware_enrichment  # noqa: E402

small = PATHS["processed_small"]
alerts = pd.read_csv(small / "compound_filter_alerts.csv")
edges = pd.read_csv(small / "chaihu_compound_gene_edges.csv")
pain = set(pd.read_csv(small / "external_pain_gene_set.csv")["gene"])
deg = set(pd.read_csv(small / "gse89632_significant_DEG.csv")["gene"])

REACTIVE = {"Michael_acceptor_1", "chinone_1", "quinone_A(370)",
            "imine_one_A(321)", "diketo_group"}
SCENARIOS = {
    "baseline": set(),
    "reactive_excluded": REACTIVE,
    "isolated_alkene_excluded": {"isolated_alkene"},
    "catechol_excluded": {"catechol", "catechol_A(92)"},
    "cumarine_excluded": {"cumarine"},
    "all_brenk_excluded": None,  # every non-empty alert
}


def alert_types(raw):
    if not isinstance(raw, str) or not raw.strip():
        return set()
    return {t.strip() for t in raw.split(";") if t.strip()}


def run(exclude):
    if exclude is None:
        keep = {r["MOL_ID"] for _, r in alerts.iterrows()
                if not alert_types(r["brenk_alert"])}
    else:
        keep = {r["MOL_ID"] for _, r in alerts.iterrows()
                if not (alert_types(r["brenk_alert"]) & exclude)}
    sub = edges[edges["MOL_ID"].isin(keep)]
    universe = set(sub["gene"])
    degree = sub.groupby("gene").size().to_dict()
    rng = np.random.default_rng(20261116)
    result = degree_aware_enrichment(
        universe, deg & universe, pain, degree, rng=rng)
    return {
        "compounds_kept": len(keep),
        "universe_size": result["universe_size"],
        "positive_in_universe": result["positive_in_universe"],
        "observed_overlap": result["observed_overlap"],
        "corrected_enrichment": result["corrected_enrichment"],
        "degree_weighted_empirical_p": result["degree_weighted_empirical_p"],
    }


rows = []
for name, exclude in SCENARIOS.items():
    row = {"scenario": name, **run(exclude)}
    rows.append(row)

table = pd.DataFrame(rows)
out = small / "brenk_sensitivity.csv"
table.to_csv(out, index=False, encoding="utf-8")
print(table.to_string(index=False))
print("saved:", out)
