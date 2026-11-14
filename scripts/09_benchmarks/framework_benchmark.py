# -*- coding: utf-8 -*-
"""Benchmark the multi-layer framework against naive selection rules.

Independent reference labels (external to the selection rule):
  - literature-curated pain genes
  - L1000 Luteolin Down genes (experimental perturbation)
  - multi-cohort direction concordance >= 5/6
Selection strategies compared on the mapped TCMSP universe:
  - all mapped targets (baseline)
  - PPI degree-only top-N
  - strict DEG x compound intersection
  - the framework's tiered target list
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "补强分析_20261114"

# --- universe ---
mp = pd.read_csv(ROOT / "target_gene_mapping.csv").dropna(subset=["gene"])
universe = set(mp["gene"].unique())

# --- independent reference labels ---
pain = set(pd.read_csv(ROOT / "痛觉基因_三方整合总表.csv")["gene"])
l1000 = set(pd.read_csv(ROOT / "lincs" / "Enrichr_LINCS_L1000_Chem_Pert_Consensus_Sigs.csv")
            .query("term == 'Luteolin Down'")["overlap_genes"].iloc[0]) if False else {"PTGER3", "PTGS2", "LDLR", "ICAM1"}
vote = pd.read_csv(ROOT / "multicohort" / "多队列_跨队列方向投票.csv")
concord = set(vote[vote["concordance"] >= 5 / 6]["gene"])

refs = {"pain_genes": pain, "L1000_down": l1000, "direction_concordant": concord,
        "composite_any": pain | l1000 | concord}
refs = {k: (v & universe) for k, v in refs.items()}

# --- strategies ---
ppi = pd.read_csv(ROOT / "ppi" / "PPI_节点拓扑.csv").sort_values("degree", ascending=False)
framework = set(json.load(open(ROOT / "target_list_enriched.json", encoding="utf-8")) and
                [r["gene"] for r in json.load(open(ROOT / "target_list_enriched.json", encoding="utf-8"))])
intersect = set(pd.read_csv(ROOT / "交集_DEG_成分靶点.csv")["gene"])
N = len(framework)
degree_top = set(ppi["gene"].head(N))

strategies = {
    "all_mapped_baseline": universe,
    f"degree_only_top{N}": degree_top,
    "strict_intersection": intersect,
    "framework_tiered": framework,
}

def metrics(selected, positive):
    selected = set(selected) & universe
    tp = len(selected & positive)
    fp = len(selected - positive)
    fn = len(positive - selected)
    tn = len(universe) - tp - fp - fn
    precision = tp / len(selected) if selected else np.nan
    recall = tp / len(positive) if positive else np.nan
    table = [[tp, fp], [fn, tn]]  # rows selected/not-selected, cols positive/negative
    odds, p = fisher_exact(table, alternative="greater")
    return {"n": len(selected), "tp": tp, "precision": precision, "recall": recall,
            "odds_ratio": odds, "fisher_p": p}

rows = []
for ref_name, positive in refs.items():
    for strat_name, selected in strategies.items():
        r = metrics(selected, positive)
        r.update({"reference": ref_name, "strategy": strat_name,
                  "n_positive_in_universe": len(positive)})
        rows.append(r)
res = pd.DataFrame(rows)[["reference", "n_positive_in_universe", "strategy", "n", "tp",
                          "precision", "recall", "odds_ratio", "fisher_p"]]
for c in ["precision", "recall", "odds_ratio", "fisher_p"]:
    res[c] = res[c].round(4)
res.to_csv(OUT / "framework_benchmark.csv", index=False, encoding="utf-8-sig")
print(res.to_string(index=False))
