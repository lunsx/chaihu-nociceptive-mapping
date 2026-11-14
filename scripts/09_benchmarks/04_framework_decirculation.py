# -*- coding: utf-8 -*-
"""De-circularized re-evaluation of the framework benchmark.

The original benchmark (framework_benchmark.py) used a composite reference that
included L1000 Luteolin-Down genes and multi-cohort direction concordance — both
of which are also layers the framework itself uses for tiering. That makes part
of the precision/recall self-fulfilling.

This script removes the overlap in four ways:
  1. primary reference = literature-curated pain genes only (no shared provenance)
  2. leave-one-source-out references (drop pain / L1000 / direction in turn)
  3. permutation null: 10,000 random 60-gene draws from the mapped universe
  4. degree-matched control: highest-degree genes that the framework did NOT
     select, to check the framework is not merely a degree cut-off

Outputs (in 补强分析_20261114/):
  framework_decirculation.csv        metrics table
  framework_permutation.json         permutation null summary
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "补强分析_20261114"
RNG = np.random.default_rng(20261114)

# --- universe -------------------------------------------------------------
mp = pd.read_csv(ROOT / "target_gene_mapping.csv").dropna(subset=["gene"])
universe = set(mp["gene"].unique())

# --- references -----------------------------------------------------------
pain = set(pd.read_csv(ROOT / "痛觉基因_三方整合总表.csv")["gene"])
lin = pd.read_csv(ROOT / "lincs" / "Enrichr_LINCS_L1000_Chem_Pert_Consensus_Sigs.csv")


def parse_gene_list(s):
    """Enrichr stores overlap_genes as a string that looks like a Python list."""
    return {g.strip().strip("'\"") for g in str(s).strip("[]").split(",") if g.strip()}


l1000 = parse_gene_list(lin.query("term == 'Luteolin Down'")["overlap_genes"].iloc[0])
l1000_up = parse_gene_list(lin.query("term == 'Luteolin Up'")["overlap_genes"].iloc[0])
vote = pd.read_csv(ROOT / "multicohort" / "多队列_跨队列方向投票.csv")
concord = set(vote[vote["concordance"] >= 5 / 6]["gene"])

# genes the framework let in (partly) because they are already known pain genes:
# these cannot serve as an independent reference for the framework either
tgt = json.load(open(ROOT / "target_list_enriched.json", encoding="utf-8"))
pain_seeded = {r["gene"] for r in tgt if "痛觉基因集" in r["sources"]}

refs = {
    "pain_only": pain,                              # external, but see pain_seeded
    "pain_external_only": pain - pain_seeded,        # strictly independent of the rule
    "l1000_only": l1000,                            # experimental, same compound
    "direction_only": concord,                      # same cohorts
    "leave_out_pain": l1000 | concord,
    "leave_out_l1000": pain | concord,
    "leave_out_direction": pain | l1000,
    "composite_any": pain | l1000 | concord,        # original (circular) reference
}
refs = {k: (v & universe) for k, v in refs.items()}

# --- strategies -----------------------------------------------------------
framework = {r["gene"] for r in tgt}
framework_decirc = (framework - pain_seeded) & universe   # selection compared against the null
core = {r["gene"] for r in tgt if r["tier"] in {"T0_上游TF", "T1_核心"}}
strong = {r["gene"] for r in tgt if r["tier"] in {"T0_上游TF", "T1_核心", "T2_强支撑"}}
ppi = pd.read_csv(ROOT / "ppi" / "PPI_节点拓扑.csv").sort_values("degree", ascending=False)
ppi["gene"] = ppi["gene"].astype(str)
degree_top = set(ppi["gene"].head(len(framework)))
degree_not_selected = set(ppi[~ppi["gene"].isin(framework)]["gene"].head(len(framework)))
strict = set(pd.read_csv(ROOT / "交集_DEG_成分靶点.csv")["gene"])

strategies = {
    "all_mapped_baseline": universe,
    f"degree_top{len(framework)}": degree_top,
    f"degree_top{len(framework)}_not_in_framework": degree_not_selected,
    "strict_intersection47": strict,
    "framework_all": framework,
    "framework_minus_pain_seeded": framework_decirc,
    "framework_core6": core,
    "framework_strong15": strong,
}


def metrics(selected, positive):
    selected = set(selected) & universe
    tp = len(selected & positive)
    fp = len(selected - positive)
    fn = len(positive - selected)
    tn = len(universe) - tp - fp - fn
    precision = tp / len(selected) if selected else np.nan
    recall = tp / len(positive) if positive else np.nan
    table = [[tp, fp], [fn, tn]]
    odds, p = fisher_exact(table, alternative="greater")
    return dict(n=len(selected), tp=tp, precision=precision, recall=recall,
                odds_ratio=odds, fisher_p=p)


rows = []
for ref_name, positive in refs.items():
    for strat_name, selected in strategies.items():
        r = metrics(selected, positive)
        r.update(reference=ref_name, strategy=strat_name,
                 n_positive_in_universe=len(positive))
        rows.append(r)
res = pd.DataFrame(rows)[["reference", "n_positive_in_universe", "strategy", "n", "tp",
                          "precision", "recall", "odds_ratio", "fisher_p"]]
for c in ["precision", "recall", "odds_ratio", "fisher_p"]:
    res[c] = pd.to_numeric(res[c], errors="coerce").round(4)
res.to_csv(OUT / "framework_decirculation.csv", index=False, encoding="utf-8-sig")

# --- permutation null for the framework against the pain-only reference ----
pos = refs["pain_external_only"]
pool = np.array(sorted(universe))
obs_tp = len(framework_decirc & pos)
k = len(framework_decirc)
obs_prec = obs_tp / k
null_prec, null_tp = [], []
for _ in range(10000):
    s = set(RNG.choice(pool, size=k, replace=False))
    t = len(s & pos)
    null_tp.append(t)
    null_prec.append(t / k)
null_prec = np.array(null_prec)
emp_p = float((np.sum(null_prec >= obs_prec) + 1) / (len(null_prec) + 1))
# expected overlap purely by chance
exp_overlap = k * len(pos) / len(universe)
perm = dict(n_perm=10000, universe=len(universe), n_selected=k,
            n_positive=len(pos), observed_overlap=obs_tp,
            observed_precision=round(float(obs_prec), 4),
            expected_overlap_by_chance=round(float(exp_overlap), 2),
            enrichment_vs_chance=round(float(obs_tp / exp_overlap), 2) if exp_overlap else None,
            null_precision_mean=round(float(null_prec.mean()), 4),
            null_precision_p95=round(float(np.percentile(null_prec, 95)), 4),
            null_precision_max=round(float(null_prec.max()), 4),
            empirical_p=round(emp_p, 5),
            reference="pain_external_only",
            selection="framework_minus_pain_seeded",
            note="permutation floor is 1/10001; report as p<1e-4",
            seed=20261114)
(OUT / "framework_permutation.json").write_text(
    json.dumps(perm, ensure_ascii=False, indent=1), encoding="utf-8")

print(res.to_string(index=False))
print()
print("permutation null vs pain_only reference:")
print(json.dumps(perm, ensure_ascii=False, indent=1))
