# -*- coding: utf-8 -*-
"""Round-6 open items (#130):

1. Gene-set-level permutation: empirical p-value for every reference x
   strategy combination (10,000 random draws of the same size).
2. Direction 3x3 criterion: predicted direction (target list) vs observed
   multi-cohort direction, classified as up/down/uncertain on both axes.
3. Docking ranks split by scoring function and within-panel z-scores.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析")
OUT = ROOT / "补强分析_20261114"
RNG = np.random.default_rng(20261114)

# ---------------------------------------------------------------- references
mp = pd.read_csv(ROOT / "target_gene_mapping.csv").dropna(subset=["gene"])
universe = set(mp["gene"].unique())
pain = set(pd.read_csv(ROOT / "痛觉基因_三方整合总表.csv")["gene"])
lin = pd.read_csv(ROOT / "lincs" / "Enrichr_LINCS_L1000_Chem_Pert_Consensus_Sigs.csv")


def parse_gene_list(s):
    return {g.strip().strip("'\"") for g in str(s).strip("[]").split(",") if g.strip()}


l1000 = parse_gene_list(lin.query("term == 'Luteolin Down'")["overlap_genes"].iloc[0])
l1000_up = parse_gene_list(lin.query("term == 'Luteolin Up'")["overlap_genes"].iloc[0])
vote0 = pd.read_csv(ROOT / "multicohort" / "多队列_跨队列方向投票.csv")
concord = set(vote0[vote0["concordance"] >= 5 / 6]["gene"])
tgt = json.load(open(ROOT / "target_list_enriched.json", encoding="utf-8"))
pain_seeded = {r["gene"] for r in tgt if "痛觉基因集" in r["sources"]}

refs = {
    "pain_only": pain,
    "pain_external_only": pain - pain_seeded,
    "l1000_only": l1000,
    "direction_only": concord,
    "leave_out_pain": l1000 | concord,
    "leave_out_l1000": pain | concord,
    "leave_out_direction": pain | l1000,
    "composite_any": pain | l1000 | concord,
}
refs = {k: (v & universe) for k, v in refs.items()}

framework = {r["gene"] for r in tgt}
framework_decirc = (framework - pain_seeded) & universe
core = {r["gene"] for r in tgt if r["tier"] in {"T0_上游TF", "T1_核心"}}
strong = {r["gene"] for r in tgt if r["tier"] in {"T0_上游TF", "T1_核心", "T2_强支撑"}}
ppi = pd.read_csv(ROOT / "ppi" / "PPI_节点拓扑.csv").sort_values("degree", ascending=False)
ppi["gene"] = ppi["gene"].astype(str)
degree_top = set(ppi["gene"].head(len(framework)))
degree_not_selected = set(ppi[~ppi["gene"].isin(framework)]["gene"].head(len(framework)))
strict = set(pd.read_csv(ROOT / "交集_DEG_成分靶点.csv")["gene"])
strategies = {
    f"degree_top{len(framework)}": degree_top,
    f"degree_top{len(framework)}_not_in_framework": degree_not_selected,
    "strict_intersection47": strict,
    "framework_all": framework,
    "framework_minus_pain_seeded": framework_decirc,
    "framework_core6": core,
    "framework_strong15": strong,
}
strategies = {k: (v & universe) for k, v in strategies.items()}

# ------------------------------------------------- 1. gene-set permutations
pool = np.array(sorted(universe))
perm_rows = []
null_cache = {}
for ref_name, positive in refs.items():
    pos_idx = np.array([g in positive for g in pool])
    for strat_name, selected in strategies.items():
        k = len(selected)
        obs_tp = len(selected & positive)
        if (ref_name, k) not in null_cache:
            null_tp = np.empty(10000, dtype=np.int16)
            for i in range(10000):
                idx = RNG.choice(len(pool), size=k, replace=False)
                null_tp[i] = int(pos_idx[idx].sum())
            null_cache[(ref_name, k)] = null_tp
        null_tp = null_cache[(ref_name, k)]
        emp_p = (np.sum(null_tp >= obs_tp) + 1) / (len(null_tp) + 1)
        table = [[obs_tp, k - obs_tp],
                 [len(positive) - obs_tp,
                  len(universe) - k - len(positive) + obs_tp]]
        fisher_p = fisher_exact(table, alternative="greater")[1]
        perm_rows.append({
            "reference": ref_name, "strategy": strat_name, "set_size": k,
            "reference_size": len(positive), "observed_overlap": obs_tp,
            "null_mean_overlap": round(float(null_tp.mean()), 3),
            "null_p95_overlap": float(np.percentile(null_tp, 95)),
            "empirical_p": round(float(emp_p), 5), "fisher_p": round(fisher_p, 5),
        })
perm_df = pd.DataFrame(perm_rows)
perm_df.to_csv(OUT / "framework_gene_set_permutation.csv", index=False, encoding="utf-8-sig")
print("=== gene-set permutation (empirical p<0.05) ===")
print(perm_df[perm_df["empirical_p"] < 0.05].to_string(index=False))

# ------------------------------------------------- 2. direction 3x3 matrix
def three_class(d):
    if d in ("up", "down"):
        return d
    return "uncertain"


obs = vote0.set_index("gene")["direction"].to_dict()
obs_conc = vote0.set_index("gene")["concordance"].to_dict()
matrix = pd.DataFrame(0, index=["pred_up", "pred_down", "pred_uncertain"],
                      columns=["obs_up", "obs_down", "obs_uncertain"])
detail = []
for r in tgt:
    g = r["gene"]
    if g not in obs:
        continue
    p = three_class(r.get("direction"))
    o = three_class(obs[g])
    if o != "uncertain" and obs_conc.get(g, 0) < 5 / 6:
        o = "uncertain"
    matrix.loc[f"pred_{p}", f"obs_{o}"] += 1
    detail.append({"gene": g, "predicted": p, "observed": o,
                   "observed_concordance": round(float(obs_conc.get(g, 0)), 3)})
matrix.to_csv(OUT / "direction_3x3_matrix.csv", encoding="utf-8-sig")
pd.DataFrame(detail).to_csv(OUT / "direction_3x3_detail.csv", index=False, encoding="utf-8-sig")
diag = sum(matrix.iloc[i, i] for i in range(3))
total = int(matrix.values.sum())
print("\n=== direction 3x3 ===")
print(matrix)
print(f"agreement: {diag}/{total} = {diag / total:.3f}")

# --------------------------------------- 3. ranks by scoring + panel z-score
res = pd.read_csv(r"C:\dockwork\ens2\results.csv")
avg = (res.groupby(["receptor", "scoring", "ligand"])["affinity"]
       .mean().reset_index())
rank_rows = []
z_rows = []
for (rec, sc), gdf in avg.groupby(["receptor", "scoring"]):
    gdf = gdf.sort_values("affinity")  # more negative first
    for rank, (_, rr) in enumerate(gdf.iterrows(), 1):
        rank_rows.append({"receptor": rec, "scoring": sc, "ligand": rr["ligand"],
                          "mean_affinity": round(rr["affinity"], 3), "rank": rank})
    mu, sd = gdf["affinity"].mean(), gdf["affinity"].std(ddof=0)
    for _, rr in gdf.iterrows():
        # positive z = stronger (more negative) than the panel mean
        z = -(rr["affinity"] - mu) / sd if sd > 0 else 0.0
        z_rows.append({"receptor": rec, "scoring": sc, "ligand": rr["ligand"],
                       "mean_affinity": round(rr["affinity"], 3),
                       "panel_z": round(float(z), 3)})
rank_df = pd.DataFrame(rank_rows)
z_df = pd.DataFrame(z_rows)
rank_df.to_csv(OUT / "docking_ranks_by_scoring.csv", index=False, encoding="utf-8-sig")
z_df.to_csv(OUT / "docking_panel_zscore.csv", index=False, encoding="utf-8-sig")
print("\n=== panel z-scores ===")
print(z_df.to_string(index=False))
