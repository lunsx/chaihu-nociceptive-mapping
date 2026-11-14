# -*- coding: utf-8 -*-
"""Summarise the multi-conformer x multi-scoring x multi-seed docking ensemble.

Design: 4 PTGS2 structures (5KIR, 6COX, 1CX2, 3PGH) x 6 ligands
(5 formula constituents + celecoxib control) x 2 scoring functions (vina,
vinardo) x 3 seeds (42, 2026, 7), exhaustiveness 8, box 22 A.

What is reported:
  1. mean +/- SD affinity per receptor/ligand/scoring, and the seed-to-seed spread
     (how much a single seed could mislead);
  2. ligand ranking per receptor/scoring, and how often luteolin stays top-1;
  3. rank concordance across structures (Spearman, Kendall) for vina and vinardo;
  4. vina vs vinardo agreement (paired, per receptor/ligand);
  5. where the positive control (celecoxib) lands - a sanity anchor, not proof.

Output: docking_ensemble_summary.csv, docking_ensemble_ranks.csv,
        docking_ensemble_stats.json
"""
import json
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau, wilcoxon

SRC = Path(r"C:\dockwork\ens2\results.csv")
OUT = Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\补强分析_20261114")

LIG_ORDER = ["luteolin", "quercetin", "kaempferol", "cyperotundone",
             "alpha_cyperone", "celecoxib"]


def main():
    df = pd.read_csv(SRC)
    df["ligand"] = df["ligand"].str.replace("-", "_", regex=False)
    df = df[df["ligand"].isin(LIG_ORDER)]
    print("rows:", len(df), "receptors:", sorted(df.receptor.unique()),
          "scorings:", sorted(df.scoring.unique()))

    # 1. per receptor/ligand/scoring summary
    agg = (df.groupby(["receptor", "ligand", "scoring"])["affinity"]
           .agg(["mean", "std", "min", "max", "count"]).reset_index())
    agg["mean"] = agg["mean"].round(3)
    agg["std"] = agg["std"].round(3)
    agg["seed_spread"] = (agg["max"] - agg["min"]).round(3)
    agg = agg.rename(columns={"std": "sd_across_seeds"})
    agg.to_csv(OUT / "docking_ensemble_summary.csv", index=False, encoding="utf-8-sig")

    # 2. rankings
    ranks = []
    for (rec, sc), sub in agg.groupby(["receptor", "scoring"]):
        s = sub.set_index("ligand")["mean"].sort_values()
        for rank, (lig, val) in enumerate(s.items(), start=1):
            ranks.append(dict(receptor=rec, scoring=sc, rank=rank, ligand=lig,
                              mean_affinity=round(float(val), 3)))
    rank_df = pd.DataFrame(ranks)
    rank_df.to_csv(OUT / "docking_ensemble_ranks.csv", index=False, encoding="utf-8-sig")

    top1 = rank_df[rank_df["rank"] == 1].groupby("ligand").size().to_dict()
    luteolin_rank = (rank_df[rank_df["ligand"] == "luteolin"]
                     .groupby("scoring")["rank"].agg(["min", "median", "max"])
                     .to_dict("index"))
    celecoxib_rank = (rank_df[rank_df["ligand"] == "celecoxib"]
                      .groupby("scoring")["rank"].agg(["min", "median", "max"])
                      .to_dict("index"))

    # 3. cross-structure concordance
    conc = {}
    for sc in ["vina", "vinardo"]:
        wide = (agg[agg.scoring == sc].pivot(index="receptor", columns="ligand",
                                             values="mean").reindex(columns=LIG_ORDER))
        pairs = []
        for a, b in itertools.combinations(wide.index, 2):
            rho, p = spearmanr(wide.loc[a], wide.loc[b])
            tau, _ = kendalltau(wide.loc[a], wide.loc[b])
            pairs.append(dict(scoring=sc, r1=a, r2=b, spearman=round(float(rho), 3),
                              spearman_p=round(float(p), 4), kendall=round(float(tau), 3)))
        conc[sc] = pairs
    pd.DataFrame([p for v in conc.values() for p in v]).to_csv(
        OUT / "docking_structure_concordance.csv", index=False, encoding="utf-8-sig")

    # 4. vina vs vinardo
    pair_vv = (agg.pivot_table(index=["receptor", "ligand"], columns="scoring",
                               values="mean"))
    rho_vv, p_vv = spearmanr(pair_vv["vina"], pair_vv["vinardo"])
    diff = (pair_vv["vinardo"] - pair_vv["vina"])
    try:
        w_p = float(wilcoxon(pair_vv["vina"], pair_vv["vinardo"]).pvalue)
    except Exception:  # noqa: BLE001
        w_p = None

    # 5. control anchor
    cel = agg[agg.ligand == "celecoxib"].set_index(["receptor", "scoring"])["mean"]
    lut = agg[agg.ligand == "luteolin"].set_index(["receptor", "scoring"])["mean"]
    anchor = {f"{r}|{s}": round(float(cel.loc[(r, s)] - lut.loc[(r, s)]), 3)
              for (r, s) in cel.index}

    stats = dict(
        n_runs=int(len(df)),
        receptors=sorted(df.receptor.unique()),
        scorings=sorted(df.scoring.unique()),
        seeds=sorted(int(x) for x in df.seed.unique()),
        top1_counts=top1,
        luteolin_rank=luteolin_rank,
        celecoxib_rank=celecoxib_rank,
        vina_vs_vinardo_spearman=round(float(rho_vv), 3),
        vina_vs_vinardo_p=round(float(p_vv), 5),
        vina_vs_vinardo_wilcoxon_p=(round(w_p, 5) if w_p is not None else None),
        vinardo_minus_vina_mean=round(float(diff.mean()), 3),
        vinardo_minus_vina_range=[round(float(diff.min()), 3), round(float(diff.max()), 3)],
        celecoxib_minus_luteolin=anchor,
        mean_seed_spread=round(float(agg["seed_spread"].mean()), 3),
        max_seed_spread=round(float(agg["seed_spread"].max()), 3),
        cross_structure=conc,
        caveat=("docking scores are not binding free energies; the ensemble tests "
                "robustness of ranking, not affinity"),
    )
    (OUT / "docking_ensemble_stats.json").write_text(
        json.dumps(stats, ensure_ascii=False, indent=1), encoding="utf-8")

    print(agg.to_string(index=False))
    print()
    print(json.dumps({k: stats[k] for k in
                      ["top1_counts", "luteolin_rank", "celecoxib_rank",
                       "vina_vs_vinardo_spearman", "vina_vs_vinardo_wilcoxon_p",
                       "mean_seed_spread", "max_seed_spread"]},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
