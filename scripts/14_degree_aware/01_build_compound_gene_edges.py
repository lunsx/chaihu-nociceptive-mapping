# -*- coding: utf-8 -*-
"""Build filtered compound-gene edge lists and gene degrees for both formulas.

The degree-aware de-circularization tests need, for every gene, the number of
pharmacokinetic-filtered compounds that map to it. The curated repository does
not ship the large TCMSP HTML-derived master tables, so this script reads them
from the project work root (see scripts/config.py), keeps compounds that pass
OB >= 30 and DL >= 0.18, joins target names to human genes, and writes small
self-contained edge and degree tables.

Inputs (project work root, outside this repository):
  Chaihu: TCMSP_全方成分主表.csv, TCMSP_成分靶点主表.csv, target_gene_mapping.csv
  Jinling: xiaoyao_case/Jinling_成分主表.csv,
           xiaoyao_case/Jinling_成分靶点主表.csv,
           data/processed_small/jinling_gene_mapping.csv (authoritative map)

Outputs (data/processed_small and src/decimap/data):
  <formula>_compound_gene_edges.csv  columns MOL_ID,gene (distinct pairs)
  <formula>_gene_degree.csv          columns gene,degree
"""
import sys
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

PROJECT = PATHS["project_root"]
REPO = PATHS["repo_root"]
OUT = PATHS["processed_small"]
PKG = REPO / "src" / "decimap" / "data"

CONFIGS = {
    "chaihu": {
        "master": PROJECT / "TCMSP_全方成分主表.csv",
        "edges": PROJECT / "TCMSP_成分靶点主表.csv",
        "gene_map": PROJECT / "target_gene_mapping.csv",
        "universe": "chaihu_filtered_gene_universe.csv",
    },
    "jinling": {
        "master": PROJECT / "xiaoyao_case" / "Jinling_成分主表.csv",
        "edges": PROJECT / "xiaoyao_case" / "Jinling_成分靶点主表.csv",
        "gene_map": OUT / "jinling_gene_mapping.csv",
        "universe": "jinling_filtered_gene_universe.csv",
    },
}


def build(cfg):
    master = pd.read_csv(cfg["master"])
    pass_mols = set(master[master["pass_filter"] == True]["MOL_ID"])
    gm = pd.read_csv(cfg["gene_map"]).dropna(subset=["gene"])
    gm = gm[gm["gene"].astype(str).str.strip() != ""]
    n2g = dict(zip(gm["target_name"], gm["gene"]))

    edges = pd.read_csv(cfg["edges"])
    edges = edges[edges["MOL_ID"].isin(pass_mols)].copy()
    edges["gene"] = edges["target_name"].map(n2g)
    cg = edges[["MOL_ID", "gene"]].dropna().drop_duplicates().reset_index(drop=True)
    cg = cg.sort_values(["MOL_ID", "gene"]).reset_index(drop=True)

    degree = (cg.groupby("gene")["MOL_ID"].nunique()
              .sort_values(ascending=False).rename("degree").reset_index()
              .rename(columns={"gene": "gene"}))
    return cg, degree


def main():
    for key, cfg in CONFIGS.items():
        cg, degree = build(cfg)
        edge_name = f"{key}_compound_gene_edges.csv"
        degree_name = f"{key}_gene_degree.csv"
        for dest in (OUT, PKG):
            cg.to_csv(dest / edge_name, index=False, encoding="utf-8")
            degree.to_csv(dest / degree_name, index=False, encoding="utf-8")

        committed = set(pd.read_csv(OUT / cfg["universe"])["gene"])
        rebuilt = set(degree["gene"])
        assert rebuilt == committed, (key, rebuilt ^ committed)
        print(f"{key}: edges={len(cg)} genes={len(degree)} "
              f"universe matches committed={len(committed)}")


if __name__ == "__main__":
    main()
