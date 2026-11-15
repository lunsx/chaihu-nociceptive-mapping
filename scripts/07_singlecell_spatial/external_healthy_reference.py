# -*- coding: utf-8 -*-
"""External healthy human liver reference: mapping and suspension sensitivity.

This script documents the independent healthy liver atlas (Aizarani integrated,
CELLxGENE 135bdd30) used as an external normal reference. It is fully
descriptive: the healthy atlas and the MASLD collection differ in platform,
processing, and annotation and share no matched controls, so no cross-dataset
p-values are computed.

It reads the source h5ad (path resolved from config / environment, no
hard-coded machine paths) and writes three redistributable tables to
``data/processed_small/``:

  1. healthy_celltype_mapping.csv
       12 author-defined cell types, their mapped MASLD label, cell counts,
       cell/nucleus suspension counts, and comparability notes.
  2. healthy_gene_pct_by_author_type.csv
       percent of cells with a nonzero raw count for every author type and
       target gene, overall and within each suspension type.
  3. external_healthy_suspension_sensitivity.csv
       macrophage/Kupffer lineage percentages split by suspension type, testing
       whether the low healthy expression is an artifact of nuclear preparations.

Run:
  python scripts/07_singlecell_spatial/external_healthy_reference.py
"""
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

OUT = PATHS["processed_small"]
H5AD = Path(os.environ.get(
    "CHAIHU_HEALTHY_H5AD",
    PATHS["healthy_liver"] / "healthy_liver_integrated.h5ad",
))

GENES = [
    "IL1B", "TNF", "PTGS2", "IL6", "IL1A", "PTGES", "PTGER1", "PTGER3",
    "CALCRL", "TACR1", "TRPV1", "TNFRSF1A", "IL6R", "IL6ST", "ICAM1",
    "SPP1", "FOS",
]

# 12 author-defined cell types -> comparable MASLD label.
AUTHOR_TO_LABEL = {
    "Bcells": "B",
    "CentralHep": "Hepatocyte",
    "Cholangiocyte": "Cholangiocyte",
    "Erythroid": "Erythroid",
    "InfMac": "Macrophage/Kupffer",
    "InterHep": "Hepatocyte",
    "NKTcell": "T/NK",
    "NonInfMac": "Macrophage/Kupffer",
    "PortalEndo": "Endothelial",
    "PortalHep": "Hepatocyte",
    "Stellate": "Hepatic stellate",
    "cvLSECs": "Endothelial",
}

# MASLD labels that have no counterpart in the healthy atlas.
ABSENT_IN_HEALTHY = ["Mast", "Plasma"]


def _categorical(group, name):
    g = group[name]
    categories = g["categories"].asstr()[:]
    codes = g["codes"][:]
    return categories, codes


def _pct_nonzero(matrix, gene_index):
    """Percentage of rows with a nonzero value in one column."""
    if matrix.shape[0] == 0:
        return np.nan
    col = np.asarray(matrix[:, gene_index].sum(axis=1)).ravel()
    return float((col > 0).mean() * 100)


def main():
    with h5py.File(H5AD, "r") as f:
        sus_cats, sus_codes = _categorical(f["obs"], "suspension_type")
        au_cats, au_codes = _categorical(f["obs"], "author_cell_type")

        feat_cats = f["raw"]["var"]["feature_name"]["categories"].asstr()[:]
        feat_codes = f["raw"]["var"]["feature_name"]["codes"][:]
        features = feat_cats[feat_codes]

        x = f["raw"]["X"]
        matrix = sparse.csr_matrix(
            (x["data"][:], x["indices"][:], x["indptr"][:]),
            shape=(len(sus_codes), len(features)),
        )

    obs = pd.DataFrame({
        "suspension": sus_cats[sus_codes],
        "author_type": au_cats[au_codes],
    })
    gene_index = {
        g: int(np.where(features == g)[0][0])
        for g in GENES if g in set(features)
    }

    # --- 1. cell-type mapping -------------------------------------------------
    mapping_rows = []
    for author_type in au_cats:
        type_mask = (obs["author_type"] == author_type).to_numpy()
        n_total = int(type_mask.sum())
        n_cell = int(((obs["suspension"] == "cell") & type_mask).sum())
        n_nucleus = int(((obs["suspension"] == "nucleus") & type_mask).sum())
        mapped = AUTHOR_TO_LABEL.get(author_type, "")
        note = "healthy-only; no MASLD counterpart" if mapped == "Erythroid" else ""
        mapping_rows.append({
            "author_cell_type": author_type,
            "mapped_MASLD_label": mapped,
            "n_cells": n_total,
            "n_cell_suspension": n_cell,
            "n_nucleus_suspension": n_nucleus,
            "note": note,
        })
    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(OUT / "healthy_celltype_mapping.csv", index=False)

    # --- 2. gene percentages by author type -----------------------------------
    pct_rows = []
    for author_type in au_cats:
        type_mask = (obs["author_type"] == author_type).to_numpy()
        cell_mask = type_mask & (obs["suspension"] == "cell").to_numpy()
        nucleus_mask = type_mask & (obs["suspension"] == "nucleus").to_numpy()
        base = {
            "author_cell_type": author_type,
            "mapped_MASLD_label": AUTHOR_TO_LABEL.get(author_type, ""),
            "n_cells": int(type_mask.sum()),
            "n_cell_suspension": int(cell_mask.sum()),
            "n_nucleus_suspension": int(nucleus_mask.sum()),
        }
        for gene in GENES:
            if gene not in gene_index:
                base[gene + "_pct_all"] = np.nan
                base[gene + "_pct_cell"] = np.nan
                base[gene + "_pct_nucleus"] = np.nan
                continue
            gi = gene_index[gene]
            base[gene + "_pct_all"] = round(_pct_nonzero(matrix[type_mask], gi), 3)
            base[gene + "_pct_cell"] = round(_pct_nonzero(matrix[cell_mask], gi), 3)
            base[gene + "_pct_nucleus"] = round(_pct_nonzero(matrix[nucleus_mask], gi), 3)
        pct_rows.append(base)
    pct = pd.DataFrame(pct_rows)
    pct.to_csv(OUT / "healthy_gene_pct_by_author_type.csv", index=False)

    # --- 3. suspension sensitivity for the macrophage lineage -----------------
    lineage_types = ["InfMac", "NonInfMac"]
    sens_rows = []
    for lineage in lineage_types + ["InfMac+NonInfMac"]:
        if lineage == "InfMac+NonInfMac":
            mask = obs["author_type"].isin(lineage_types).to_numpy()
        else:
            mask = (obs["author_type"] == lineage).to_numpy()
        for suspension in ["cell", "nucleus", "ALL"]:
            sub_mask = mask if suspension == "ALL" else (
                mask & (obs["suspension"] == suspension).to_numpy()
            )
            row = {
                "lineage": lineage,
                "suspension": suspension,
                "n_cells": int(sub_mask.sum()),
            }
            for gene in GENES:
                if gene not in gene_index:
                    row[gene] = np.nan
                    continue
                row[gene] = round(_pct_nonzero(matrix[sub_mask], gene_index[gene]), 3)
            sens_rows.append(row)
    sens = pd.DataFrame(sens_rows)
    sens.to_csv(OUT / "external_healthy_suspension_sensitivity.csv", index=False)

    print("wrote healthy_celltype_mapping.csv rows:", len(mapping))
    print("wrote healthy_gene_pct_by_author_type.csv rows:", len(pct))
    print("wrote external_healthy_suspension_sensitivity.csv rows:", len(sens))
    print("MASLD labels absent in healthy atlas:", ABSENT_IN_HEALTHY)


if __name__ == "__main__":
    main()
