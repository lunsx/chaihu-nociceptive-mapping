# -*- coding: utf-8 -*-
"""GSE344087 complete re-analysis: all 20 samples (19 MASLD + 1 healthy control).

The earlier run used only 7 complete MASLD samples cut from a partially downloaded
RAW.tar and had no control. Per-sample h5 files are individually available from
GEO (8-14 MB each), so the whole series can be analysed.

Honest statistics: there is ONE control donor. Cell-level tests across 19 vs 1
donor would be pseudoreplicated, so the output is donor-level descriptive:
per-donor percentages, MASLD median (IQR) versus the single control value.
No p-values are attached to the disease contrast.

Outputs (C:\\dockwork\\masld):
  masld_full_localization.csv   cell_type x gene, per condition and per donor
  masld_full_composition.csv    cell-type fractions per donor and per condition
  masld_full_summary.json       QC, clusters, annotation map, caveats
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

W = Path(__file__).resolve().parent
H5_OLD = W / "h5"
H5_NEW = W / "gsm"
sc.settings.n_jobs = 4

TARGETS = ["IL1B", "TNF", "PTGS2", "IL6", "IL1A", "PTGES", "PTGER1", "PTGER3",
           "CALCRL", "TACR1", "TRPV1", "TNFRSF1A", "IL6R", "IL6ST", "ICAM1",
           "SPP1", "FOS"]

MARKERS = {
    "Hepatocyte": ["ALB", "TTR", "APOA1", "APOB"],
    "Endothelial": ["PECAM1", "KDR", "CLEC4G", "CLEC14A"],
    "Cholangiocyte": ["KRT19", "EPCAM", "KRT7"],
    "Macrophage/Kupffer": ["CD68", "LYZ", "C1QA", "MARCO", "C1QB"],
    "T/NK": ["NKG7", "CD3D", "CD3E", "GNLY", "CD8A"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "Plasma": ["MZB1", "JCHAIN", "IGHG1"],
    "Hepatic stellate": ["COL1A1", "ACTA2", "RGS5", "COL3A1"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
}


def sample_meta(path):
    """GSM#####_{condition}_{donor}_LIVER_... -> (condition, donor, gsm)."""
    parts = path.name.split("_")
    gsm = parts[0]
    cond = "HC" if "HC" in parts else "MASLD"
    donor = next((p for p in parts if p.startswith("S") and p[1:].isdigit()), None)
    if donor is None:
        donor = "S001" if "9970276" in gsm else gsm
    return cond, donor, gsm


def main():
    files = sorted(H5_OLD.glob("*.h5")) + sorted(H5_NEW.glob("*.h5"))
    print("samples:", len(files))
    adatas = []
    meta = []
    for f in files:
        cond, donor, gsm = sample_meta(f)
        a = sc.read_10x_h5(str(f))
        a.var_names_make_unique()
        a.obs["sample"] = f"{cond}_{donor}"
        a.obs["donor"] = donor
        a.obs["condition"] = cond
        meta.append(dict(file=f.name, gsm=gsm, condition=cond, donor=donor,
                         n_cells_raw=int(a.n_obs), n_genes=int(a.n_vars)))
        adatas.append(a)
    adata = sc.concat(adatas, join="outer", fill_value=0, label="batch",
                      keys=[a.obs["sample"].iloc[0] for a in adatas])
    print("merged", adata.shape)
    adata.obs["sample"] = adata.obs["batch"].astype(str)
    adata.obs["condition"] = [s.split("_")[0] for s in adata.obs["sample"]]
    adata.obs["donor"] = [s.split("_")[1] for s in adata.obs["sample"]]

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < 20].copy()
    print(f"QC: {before} -> {adata.n_obs} cells, {adata.n_vars} genes")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat")
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack", random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2,
                 random_state=0)
    print("clusters:", adata.obs.leiden.nunique())

    present = {k: [g for g in v if g in adata.raw.var_names] for k, v in MARKERS.items()}
    scoremat = {}
    for ctype, genes in present.items():
        if genes:
            sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{ctype}",
                              random_state=0)
            scoremat[ctype] = adata.obs[f"score_{ctype}"]
    scores = pd.DataFrame(scoremat)
    cluster_map = {cl: scores.loc[adata.obs.leiden == cl].mean().idxmax()
                   for cl in adata.obs.leiden.cat.categories}
    adata.obs["cell_type"] = pd.Categorical(adata.obs.leiden.map(cluster_map))

    raw = adata.raw.to_adata()

    # --- donor-level percentages ---
    rows = []
    for ct in adata.obs.cell_type.cat.categories:
        for cond in ["MASLD", "HC"]:
            sel = (adata.obs.cell_type == ct) & (adata.obs.condition == cond)
            if sel.sum() == 0:
                continue
            for g in TARGETS:
                if g not in raw.var_names:
                    continue
                per_donor = {}
                for donor in sorted(adata.obs.loc[sel, "donor"].unique()):
                    dsel = sel & (adata.obs.donor == donor)
                    if dsel.sum() < 20:
                        continue
                    x = raw[dsel, g].X
                    x = (np.asarray(x.todense()).ravel()
                         if hasattr(x, "todense") else np.asarray(x).ravel())
                    per_donor[donor] = round(float(100 * (x > 0).mean()), 2)
                x_all = raw[sel, g].X
                x_all = (np.asarray(x_all.todense()).ravel()
                         if hasattr(x_all, "todense") else np.asarray(x_all).ravel())
                vals = np.array(list(per_donor.values()), dtype=float)
                rows.append(dict(
                    cell_type=ct, condition=cond, gene=g, n_cells=int(sel.sum()),
                    n_donors=len(per_donor),
                    pct_all_cells=round(float(100 * (x_all > 0).mean()), 2),
                    pct_donor_median=(round(float(np.median(vals)), 2) if len(vals) else None),
                    pct_donor_min=(round(float(vals.min()), 2) if len(vals) else None),
                    pct_donor_max=(round(float(vals.max()), 2) if len(vals) else None),
                    _donor_values=";".join(f"{k}:{v}" for k, v in per_donor.items()),
                ))
    loc = pd.DataFrame(rows).drop(columns=["_donor_values"])
    loc.to_csv(W / "masld_full_localization.csv", index=False, encoding="utf-8-sig")

    # wide contrast table
    wide = loc.pivot_table(index=["cell_type", "gene"],
                           columns="condition",
                           values="pct_donor_median", aggfunc="first").reset_index()
    wide.columns = [c if c in ("cell_type", "gene") else f"pct_{c}" for c in wide.columns]
    wide["delta_MASLD_minus_HC"] = (wide.get("pct_MASLD", np.nan)
                                    - wide.get("pct_HC", np.nan)).round(2)
    wide.to_csv(W / "masld_full_contrast.csv", index=False, encoding="utf-8-sig")

    # --- composition ---
    comp = (adata.obs.groupby(["donor", "condition", "cell_type"], observed=True)
            .size().rename("n_cells").reset_index())
    comp["fraction"] = (comp["n_cells"]
                        / comp.groupby("donor")["n_cells"].transform("sum")).round(4)
    comp.to_csv(W / "masld_full_composition.csv", index=False, encoding="utf-8-sig")

    summary = dict(
        total_cells=int(adata.n_obs), total_genes=int(adata.n_vars),
        n_samples=int(adata.obs["sample"].nunique()),
        conditions=adata.obs.condition.value_counts().to_dict(),
        n_cell_types=int(adata.obs.cell_type.nunique()),
        cluster_map=cluster_map, samples=meta,
        caveat=("single healthy-control donor; disease contrast is descriptive "
                "donor-level, no p-values, no donor-level inference"),
    )
    (W / "masld_full_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in
                      ["total_cells", "total_genes", "n_samples", "conditions"]},
                     ensure_ascii=False, indent=1))

    adata.obs.drop(columns=[c for c in adata.obs.columns if c.startswith("score_")],
                   inplace=True, errors="ignore")
    adata.write_h5ad(W / "masld_full_processed.h5ad")
    print("DONE")


if __name__ == "__main__":
    main()
