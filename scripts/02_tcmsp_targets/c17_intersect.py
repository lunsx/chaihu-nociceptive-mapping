# -*- coding: utf-8 -*-
import os, pandas as pd
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
mp = pd.read_csv(os.path.join(OUT,"target_gene_mapping.csv")).dropna(subset=["gene"])
gmap = dict(zip(mp["target_name"], mp["gene"]))

T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv"))
M = pd.read_csv(os.path.join(OUT,"TCMSP_全方成分主表.csv"))
fmol = set(M[M["pass_filter"]]["MOL_ID"])

# attach human gene to every pair; keep filtered-molecule universe
T["gene"] = T["target_name"].map(gmap)
FT = T[T["MOL_ID"].isin(fmol) & T["gene"].notna()].copy()
print("filtered-molecule pairs with human gene:", len(FT),
      "| unique human genes:", FT["gene"].nunique())

# GEO DEG (significant only)
deg = pd.read_csv(os.path.join(OUT,"GSE89632_DEG_NASH_vs_HC.csv"))
degcol = "gene"
LOGFC, ADJP = "logFC_NASH_vs_HC", "adj_p_BH"
sig = deg[(deg[ADJP] < 0.05) & (deg[LOGFC].abs() >= 0.5)].copy()
DEG = set(sig[degcol])
print("significant DEG:", len(DEG))

# pain gene universe (9 core + pathway)
PAIN = {"TRPV1","PTGS2","TNF","IL1B","TRPA1","TACR1","CALCRL","NGF","SCN9A","IL6",
        "PTGS1","SCN10A","TRPV3","TRPC3","IL1A","IL1R1","SCN5A","PTGES","PTGER3"}
formula_genes = set(FT["gene"])
print("\n=== pain genes hit by filtered formula compounds ===")
for g in sorted(PAIN & formula_genes):
    print("  ", g)
print("\n=== three-way: DEG ∩ pain ∩ formula ===")
three = sorted(DEG & PAIN & formula_genes)
print(three)
print("\n=== DEG ∩ formula (all disease-relevant human targets) count:", len(DEG & formula_genes))
# detail table for the three-way
rows = []
for g in sorted(DEG & formula_genes):
    sub = FT[FT["gene"] == g]
    drow = sig[sig[degcol] == g]
    rows.append({
        "gene": g,
        "in_pain_set": g in PAIN,
        "logFC": float(drow[LOGFC].iloc[0]) if len(drow) else None,
        "adj_p": float(drow[ADJP].iloc[0]) if len(drow) else None,
        "n_compounds": sub["MOL_ID"].nunique(),
        "compounds": "; ".join(sorted(sub["molecule_name"].unique())[:8]),
        "herbs": "; ".join(sorted(set(",".join(sub["herbs"]).split(",")))),
        "any_validated": bool(sub["validated"].any()),
    })
R = pd.DataFrame(rows).sort_values(["in_pain_set","adj_p"], ascending=[False, True])
R.to_csv(os.path.join(OUT,"交集_DEG_成分靶点.csv"), index=False, encoding="utf-8-sig")
print("\nsaved", len(R), "DEG∩formula genes; pain ones on top.")
print(R[R["in_pain_set"]].to_string(index=False))
