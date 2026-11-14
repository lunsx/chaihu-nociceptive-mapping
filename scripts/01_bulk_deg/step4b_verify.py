# -*- coding: utf-8 -*-
import gzip, os, numpy as np, pandas as pd
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
SM = os.path.join(OUT, "GSE89632_series_matrix.txt.gz")
with gzip.open(SM, "rt", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()

def char_rows():
    rows={}
    for ln in lines:
        if ln.startswith("!Sample_characteristics_ch1"):
            vals=[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            label=vals[0].split(":")[0].strip()
            rows[label]=[v.split(":",1)[1].strip() if ":" in v else v for v in vals]
    return rows
C=char_rows()
diag=C["diagnosis"]
df=pd.DataFrame({"diag":diag,"bmi":pd.to_numeric(C["body mass index (kg/m2)"],errors="coerce"),
                 "age":pd.to_numeric(C["age (y)"],errors="coerce")})
print(df.groupby("diag").agg(n=("diag","size"),bmi_mean=("bmi","mean"),bmi_sd=("bmi","std"),
                             age_mean=("age","mean")))
print("\ndiabetes by group:")
print(pd.crosstab(pd.Series(diag), pd.Series(C["gender"])))
print(pd.crosstab(pd.Series(diag), pd.Series(C["diabetes"])))

# expression per group for key genes
tb=[]; inside=False
for ln in lines:
    if ln.startswith("!series_matrix_table_begin"): inside=True; continue
    if ln.startswith("!series_matrix_table_end"): break
    if inside: tb.append([x.strip('"') for x in ln.rstrip("\n").split("\t")])
expr=pd.DataFrame(tb[1:],columns=tb[0]).set_index("ID_REF").apply(pd.to_numeric)
expr.columns=diag
for g in ["ILMN_1763782"]:
    pass
# map probe via symbol map
pmap=pd.read_csv(os.path.join(OUT,"probe_gene_map.csv"),dtype=str)
for sym in ["PTGS2","IL6","TRPV1","SOCS3"]:
    pr=pmap[pmap["symbol"]==sym]["probe_id"].iloc[0]
    row=expr.loc[pr]
    print(f"\n{sym} ({pr}) per-group summary:")
    for grp in ["HC","SS","NASH"]:
        v=row[[d==grp for d in diag]].values
        print(f"  {grp:4s} mean={v.mean():.2f} min={v.min():.2f} max={v.max():.2f}")
