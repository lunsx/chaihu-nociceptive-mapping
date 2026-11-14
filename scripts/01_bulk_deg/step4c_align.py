# -*- coding: utf-8 -*-
import gzip, os, numpy as np, pandas as pd
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
SM = os.path.join(OUT, "GSE89632_series_matrix.txt.gz")
with gzip.open(SM,"rt",encoding="utf-8",errors="replace") as f:
    lines=f.readlines()
rows=[]
for ln in lines:
    if ln.startswith("!Sample_characteristics_ch1"):
        vals=[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        rows.append((vals[0].split(":")[0].strip(), vals))
print("field counts per characteristic row:")
for lab,vals in rows:
    print(f"  {lab!r:35s} n={len(vals)}")
# diagnosis and bmi
diag=[x.split(":",1)[1].strip() for x in dict(rows)["diagnosis"]]
bmi_raw=dict(rows)["body mass index (kg/m2)"]
bmi =[x.split(":",1)[1].strip() if ":" in x else "" for x in bmi_raw]
bmi=pd.to_numeric(pd.Series(bmi),errors="coerce")
d=pd.DataFrame({"diag":diag,"bmi":bmi})
print("\nHC BMI raw values:\n", d[d.diag=="HC"]["bmi"].tolist())
print("\nNASH BMI raw values:\n", d[d.diag=="NASH"]["bmi"].tolist())

# all probes for PTGS2
pmap=pd.read_csv(os.path.join(OUT,"probe_gene_map.csv"),dtype=str)
tb=[]; inside=False
for ln in lines:
    if ln.startswith("!series_matrix_table_begin"): inside=True; continue
    if ln.startswith("!series_matrix_table_end"): break
    if inside: tb.append([x.strip('"') for x in ln.rstrip("\n").split("\t")])
expr=pd.DataFrame(tb[1:],columns=tb[0]).set_index("ID_REF").apply(pd.to_numeric)
expr.columns=diag
print("\nall probes mapping to PTGS2:")
for pr in pmap[pmap.symbol=="PTGS2"]["probe_id"]:
    if pr in expr.index:
        hc=expr.loc[pr,[x=="HC" for x in diag]].mean()
        na=expr.loc[pr,[x=="NASH" for x in diag]].mean()
        print(f"  {pr}: HC={hc:.2f} NASH={na:.2f} logFC={na-hc:+.2f}")
