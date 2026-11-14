# -*- coding: utf-8 -*-
import gzip, os, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
SM = os.path.join(OUT, "GSE89632_series_matrix.txt.gz")
MAP = os.path.join(OUT, "probe_gene_map.csv")

# ---------- 1. parse series matrix + groups ----------
with gzip.open(SM, "rt", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()
def meta(prefix):
    for ln in lines:
        if ln.startswith(prefix):
            return [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
gsm = meta("!Sample_geo_accession")
diag_rows = [[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]][0]
             for ln in lines if ln.startswith("!Sample_characteristics_ch1")]
# find diagnosis characteristic row
diag = None
for ln in lines:
    if ln.startswith("!Sample_characteristics_ch1"):
        vals = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        if all(v.startswith("diagnosis:") for v in vals):
            diag = [v.split(":",1)[1].strip() for v in vals]; break
assert diag and len(diag) == len(gsm)

# expression table
tb = []
inside = False
for ln in lines:
    if ln.startswith("!series_matrix_table_begin"): inside=True; continue
    if ln.startswith("!series_matrix_table_end"): break
    if inside: tb.append([x.strip('"') for x in ln.rstrip("\n").split("\t")])
header = tb[0]
expr = pd.DataFrame(tb[1:], columns=header).set_index("ID_REF")
expr = expr.apply(pd.to_numeric)
expr.columns = gsm
print("expression matrix:", expr.shape)
print("groups:", pd.Series(diag, index=gsm).value_counts().to_dict())

# ---------- 2. annotate probes ----------
pmap = pd.read_csv(MAP, dtype=str).set_index("probe_id")
expr = expr.join(pmap[["symbol","entrez_id"]], how="left")
# drop controls / missing symbols
mask = expr["symbol"].notna() & ~expr["symbol"].str.startswith("ERCC") & (expr["symbol"].str.strip()!="")
expr = expr[mask].copy()
print("probes after removing controls:", len(expr))

# ---------- 3. collapse probes to gene (keep highest mean expression) ----------
sample_cols = gsm
expr["mean_all"] = expr[sample_cols].mean(axis=1)
expr = expr.sort_values("mean_all", ascending=False)
expr = expr.drop_duplicates(subset="symbol", keep="first")
expr = expr.set_index("symbol")
print("genes after collapse:", len(expr))

# ---------- 4. NASH vs HC differential expression ----------
nash = [g for g,d in zip(gsm,diag) if d=="NASH"]
hc   = [g for g,d in zip(gsm,diag) if d=="HC"]
Xn = expr[nash].values; Xh = expr[hc].values
logFC = Xn.mean(axis=1) - Xh.mean(axis=1)          # data already log2
tstat, pval = stats.ttest_ind(Xn, Xh, axis=1, equal_var=False, nan_policy="omit")
padj = multipletests(pval, method="fdr_bh")[1]

res = pd.DataFrame({
    "gene": expr.index,
    "logFC_NASH_vs_HC": logFC,
    "mean_NASH": Xn.mean(axis=1),
    "mean_HC": Xh.mean(axis=1),
    "t_stat": tstat,
    "p_value": pval,
    "adj_p_BH": padj,
}).sort_values("p_value")
res.to_csv(os.path.join(OUT,"GSE89632_DEG_NASH_vs_HC.csv"), index=False, encoding="utf-8-sig")

sig = res[(res["adj_p_BH"]<0.05) & (res["logFC_NASH_vs_HC"].abs()>=0.5)]
up = sig[sig["logFC_NASH_vs_HC"]>0]; down = sig[sig["logFC_NASH_vs_HC"]<0]
print(f"\nDEG (adj p<0.05 & |logFC|>=0.5): {len(sig)}  up={len(up)} down={len(down)}")
print("top 15 by p:\n", res.head(15)[["gene","logFC_NASH_vs_HC","p_value","adj_p_BH"]].to_string(index=False))
