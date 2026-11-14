# -*- coding: utf-8 -*-
import gzip, os, numpy as np, pandas as pd
from scipy import stats
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"

# ---- data sources ----
mp = pd.read_csv(os.path.join(OUT,"target_gene_mapping.csv")).dropna(subset=["gene"])
gmap = dict(zip(mp["target_name"], mp["gene"]))
T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv")); T["gene"] = T["target_name"].map(gmap)
M = pd.read_csv(os.path.join(OUT,"TCMSP_全方成分主表.csv"))
fmol = set(M[M["pass_filter"]]["MOL_ID"])
deg = pd.read_csv(os.path.join(OUT,"GSE89632_DEG_NASH_vs_HC.csv"))

# ---- GSE130970 NAS + TPM ----
with gzip.open(os.path.join(OUT,"GSE130970_series_matrix.txt.gz"),"rt",errors="replace") as f:
    lines=f.readlines()
names=None; data={}
for ln in lines:
    if ln.startswith("!Sample_title"):
        names=[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
    if ln.startswith("!Sample_characteristics_ch1"):
        vals=[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        lab=vals[0].split(":")[0].strip()
        data[lab]=[v.split(":",1)[1].strip() if ":" in v else v for v in vals]
meta=pd.DataFrame({"sample":names})
meta["nas"]=pd.to_numeric(data["nafld activity score"],errors="coerce")
tpm=pd.read_csv(os.path.join(OUT,"GSE130970_TPM.csv.gz")).set_index("entrez_id")
tpm=tpm[[c for c in tpm.columns if c in set(names)]]
logt=np.log2(tpm+0.5)
nas=meta.set_index("sample").loc[tpm.columns,"nas"].values
pmap=pd.read_csv(os.path.join(OUT,"probe_gene_map.csv"),dtype=str)
pmap["entrez_id"]=pd.to_numeric(pmap["entrez_id"],errors="coerce")
s2e=pmap.dropna(subset=["entrez_id"]).drop_duplicates("symbol").set_index("symbol")["entrez_id"]

PAIN=["TRPV1","PTGS2","TNF","IL1B","IL6","IL1A","PTGES","PTGER3","PTGS1","TRPA1",
      "TACR1","CALCRL","NGF","SCN9A","SCN5A","SCN10A","TRPV3","TRPC3","IL1R1"]
FT=T[T.MOL_ID.isin(fmol) & T.gene.notna()]
fgenes=set(FT.gene)
allgenes=set(T.gene.dropna())
rows=[]
for g in PAIN:
    d=deg[deg.gene==g]
    # GSE130970
    r130=p130=np.nan
    eid=s2e.get(g)
    if pd.notna(eid) and int(eid) in logt.index:
        r130,p130=stats.spearmanr(logt.loc[int(eid)].values,nas)
    # formula compounds
    fsub=FT[FT.gene==g]; asub=T[(T.gene==g)]
    fcomp=sorted(fsub.molecule_name.unique())
    rows.append({
        "gene":g,
        "formula_target_filtered": g in fgenes,
        "formula_target_any": g in allgenes,
        "n_filtered_compounds": len(fcomp),
        "filtered_compounds": "; ".join(fcomp[:6]),
        "herbs": "; ".join(sorted(set(",".join(fsub.herbs).split(",")))) if len(fsub) else "",
        "validated": bool(fsub.validated.any()) if len(fsub) else False,
        "GSE89632_logFC": float(d.logFC_NASH_vs_HC.iloc[0]) if len(d) else np.nan,
        "GSE89632_adj_p": float(d.adj_p_BH.iloc[0]) if len(d) else np.nan,
        "GSE89632_sig_DEG": bool(len(d) and d.adj_p_BH.iloc[0]<0.05 and abs(d.logFC_NASH_vs_HC.iloc[0])>=0.5),
        "GSE130970_r": round(float(r130),3) if not np.isnan(r130) else np.nan,
        "GSE130970_p": float(p130) if not np.isnan(p130) else np.nan,
    })
R=pd.DataFrame(rows)
R.to_csv(os.path.join(OUT,"痛觉基因_三方整合总表.csv"),index=False,encoding="utf-8-sig")
pd.set_option("display.width",250); pd.set_option("display.max_columns",30)
print(R[["gene","formula_target_filtered","n_filtered_compounds","GSE89632_logFC","GSE89632_sig_DEG","GSE130970_r","GSE130970_p","validated"]].to_string(index=False))
# TRPV1 compounds (all, incl non-filtered)
print("\nTRPV1 compounds (all 938):")
print(T[(T.gene=="TRPV1")][["molecule_name","MOL_ID","herbs","validated"]].to_string(index=False))
