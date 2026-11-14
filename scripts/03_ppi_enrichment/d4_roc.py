# -*- coding: utf-8 -*-
import gzip, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
plt.rcParams["font.sans-serif"] = ["SimHei"]; plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
RC = os.path.join(OUT, "roc"); os.makedirs(RC, exist_ok=True)

with gzip.open(os.path.join(OUT,"GSE89632_series_matrix.txt.gz"),"rt",errors="replace") as f:
    lines=f.readlines()
def meta(prefix):
    for ln in lines:
        if ln.startswith(prefix): return [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
gsm=meta("!Sample_geo_accession")
diag=None
for ln in lines:
    if ln.startswith("!Sample_characteristics_ch1"):
        vals=[x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        if all(v.startswith("diagnosis:") for v in vals):
            diag=[v.split(":",1)[1].strip() for v in vals]; break
tb=[]; inside=False
for ln in lines:
    if ln.startswith("!series_matrix_table_begin"): inside=True; continue
    if ln.startswith("!series_matrix_table_end"): break
    if inside: tb.append([x.strip('"') for x in ln.rstrip("\n").split("\t")])
expr=pd.DataFrame(tb[1:],columns=tb[0]).set_index("ID_REF").apply(pd.to_numeric)
expr.columns=gsm
pmap=pd.read_csv(os.path.join(OUT,"probe_gene_map.csv"),dtype=str).set_index("probe_id")
expr=expr.join(pmap[["symbol"]],how="left")
m=expr.symbol.notna() & ~expr.symbol.str.startswith("ERCC") & (expr.symbol.str.strip()!="")
expr=expr[m].copy()
expr["mean_all"]=expr[gsm].mean(axis=1)
expr=expr.sort_values("mean_all",ascending=False).drop_duplicates("symbol").set_index("symbol")

nash=[g for g,d in zip(gsm,diag) if d=="NASH"]; hc=[g for g,d in zip(gsm,diag) if d=="HC"]
keep=nash+hc
y=np.array([1 if g in nash else 0 for g in keep])

candidates=["IL1B","TNF","PTGS2","IL6","IL1A","PTGES","PTGER3","AKT1","TP53","STAT3","SCN5A"]
candidates=[g for g in candidates if g in expr.index]

def boot_auc(score,y,n_boot=1000,seed=11):
    rng=np.random.default_rng(seed); point=auc(*roc_curve(y,score)[:2]); vals=[]
    for _ in range(n_boot):
        idx=np.concatenate([rng.choice(np.where(y==c)[0],len(np.where(y==c)[0]),replace=True) for c in [0,1]])
        vals.append(auc(*roc_curve(y[idx],score[idx])[:2]))
    return point,np.percentile(vals,[2.5,97.5])

rows=[]; curves={}
for g in candidates:
    raw=expr.loc[g,keep].values.astype(float)
    logfc=raw[y==1].mean()-raw[y==0].mean()
    direction="NASH上调" if logfc>0 else "NASH下调"
    a,ci=boot_auc(raw,y)
    oriented=max(a,1-a)                      # separation ability regardless of direction
    # orient score so curve sits above diagonal
    s_orient=raw if a>=0.5 else -raw
    fpr,tpr,_=roc_curve(y,s_orient)
    curves[g]=(fpr,tpr,oriented,direction)
    rows.append({"gene":g,"AUC_原始(预测NASH)":round(a,3),"区分度AUC":round(oriented,3),
                 "方向":direction,"CI_low":round(min(ci[0],1-ci[1]),3),
                 "CI_high":round(max(ci[1],1-ci[0]),3)})
R=pd.DataFrame(rows).sort_values("区分度AUC",ascending=False)
R.to_csv(os.path.join(RC,"ROC_AUC汇总.csv"),index=False,encoding="utf-8-sig")
print(R.to_string(index=False))

core=["IL1B","TNF","PTGS2","IL6","IL1A","PTGES","PTGER3"]
fig,ax=plt.subplots(figsize=(9,8))
for i,g in enumerate([c for c in core if c in curves]):
    fpr,tpr,a,direction=curves[g]
    ax.plot(fpr,tpr,lw=2,color=plt.cm.tab10.colors[i],
            label=f"{g} (AUC={a:.2f}, {direction})")
ax.plot([0,1],[0,1],"--",color="gray",lw=1)
ax.set_xlabel("False Positive Rate"); ax.set_ylabel("True Positive Rate")
ax.set_title("核心靶点区分 NASH 与 HC 的 ROC（已按方向正向化；GSE89632）",fontsize=12)
ax.legend(loc="lower right",fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(RC,"fig_ROC曲线.png"),dpi=200); plt.close()

fig,ax=plt.subplots(figsize=(10,6))
Rp=R.sort_values("区分度AUC")
ax.barh(Rp.gene,Rp["区分度AUC"],
        xerr=[Rp["区分度AUC"]-Rp.CI_low,Rp.CI_high-Rp["区分度AUC"]],
        color=["#C0392B" if g in core else "#4A90D9" for g in Rp.gene],capsize=4)
ax.axvline(0.5,color="gray",ls="--"); ax.set_xlim(0.4,1.0)
ax.set_xlabel("区分度 AUC（误差棒为 95% bootstrap CI）")
ax.set_title("各靶点区分 NASH 与 HC 的能力（不分方向）",fontsize=13)
plt.tight_layout(); plt.savefig(os.path.join(RC,"fig_AUC条形.png"),dpi=200); plt.close()
print("\nfiles:",os.listdir(RC))
