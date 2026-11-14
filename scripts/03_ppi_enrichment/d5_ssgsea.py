# -*- coding: utf-8 -*-
import gzip, os, json, socket, numpy as np, pandas as pd, requests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
plt.rcParams["font.sans-serif"]=["SimHei"]; plt.rcParams["axes.unicode_minus"]=False

OUT=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
IM=os.path.join(OUT,"immune"); os.makedirs(IM,exist_ok=True)
PX={"http":None,"https":None}

# ---- 1. get CellMarker GMT (cached) ----
cache=os.path.join(IM,"CellMarker_Augmented_2021.gmt.json")
if os.path.exists(cache):
    lib=json.load(open(cache,encoding="utf-8"))
else:
    url="https://maayanlab.cloud/Enrichr/geneSetLibrary"
    r=requests.get(url,params={"mode":"text","libraryName":"CellMarker_Augmented_2021"},
                   proxies=PX,timeout=60)
    r.raise_for_status()
    lib={}
    for line in r.text.splitlines():
        p=line.split("\t")
        if len(p)>=3: lib[p[0]]=[g for g in p[2:] if g]
    json.dump(lib,open(cache,"w",encoding="utf-8"))
print("CellMarker terms:",len(lib))

# ---- 2. aggregate into broad immune lineages ----
def match(prefix):
    s=prefix.lower()
    if "natural killer" in s or s.startswith("nk") or " nk " in s: return "NK细胞"
    if "neutrophil" in s: return "中性粒细胞"
    if "mast" in s: return "肥大细胞"
    if "dendritic" in s: return "树突状细胞"
    if "monocyte" in s: return "单核细胞"
    if "macrophage" in s: return "巨噬细胞"
    if "plasma" in s or "plasmablast" in s: return "浆细胞"
    if "b cell" in s or "b lymphocyte" in s or s=="b cell": return "B细胞"
    if any(k in s for k in ["t cell","t lymphocyte","cd4","cd8","helper t","regulatory t",
                            "gamma delta","cytotoxic t","naive t","memory t","treg","t helper"]): return "T细胞"
    return None
raw={}; counts={}
for term,genes in lib.items():
    lin=match(term.split(":")[0])
    if not lin: continue
    raw.setdefault(lin,set()).update(genes)
    for g in genes: counts[(lin,g)]=counts.get((lin,g),0)+1
# keep markers recurring in >=2 independent CellMarker terms (robust core markers)
lineage={}
for lin,allg in raw.items():
    keep=sorted(g for g in allg if counts[(lin,g)]>=2)
    lineage[lin]=keep if len(keep)>=10 else sorted(allg)   # fallback if too strict
lineage={k:sorted(v) for k,v in lineage.items()}
for k,v in lineage.items(): print(f"  {k}: {len(v)} markers (recurrence>=2)")
assert len(lineage)>=8, "too few lineages"
json.dump(lineage,open(os.path.join(IM,"免疫标记集_已整理.json"),"w",encoding="utf-8"))

# ---- 3. expression matrix ----
with gzip.open(os.path.join(OUT,"GSE89632_series_matrix.txt.gz"),"rt",errors="replace") as f:
    lines=f.readlines()
def meta(p):
    for ln in lines:
        if ln.startswith(p): return [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
gsm=meta("!Sample_geo_accession"); diag=None
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
expr=pd.DataFrame(tb[1:],columns=tb[0]).set_index("ID_REF").apply(pd.to_numeric); expr.columns=gsm
pmap=pd.read_csv(os.path.join(OUT,"probe_gene_map.csv"),dtype=str).set_index("probe_id")
expr=expr.join(pmap[["symbol"]],how="left")
m=expr.symbol.notna() & ~expr.symbol.str.startswith("ERCC") & (expr.symbol.str.strip()!="")
expr=expr[m].copy(); expr["mean_all"]=expr[gsm].mean(axis=1)
expr=expr.sort_values("mean_all",ascending=False).drop_duplicates("symbol").set_index("symbol")
nash=[g for g,d in zip(gsm,diag) if d=="NASH"]; hc=[g for g,d in zip(gsm,diag) if d=="HC"]
print("samples NASH/HC:",len(nash),len(hc))

# ---- 4. ssGSEA via gseapy (validated implementation) ----
import gseapy as gp
ss=gp.ssgsea(expr[gsm], gene_sets=lineage, sample_norm_method="rank",
             no_plot=True, outdir=None, min_size=3, max_size=6000, verbose=False)
wide=ss.res2d.pivot(index="Term", columns="Name", values="NES")
scores=wide.reindex(columns=gsm)             # NES matrix, lineages x samples
print("ssGSEA score matrix:",scores.shape,"| lineages:",list(scores.index))
assert scores.shape[0]>=8 and scores.shape[1]==len(gsm)
scores.to_csv(os.path.join(IM,"ssGSEA_免疫细胞评分.csv"),encoding="utf-8-sig")

# ---- 5. NASH vs HC comparison ----
res=[]
for lin in scores.index:
    a=scores.loc[lin,nash].dropna(); b=scores.loc[lin,hc].dropna()
    t,p=stats.ttest_ind(a,b,equal_var=False)
    res.append({"免疫细胞":lin,"NASH均值":round(a.mean(),3),"HC均值":round(b.mean(),3),
                "差值":round(a.mean()-b.mean(),3),"t":round(t,2),"p_value":p})
R=pd.DataFrame(res); R["adj_p_BH"]=multipletests(R.p_value,method="fdr_bh")[1]
R=R.sort_values("p_value"); R.to_csv(os.path.join(IM,"免疫浸润_NASH_vs_HC.csv"),index=False,encoding="utf-8-sig")
print("\n",R.to_string(index=False))

# ---- 6. correlation of immune scores with core targets ----
core=["IL1B","TNF","PTGS2","IL6","IL1A"]
cors=[]
for lin in scores.index:
    for g in core:
        if g in expr.index:
            r,p=stats.spearmanr(scores.loc[lin,gsm],expr.loc[g,gsm])
            cors.append({"免疫细胞":lin,"靶点":g,"spearman_r":round(r,3),"p":round(p,4)})
C=pd.DataFrame(cors); C.to_csv(os.path.join(IM,"免疫评分_靶点相关性.csv"),index=False,encoding="utf-8-sig")

# ---- figures ----
fig,ax=plt.subplots(figsize=(11,6))
Rp=R.sort_values("差值")
cols=["#C0392B" if q<0.05 else "#95A5A6" for q in Rp.adj_p_BH]
ax.barh(Rp["免疫细胞"],Rp["差值"],color=cols)
ax.axvline(0,color="k",lw=0.8); ax.set_xlabel("ssGSEA 评分差值 (NASH - HC)")
ax.set_title("免疫细胞浸润在 NASH 与 HC 间的差异（红=adj p<0.05）",fontsize=13)
plt.tight_layout(); plt.savefig(os.path.join(IM,"fig_免疫差异.png"),dpi=200); plt.close()

# heatmap of scores
fig,ax=plt.subplots(figsize=(13,6))
Z=scores[gsm].values
Zm=(Z.T-np.nanmean(Z,axis=1)).T; Zs=np.divide(Zm,np.nanstd(Zm,axis=1,keepdims=True))
im=ax.imshow(Zs,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
ax.set_yticks(range(len(scores.index))); ax.set_yticklabels(scores.index)
ax.set_xticks(range(len(gsm))); ax.set_xticklabels(gsm,rotation=90,fontsize=6)
ax.set_title("各样本免疫细胞 ssGSEA 评分（标准化）",fontsize=13); fig.colorbar(im,ax=ax)
plt.tight_layout(); plt.savefig(os.path.join(IM,"fig_免疫热图.png"),dpi=200); plt.close()
print("\nfiles:",os.listdir(IM))
