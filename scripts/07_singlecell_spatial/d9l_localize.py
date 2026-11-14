# -*- coding: utf-8 -*-
import anndata as ad, numpy as np, pandas as pd, csv
p=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\singlecell\MASH_liver.h5ad"
a=ad.read_h5ad(p)
X=a.X
# check normalization
sample=X[:200].toarray()
print("X max:",round(float(sample.max()),3),"mean:",round(float(sample.mean()),4))
# gene name -> column
names=list(a.var["feature_name"])
g2i={g:i for i,g in enumerate(names)}
targets=["Il1b","Tnf","Ptgs2","Il6","Il1a","Ptges","Ptger3","Trpv1","Trpa1"]
obs=a.obs.reset_index(drop=True)
celltypes=list(obs["cell_type"].cat.categories) if hasattr(obs["cell_type"],"cat") else sorted(obs["cell_type"].unique())
diseases=["normal","metabolic dysfunction-associated steatohepatitis","hepatocellular carcinoma"]
rows=[]
dense_idx=np.arange(a.n_obs)
for g in targets:
    col=np.asarray(X[:,g2i[g]].todense()).ravel()
    for ct in celltypes:
        for dis in diseases:
            m=(obs["cell_type"].values==ct)&(obs["disease"].values==dis)
            n=int(m.sum())
            if n==0: continue
            vals=col[m]
            pct=float((vals>0).mean()*100)
            mean=float(vals.mean())
            rows.append({"gene":g,"cell_type":ct,"disease":dis,"n_cells":n,
                         "pct_expressing":round(pct,2),"mean_expr":round(mean,4)})
out=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\singlecell\靶点细胞定位_长表.csv"
pd.DataFrame(rows).to_csv(out,index=False,encoding="utf-8-sig")
print("rows:",len(rows),"->",out)

# cell type composition by disease
comp=[]
for dis in diseases:
    sub=obs[obs["disease"]==dis]
    tot=len(sub)
    for ct in celltypes:
        n=int((sub["cell_type"]==ct).sum())
        comp.append({"disease":dis,"cell_type":ct,"n":n,"pct":round(n/tot*100,2)})
cout=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\singlecell\细胞类型组成.csv"
pd.DataFrame(comp).to_csv(cout,index=False,encoding="utf-8-sig")
print("composition ->",cout)

# quick console: MASH localization top pct
df=pd.DataFrame(rows)
mash=df[df["disease"]=="metabolic dysfunction-associated steatohepatitis"]
for g in targets:
    sub=mash[mash["gene"]==g].sort_values("pct_expressing",ascending=False).head(3)
    print("\n",g,"(MASH top cell types):")
    for _,r in sub.iterrows(): print("  ",r["cell_type"],r["pct_expressing"],"% mean",r["mean_expr"])
