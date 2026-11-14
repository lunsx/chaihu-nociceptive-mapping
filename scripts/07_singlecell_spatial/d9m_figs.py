# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\singlecell"
df=pd.read_csv(base+r"\靶点细胞定位_长表.csv")
targets=["Il1b","Tnf","Ptgs2","Il6","Il1a","Ptges","Ptger3","Trpv1","Trpa1"]
cts=["Kupffer cell","macrophage","monocyte","dendritic cell","neutrophil",
     "endothelial cell","cholangiocyte","hepatocyte","B cell","CD4-positive, alpha-beta T cell",
     "CD8-positive, alpha-beta T cell","natural killer cell","interneuron","neuron associated cell","hematopoietic cell"]
cts=[c for c in cts if c in set(df["cell_type"])]

def pivot(dis,col):
    sub=df[df["disease"]==dis]
    return sub.pivot(index="gene",columns="cell_type",values=col).reindex(index=targets,columns=cts)

# dot plot for MASH
meanM=pivot("metabolic dysfunction-associated steatohepatitis","mean_expr")
pctM=pivot("metabolic dysfunction-associated steatohepatitis","pct_expressing")
fig,ax=plt.subplots(figsize=(11,5.2))
for i,g in enumerate(targets):
    for j,c in enumerate(cts):
        v=meanM.iloc[i,j]; pc=pctM.iloc[i,j]
        if pd.isna(v): continue
        ax.scatter(j,i,s=10+pc*5,c=v,cmap="magma",vmin=0,vmax=3.6,edgecolors="gray",linewidths=0.3)
        if pc>5: ax.text(j,i,f"{pc:.0f}",ha="center",va="center",fontsize=7,color="white" if v>1.6 else "black")
ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts,rotation=40,ha="right",fontsize=8)
ax.set_yticks(range(len(targets))); ax.set_yticklabels(targets)
ax.set_xlim(-0.6,len(cts)-0.4); ax.set_ylim(len(targets)-0.5,-0.5)
sm=plt.cm.ScalarMappable(cmap="magma",norm=plt.Normalize(0,3.6)); sm.set_array([])
cb=fig.colorbar(sm,ax=ax,shrink=0.7); cb.set_label("mean expression (log-normalized)")
ax.set_title("Target gene localization in MASH liver (size/number = % expressing cells, color = mean)")
plt.tight_layout(); plt.savefig(base+r"\fig_单细胞定位点图.png",dpi=160); plt.close()

# MASH vs normal comparison
meanN=pivot("normal","mean_expr"); pctN=pivot("normal","pct_expressing")
rows=[]
for g in targets:
    for c in cts:
        a=meanM.loc[g,c]; b=meanN.loc[g,c]
        pa=pctM.loc[g,c]; pb=pctN.loc[g,c]
        if pd.isna(a) or pd.isna(b): continue
        rows.append({"gene":g,"cell_type":c,"mean_MASH":a,"mean_normal":b,
                     "delta_mean":round(a-b,3),"pct_MASH":pa,"pct_normal":pb,"delta_pct":round(pa-pb,2)})
cmp=pd.DataFrame(rows).sort_values("delta_mean",ascending=False)
cmp.to_csv(base+r"\MASH_vs_normal_对比.csv",index=False,encoding="utf-8-sig")
print("top increases:"); print(cmp.head(10).to_string(index=False))

# composition plot
comp=pd.read_csv(base+r"\细胞类型组成.csv")
piv=comp.pivot(index="cell_type",columns="disease",values="pct").fillna(0)
piv=piv.reindex(cts)
fig,ax=plt.subplots(figsize=(11,5))
x=np.arange(len(cts)); w=0.27
for k,d in enumerate(["normal","metabolic dysfunction-associated steatohepatitis","hepatocellular carcinoma"]):
    if d in piv.columns:
        ax.bar(x+(k-1)*w,piv[d].values,w,label=d.replace("metabolic dysfunction-associated steatohepatitis","MASH"))
ax.set_xticks(x); ax.set_xticklabels(cts,rotation=40,ha="right",fontsize=8)
ax.set_ylabel("% of cells"); ax.legend(); ax.set_title("Cell-type composition by disease group")
plt.tight_layout(); plt.savefig(base+r"\fig_细胞组成.png",dpi=160); plt.close()
print("\nfigures + comparison saved")
