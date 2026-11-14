# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd, gseapy as gp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["SimHei"]; plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
EN = os.path.join(OUT, "enrichment"); os.makedirs(EN, exist_ok=True)

# ---- gene sets ----
M = pd.read_csv(os.path.join(OUT,"TCMSP_全方成分主表.csv"))
T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv"))
mp = pd.read_csv(os.path.join(OUT,"target_gene_mapping.csv")).dropna(subset=["gene"])
T["gene"] = T["target_name"].map(dict(zip(mp.target_name, mp.gene)))
fmol = set(M[M.pass_filter].MOL_ID)
formula237 = sorted(set(T[T.MOL_ID.isin(fmol)].gene.dropna()))
deg = pd.read_csv(os.path.join(OUT,"GSE89632_DEG_NASH_vs_HC.csv"))
sig = deg[(deg.adj_p_BH<0.05)&(deg.logFC_NASH_vs_HC.abs()>=0.5)]
inter47 = sorted(set(sig.gene) & set(formula237))
print("formula237:", len(formula237), "| inter47:", len(inter47))

libs_all = gp.get_library_name(organism="Human")
kegg = next((l for l in ["KEGG_2021_Human","KEGG_2019_Human","KEGG_2016"] if l in libs_all), None)
LIBS = {"GO_BP":"GO_Biological_Process_2023","GO_CC":"GO_Cellular_Component_2023",
        "GO_MF":"GO_Molecular_Function_2023","KEGG":kegg,"Reactome":"Reactome_2022"}
LIBS = {k:v for k,v in LIBS.items() if v}
print("using libraries:", LIBS)

def run_enrich(genes, tag):
    outdir = os.path.join(EN, tag); os.makedirs(outdir, exist_ok=True)
    summary = []
    for short, lib in LIBS.items():
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=lib, organism="human",
                             outdir=None, background=20000, cutoff=0.05, no_plot=True)
            r = enr.results
            r.to_csv(os.path.join(outdir, f"{short}_{lib}.csv"), index=False, encoding="utf-8-sig")
            sig_r = r[r["Adjusted P-value"]<0.05]
            summary.append((short, lib, len(sig_r)))
        except Exception as e:
            print("  ERR", tag, lib, repr(e)[:100])
    print(f"[{tag}] significant terms:", summary)
    return outdir

d1 = run_enrich(formula237, "formula237")
d2 = run_enrich(inter47, "inter47")

# ---- custom plots: top 12 by adj p for each lib (formula237) ----
def plot_top(libshort, tag, ax):
    folder = os.path.join(EN, tag)
    fs = [f for f in os.listdir(folder) if f.startswith(libshort)]
    if not fs: return False
    r = pd.read_csv(os.path.join(folder, fs[0]))
    r = r[r["Adjusted P-value"]<0.05].sort_values("Adjusted P-value").head(12)
    if not len(r): return False
    y = np.arange(len(r))[::-1]
    ax.barh(y, -np.log10(r["Adjusted P-value"]), color="#4A90D9")
    ax.set_yticks(y); ax.set_yticklabels([t[:62] for t in r.Term], fontsize=8)
    ax.set_xlabel("-log10(adj P)"); ax.set_title(f"{libshort} Top terms", fontsize=11)
    return True

fig, axs = plt.subplots(2, 2, figsize=(15, 11))
for ax, k in zip(axs.flat, ["GO_BP","GO_MF","KEGG","Reactome"]):
    if not plot_top(k, "formula237", ax): ax.set_visible(False)
fig.suptitle("成分靶点功能富集（237 个过筛成分人类基因）", fontsize=14)
plt.tight_layout(); plt.savefig(os.path.join(EN,"fig_富集_formula237.png"), dpi=200); plt.close()

# GO CC separate
fig, ax = plt.subplots(figsize=(11,6))
if not plot_top("GO_CC","formula237",ax): ax.set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(EN,"fig_富集_GOCC.png"), dpi=200); plt.close()
print("\nfiles:", os.listdir(EN))
