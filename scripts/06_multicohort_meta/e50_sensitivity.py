# -*- coding: utf-8 -*-
"""Meta 敏感性分析：不同纳入口径下的合并 Hedges g。"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = "multicohort"
META = json.load(open(os.path.join(OUTDIR, "_meta", "cohort_metadata.json"), encoding="utf-8"))
inter = pd.read_csv("交集_DEG_成分靶点.csv")
TARGETS = list(dict.fromkeys(inter["gene"].tolist() +
                             ["TNF", "TRPV1", "TRPA1", "TACR1", "CALCRL", "NGF", "SCN9A"]))

def dedupe(df, symbol_col):
    df[symbol_col] = df[symbol_col].astype(str).str.strip()
    df = df[df[symbol_col].str.len() > 0]
    num_cols = [c for c in df.columns if c != symbol_col]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")
    return df.groupby(symbol_col, sort=False)[num_cols].mean()

def load_281797():
    d = pd.read_csv(os.path.join(OUTDIR, "GSE281797", "GSE281797_tpm.txt.gz"), sep="\t")
    d = d.rename(columns={"HGNC symbol": "symbol"})
    keep = ["symbol"] + [c for c in d.columns if c.startswith("S") and len(c) <= 5]
    return dedupe(d[keep], "symbol")

def load_260666():
    d = pd.read_csv(os.path.join(OUTDIR, "GSE260666", "GSE260666_raw_counts.txt.gz"), sep="\t")
    d = d.rename(columns={"Symbol": "symbol"})
    m = dedupe(d, "symbol")
    return m.div(m.sum(axis=0), axis=1) * 1e6

def load_147304():
    d = pd.read_csv(os.path.join(OUTDIR, "GSE147304", "GSE147304_gene_expression.txt.gz"), sep="\t")
    fpkm = [c for c in d.columns if c.endswith("_FPKM")]
    d = d[["gene_name"] + fpkm].rename(columns={"gene_name": "symbol"})
    d.columns = ["symbol"] + [c[:-5] for c in fpkm]
    return dedupe(d, "symbol")

logm = {
    "GSE281797": np.log2(load_281797().astype(float) + 1),
    "GSE260666": np.log2(load_260666().astype(float) + 1),
    "GSE147304": np.log2(load_147304().astype(float) + 1),
}

# 各口径分组定义
def g_281797_disease():
    md = META["GSE281797"]
    g0 = [s for s, v in md.items() if v["diagnosis"] == "No Pathology"]
    g1 = [s for s, v in md.items() if v["diagnosis"] in ("MASH", "MASL")]
    return g0, g1

def g_281797_mash():
    md = META["GSE281797"]
    g0 = [s for s, v in md.items() if v["diagnosis"] == "No Pathology"]
    g1 = [s for s, v in md.items() if v["diagnosis"] == "MASH"]
    return g0, g1

def g_260666():
    md = META["GSE260666"]
    return ([s for s, v in md.items() if v["group"] == "Control"],
            [s for s, v in md.items() if v["group"] == "NASH"])

def g_147304():
    cols = list(logm["GSE147304"].columns)
    return ([c for c in cols if not c.startswith("NASH")],
            [c for c in cols if c.startswith("NASH")])

SCENARIOS = {
    "S1_all3": [("GSE281797", g_281797_disease), ("GSE260666", g_260666), ("GSE147304", g_147304)],
    "S2_exclude147304": [("GSE281797", g_281797_disease), ("GSE260666", g_260666)],
    "S3_MASH_severe": [("GSE281797", g_281797_mash), ("GSE260666", g_260666)],
}

def hedges_g(a, b):
    n1, n0 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n0 - 1) * np.var(b, ddof=1)) / (n1 + n0 - 2))
    if sp == 0:
        return np.nan, np.nan
    g = (np.mean(a) - np.mean(b)) / sp * (1 - 3 / (4 * (n1 + n0) - 9))
    v = (n1 + n0) / (n1 * n0) + g ** 2 / (2 * (n1 + n0))
    return float(g), float(v)

def pool(es):
    es = [(g, v) for g, v in es if not np.isnan(g)]
    if not es:
        return np.nan, np.nan, np.nan
    g = np.array([x[0] for x in es]); v = np.array([x[1] for x in es]); w = 1 / v
    fe = np.sum(w * g) / np.sum(w)
    Q = np.sum(w * (g - fe) ** 2); df = len(es) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = 1 / (v + tau2)
    pg = np.sum(wr * g) / np.sum(wr); se = 1 / np.sqrt(np.sum(wr))
    p = 2 * (1 - stats.norm.cdf(abs(pg / se)))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return pg, p, I2

rows = []
for gene in TARGETS:
    for sc, defs in SCENARIOS.items():
        es = []
        for c, fn in defs:
            lm = logm[c]
            if gene not in lm.index:
                continue
            g0, g1 = fn()
            a = lm.loc[gene, [s for s in g1 if s in lm.columns]].dropna().values
            b = lm.loc[gene, [s for s in g0 if s in lm.columns]].dropna().values
            es.append(hedges_g(a, b))
        pg, p, I2 = pool(es)
        rows.append({"gene": gene, "scenario": sc, "k": len(es),
                     "pooled_g": pg, "p": p, "I2": I2})

sens = pd.DataFrame(rows)
sens.to_csv(os.path.join(OUTDIR, "Meta_敏感性分析.csv"), index=False, encoding="utf-8-sig")
core = ["IL1B", "TNF", "PTGS2", "IL6", "CALCRL", "TACR1", "TRPV1"]
print(sens[sens.gene.isin(core)].pivot(index="gene", columns="scenario",
                                        values="pooled_g").round(3).to_string())
print("\n--- p values ---")
print(sens[sens.gene.isin(core)].pivot(index="gene", columns="scenario",
                                        values="p").round(3).to_string())
