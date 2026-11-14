# -*- coding: utf-8 -*-
"""转录组 Meta 整合：各队列 Hedges' g -> 随机效应（DerSimonian-Laird）合并，I² 异质性，森林图。"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

cohorts = {
    "GSE281797": np.log2(load_281797().astype(float) + 1),
    "GSE260666": np.log2(load_260666().astype(float) + 1),
    "GSE147304": np.log2(load_147304().astype(float) + 1),
}

# 各队列病例对照分组（疾病 vs 对照）
def groups_281797():
    md = META["GSE281797"]
    g0 = [s for s, v in md.items() if v["diagnosis"] == "No Pathology"]
    g1 = [s for s, v in md.items() if v["diagnosis"] in ("MASH", "MASL")]
    return g0, g1

def groups_260666():
    md = META["GSE260666"]
    g0 = [s for s, v in md.items() if v["group"] == "Control"]
    g1 = [s for s, v in md.items() if v["group"] == "NASH"]
    return g0, g1

def groups_147304():
    cols = list(cohorts["GSE147304"].columns)
    g1 = [c for c in cols if c.startswith("NASH")]
    g0 = [c for c in cols if not c.startswith("NASH")]
    return g0, g1

GROUP_FN = {
    "GSE281797": groups_281797,
    "GSE260666": groups_260666,
    "GSE147304": groups_147304,
}

def hedges_g(a, b):
    """a=疾病, b=对照（log2 尺度）。返回校正后 g 与方差。"""
    n1, n0 = len(a), len(b)
    s1, s0 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / (n1 + n0 - 2))
    if sp == 0:
        return np.nan, np.nan
    g = (np.mean(a) - np.mean(b)) / sp
    J = 1 - 3 / (4 * (n1 + n0) - 9)          # 小样本校正
    g = g * J
    v = (n1 + n0) / (n1 * n0) + g ** 2 / (2 * (n1 + n0))
    return float(g), float(v)

def random_effects(es):
    """es: list of (g, v)。返回合并效应、SE、p、Q、I²、tau²。"""
    es = [(g, v) for g, v in es if not np.isnan(g) and not np.isnan(v)]
    k = len(es)
    g = np.array([x[0] for x in es]); v = np.array([x[1] for x in es])
    w = 1 / v
    fe = np.sum(w * g) / np.sum(w)
    Q = np.sum(w * (g - fe) ** 2)
    df = k - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = 1 / (v + tau2)
    pooled = np.sum(wr * g) / np.sum(wr)
    se = 1 / np.sqrt(np.sum(wr))
    z = pooled / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return {"k": k, "pooled_g": pooled, "se": se, "z": z, "p": p,
            "Q": Q, "df": df, "I2": I2, "tau2": tau2,
            "g": g, "v": v, "wr": wr}

rows = []
detail = {}
for gene in TARGETS:
    es = []
    per = []
    for c, fn in GROUP_FN.items():
        lm = cohorts[c]
        if gene not in lm.index:
            per.append((c, np.nan, np.nan, 0, 0)); continue
        g0, g1 = fn()
        a = lm.loc[gene, [s for s in g1 if s in lm.columns]].dropna().values
        b = lm.loc[gene, [s for s in g0 if s in lm.columns]].dropna().values
        g, v = hedges_g(a, b)
        es.append((g, v)); per.append((c, g, v, len(a), len(b)))
    re = random_effects(es)
    detail[gene] = (per, re)
    rows.append({"gene": gene, "k": re["k"], "pooled_hedges_g": re["pooled_g"],
                 "se": re["se"], "z": re["z"], "p": re["p"],
                 "I2_percent": re["I2"], "Q": re["Q"], "tau2": re["tau2"],
                 "direction": "up" if re["pooled_g"] > 0 else "down"})

res = pd.DataFrame(rows).sort_values(["direction", "p"])
res.to_csv(os.path.join(OUTDIR, "Meta_随机效应合并结果.csv"), index=False, encoding="utf-8-sig")

# 逐基因效应明细长表
long = []
for gene, (per, re) in detail.items():
    for c, g, v, n1, n0 in per:
        long.append({"gene": gene, "cohort": c, "hedges_g": g, "variance": v,
                     "n_case": n1, "n_control": n0,
                     "ci_low": (g - 1.96 * np.sqrt(v)) if not np.isnan(g) else np.nan,
                     "ci_high": (g + 1.96 * np.sqrt(v)) if not np.isnan(g) else np.nan})
pd.DataFrame(long).to_csv(os.path.join(OUTDIR, "Meta_各队列效应量.csv"), index=False, encoding="utf-8-sig")

# 森林图（核心基因）
core = ["IL1B", "TNF", "PTGS2", "IL6", "CALCRL", "TACR1", "TRPV1"]
figdir = os.path.join(OUTDIR, "forest"); os.makedirs(figdir, exist_ok=True)
for gene in core:
    per, re = detail[gene]
    items = [(c, g, v) for c, g, v, n1, n0 in per if not np.isnan(g)]
    items.append(("Pooled (RE)", re["pooled_g"], re["se"] ** 2))
    y = np.arange(len(items))[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.7 * len(items) + 1.5))
    for yi, (c, g, v) in zip(y, items):
        se = np.sqrt(v)
        if c.startswith("Pooled"):
            ax.errorbar(g, yi, xerr=1.96 * se, fmt="D", color="#b03a2e", ecolor="#b03a2e", capsize=4, markersize=9)
        else:
            ax.errorbar(g, yi, xerr=1.96 * se, fmt="s", color="#2471a3", ecolor="#5499c7", capsize=3, markersize=7)
        ax.text(-0.02, yi, c, ha="right", va="center", transform=ax.get_yaxis_transform(), fontsize=9)
    ax.axvline(0, color="gray", lw=1, ls="--")
    ax.set_yticks([]); ax.set_xlabel("Hedges' g (95% CI)")
    ax.set_title(f"{gene} | pooled g={re['pooled_g']:.2f}, I²={re['I2']:.0f}%, p={re['p']:.3f}")
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(figdir, f"forest_{gene}.png"), dpi=150); plt.close(fig)

print("genes meta-analyzed:", len(res))
print(res[res["gene"].isin(core)].to_string(index=False))
