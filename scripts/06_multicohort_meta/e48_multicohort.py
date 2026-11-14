# -*- coding: utf-8 -*-
"""多队列独立验证：5 个独立人肝队列，病例对照 + 分级相关 + 跨队列投票。"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = "multicohort"
META = json(open(os.path.join(OUTDIR, "_meta", "cohort_metadata.json"), encoding="utf-8")) if False else json.load(
    open(os.path.join(OUTDIR, "_meta", "cohort_metadata.json"), encoding="utf-8"))

# 靶基因：47 交集 + 痛觉基因
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
    f = os.path.join(OUTDIR, "GSE281797", "GSE281797_tpm.txt.gz")
    d = pd.read_csv(f, sep="\t")
    d = d.rename(columns={"HGNC symbol": "symbol"})
    keep = ["symbol"] + [c for c in d.columns if c.startswith("S") and len(c) <= 5]
    return dedupe(d[keep], "symbol")

def load_260666():
    f = os.path.join(OUTDIR, "GSE260666", "GSE260666_raw_counts.txt.gz")
    d = pd.read_csv(f, sep="\t").rename(columns={"Symbol": "symbol"})
    m = dedupe(d, "symbol")
    # counts -> CPM
    m = m.div(m.sum(axis=0), axis=1) * 1e6
    return m

def load_147304():
    f = os.path.join(OUTDIR, "GSE147304", "GSE147304_gene_expression.txt.gz")
    d = pd.read_csv(f, sep="\t")
    fpkm = [c for c in d.columns if c.endswith("_FPKM")]
    d = d[["gene_name"] + fpkm].rename(columns={"gene_name": "symbol"})
    d.columns = ["symbol"] + [c[:-5] for c in fpkm]
    return dedupe(d, "symbol")

def load_gct(gse, fname):
    f = os.path.join(OUTDIR, gse, fname)
    d = pd.read_csv(f, sep="\t", skiprows=2)
    d = d.rename(columns={d.columns[0]: "symbol"})
    sample_cols = [c for c in d.columns if c not in ("symbol", "gene_id")]
    return dedupe(d[["symbol"] + sample_cols], "symbol")

def to_log2(m):
    return np.log2(m.astype(float) + 1.0)

# 载入各队列（线性归一化尺度）
cohorts = {
    "GSE281797": load_281797(),
    "GSE260666": load_260666(),
    "GSE147304": load_147304(),
    "GSE193066": load_gct("GSE193066", "GSE193066_NAFLD.HUn164.gct.gz"),
    "GSE192959": load_gct("GSE192959", "GSE192959_NAFLD.TUn48.gct.gz"),
}
logm = {k: to_log2(v) for k, v in cohorts.items()}

case_rows = []   # 病例对照结果
corr_rows = []   # 分级相关结果

def case_test(cohort, comparison, g0, g1):
    lm = logm[cohort]
    g0 = [s for s in g0 if s in lm.columns]
    g1 = [s for s in g1 if s in lm.columns]
    for gene in TARGETS:
        if gene not in lm.index:
            case_rows.append({"cohort": cohort, "comparison": comparison, "gene": gene,
                              "n0": len(g0), "n1": len(g1), "log2FC": np.nan,
                              "welch_p": np.nan, "mw_p": np.nan, "present": False})
            continue
        a = lm.loc[gene, g1].dropna().values
        b = lm.loc[gene, g0].dropna().values
        lfc = float(np.mean(a) - np.mean(b))
        med_lfc = float(np.median(a) - np.median(b))
        wp = stats.ttest_ind(a, b, equal_var=False).pvalue if len(a) > 1 and len(b) > 1 else np.nan
        mp = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) >= 1 and len(b) >= 1 else np.nan
        # 离群标注：对照组存在明显高于其余的极端值（用 1.5*IQR 或最大值/中位数比）
        outlier = False
        if len(b) >= 4:
            q1, q3 = np.percentile(b, [25, 75])
            iqr = q3 - q1
            hi = q3 + 1.5 * iqr
            outlier = bool(np.any(b > hi) and np.max(b) > np.median(b) + 1.0)
        case_rows.append({"cohort": cohort, "comparison": comparison, "gene": gene,
                          "n0": len(b), "n1": len(a), "log2FC": lfc,
                          "median_log2FC": med_lfc, "control_outlier": outlier,
                          "welch_p": wp, "mw_p": mp, "present": True})

def corr_test(cohort, grade_key, label):
    lm = logm[cohort]
    md = META[cohort]
    samples = [s for s in lm.columns if s in md and md[s].get(grade_key) is not None]
    grades = np.array([md[s][grade_key] for s in samples], dtype=float)
    for gene in TARGETS:
        if gene not in lm.index:
            corr_rows.append({"cohort": cohort, "grade": label, "gene": gene,
                              "n": len(samples), "spearman_rho": np.nan,
                              "spearman_p": np.nan, "present": False})
            continue
        x = lm.loc[gene, samples].dropna().astype(float).values
        if len(x) != len(grades) or np.std(x) == 0:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(x, grades)
        corr_rows.append({"cohort": cohort, "grade": label, "gene": gene,
                          "n": len(samples), "spearman_rho": float(rho) if rho == rho else np.nan,
                          "spearman_p": float(p) if p == p else np.nan, "present": True})

# ---- GSE281797 ----
md = META["GSE281797"]
nop = [s for s, v in md.items() if v["diagnosis"] == "No Pathology"]
mash = [s for s, v in md.items() if v["diagnosis"] == "MASH"]
dis = [s for s, v in md.items() if v["diagnosis"] in ("MASH", "MASL")]
case_test("GSE281797", "MASH vs NoPath", nop, mash)
case_test("GSE281797", "Disease vs NoPath", nop, dis)
corr_test("GSE281797", "inflammation", "inflammation grade")
corr_test("GSE281797", "fibrosis", "fibrosis grade")

# ---- GSE260666 ----
md = META["GSE260666"]
ctrl = [s for s, v in md.items() if v["group"] == "Control"]
nafl = [s for s, v in md.items() if v["group"] == "NAFL"]
nash = [s for s, v in md.items() if v["group"] == "NASH"]
case_test("GSE260666", "NASH vs Control", ctrl, nash)
case_test("GSE260666", "NAFL vs Control", ctrl, nafl)
# 有序趋势
lm = logm["GSE260666"]
order = {"Control": 0, "NAFL": 1, "NASH": 2}
samples = [s for s in lm.columns if s in md]
grp = np.array([order[md[s]["group"]] for s in samples])
for gene in TARGETS:
    if gene not in lm.index:
        continue
    x = lm.loc[gene, samples].astype(float).values
    rho, p = stats.spearmanr(x, grp) if np.std(x) > 0 else (np.nan, np.nan)
    corr_rows.append({"cohort": "GSE260666", "grade": "group order (Ctrl<NAFL<NASH)",
                      "gene": gene, "n": len(samples), "spearman_rho": rho,
                      "spearman_p": p, "present": True})

# ---- GSE147304 ----
lm = logm["GSE147304"]
nash_cols = [c for c in lm.columns if c.startswith("NASH")]
norm_cols = [c for c in lm.columns if not c.startswith("NASH")]
case_test("GSE147304", "NASH vs Normal", norm_cols, nash_cols)

# ---- GSE193066 / GSE192959 ----
for c in ("GSE193066", "GSE192959"):
    corr_test(c, "fibrosis", "fibrosis stage")
    corr_test(c, "nas", "NAFLD activity score")

case_df = pd.DataFrame(case_rows)
corr_df = pd.DataFrame(corr_rows)

# BH 校正（在每个比较内、对 present 基因）
def bh(pvals):
    """标准 Benjamini-Hochberg（升序排，保证单调）。"""
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    valid = ~np.isnan(p)
    n = int(valid.sum())
    if n == 0:
        return out
    order = np.argsort(np.where(valid, p, np.inf), kind="stable")
    adj = np.empty(n)
    prev = 1.0
    for k in range(n - 1, -1, -1):
        idx = order[k]              # 第 k+1 小（1-based rank = k+1）
        rank = k + 1
        val = min(prev, p[idx] * n / rank)
        adj[k] = val
        prev = val
    for k in range(n):
        out[order[k]] = adj[k]
    return out

case_df["welch_adj"] = np.nan
for (c, cmp), sub in case_df.groupby(["cohort", "comparison"]):
    idx = sub.index
    case_df.loc[idx, "welch_adj"] = bh(sub["welch_p"].values)

corr_df["spearman_adj"] = np.nan
for (c, g), sub in corr_df.groupby(["cohort", "grade"]):
    idx = sub.index
    corr_df.loc[idx, "spearman_adj"] = bh(sub["spearman_p"].values)

case_df.to_csv(os.path.join(OUTDIR, "多队列_病例对照结果.csv"), index=False, encoding="utf-8-sig")
corr_df.to_csv(os.path.join(OUTDIR, "多队列_分级相关结果.csv"), index=False, encoding="utf-8-sig")

# ---- 跨队列投票：每个基因在独立比较中的方向 ----
# 病例对照：以 5 个独立“主比较”各算一票（取每个队列最主要的疾病对照比较）
main_case = {
    "GSE281797": "Disease vs NoPath",
    "GSE260666": "NASH vs Control",
}
INCONCLUSIVE_CASE = ("GSE147304", "NASH vs Normal")
vote = {}
for gene in TARGETS:
    signs = []
    detail = []
    incon = []
    for c, cmp in main_case.items():
        r = case_df[(case_df["cohort"] == c) & (case_df["comparison"] == cmp) & (case_df["gene"] == gene)]
        if len(r) and r.iloc[0]["present"] and not np.isnan(r.iloc[0]["log2FC"]):
            s = 1 if r.iloc[0]["log2FC"] > 0 else -1
            signs.append(s); detail.append(f"{c}:{'↑' if s>0 else '↓'}")
    # GSE147304：样本小、对照特殊、含离群，标记为证据不足，不计入一致率
    c, cmp = INCONCLUSIVE_CASE
    r = case_df[(case_df["cohort"] == c) & (case_df["comparison"] == cmp) & (case_df["gene"] == gene)]
    if len(r) and r.iloc[0]["present"]:
        s = 1 if r.iloc[0]["log2FC"] > 0 else -1
        incon.append(f"{c}:{'↑' if s>0 else '↓'}(inconclusive,n=5)")
    # 分级相关：4 个队列的纤维化 + 炎症/NAS（每个队列取一个代表分级）
    main_corr = {
        "GSE281797": "inflammation grade",
        "GSE260666": "group order (Ctrl<NAFL<NASH)",
        "GSE193066": "fibrosis stage",
        "GSE192959": "fibrosis stage",
    }
    for c, g in main_corr.items():
        r = corr_df[(corr_df["cohort"] == c) & (corr_df["grade"] == g) & (corr_df["gene"] == gene)]
        if len(r) and r.iloc[0]["present"] and not np.isnan(r.iloc[0]["spearman_rho"]):
            s = 1 if r.iloc[0]["spearman_rho"] > 0 else -1
            signs.append(s); detail.append(f"{c}:{'+' if s>0 else '-'}")
    up = sum(1 for s in signs if s > 0)
    down = len(signs) - up
    vote[gene] = {"n_votes": len(signs), "up_positive": up, "down_negative": down,
                  "concordance": (max(up, down) / len(signs) if signs else np.nan),
                  "direction": ("up" if up > down else "down" if down > up else "tie"),
                  "detail": "; ".join(detail),
                  "inconclusive": "; ".join(incon)}

vote_df = pd.DataFrame([dict(gene=g, **v) for g, v in vote.items()])
vote_df = vote_df.sort_values(["direction", "concordance", "up_positive"], ascending=[True, False, False])
vote_df.to_csv(os.path.join(OUTDIR, "多队列_跨队列方向投票.csv"), index=False, encoding="utf-8-sig")

print("cohorts loaded:", {k: v.shape for k, v in cohorts.items()})
print("case rows:", len(case_df), "| corr rows:", len(corr_df), "| vote genes:", len(vote_df))
print("\n=== 核心基因跨队列投票 ===")
core = ["PTGS2", "IL6", "IL1B", "IL1A", "PTGES", "PTGER3", "TNF",
        "TRPV1", "TRPA1", "TACR1", "CALCRL", "NGF", "SCN9A"]
for gene in core:
    if gene in vote:
        v = vote[gene]
        print(f"{gene:>7} | {v['direction']:>4} | up {v['up_positive']}/{v['n_votes']} "
              f"| concord {v['concordance']:.2f} | {v['detail']}")
