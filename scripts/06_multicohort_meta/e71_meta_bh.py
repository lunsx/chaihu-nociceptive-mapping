# -*- coding: utf-8 -*-
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np, os
R = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
m = pd.read_csv(os.path.join(R,"multicohort","Meta_随机效应合并结果.csv")).copy()
p = m.p.values; n = len(p); order = np.argsort(p); ranked = p[order]
adj = np.empty(n); prev = 1.0
for i in range(n-1,-1,-1):
    prev = min(prev, ranked[i]*n/(i+1)); adj[i] = prev
bh = np.empty(n); bh[order] = adj
m["p_bh54"] = bh
m.to_csv(os.path.join(R,"multicohort","Meta_随机效应合并_BH校正.csv"), index=False, encoding="utf-8-sig")

s = pd.read_csv(os.path.join(R,"multicohort","Meta_敏感性分析.csv"))
print("敏感性 scenario 列唯一值:", s.scenario.unique())
print("敏感性 p<0.05 数量(原始):", (s.p<0.05).sum(), "/", len(s))
# 敏感性内逐口径BH
for sc, g in s.groupby("scenario"):
    gg = g.sort_values("p")
    print(sc, "最小p:", round(gg.p.min(),4), "基因数:", len(g))
