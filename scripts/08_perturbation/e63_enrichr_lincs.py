# -*- coding: utf-8 -*-
import warnings; warnings.filterwarnings("ignore")
import requests, pandas as pd, json, os
PROX = {"http": "http://127.0.0.1:8902", "https": "http://127.0.0.1:8902"}
BASE = "https://maayanlab.cloud/Enrichr"

# 1. 构建疾病上调签名
vote = pd.read_csv(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\multicohort\多队列_跨队列方向投票.csv")
up = vote[(vote.direction == "up") & (vote.concordance >= 0.6)]["gene"].tolist()
print("disease up signature size:", len(up))
print(up)

# 2. 提交基因列表
r = requests.post(BASE + "/addList",
                  files={"list": (None, "\n".join(up)),
                         "description": (None, "liver inflammation up signature")},
                  proxies=PROX, timeout=60)
ulid = r.json()["userListId"]
print("userListId:", ulid)

libs = ["LINCS_L1000_Chem_Pert_down", "LINCS_L1000_Chem_Pert_up",
        "LINCS_L1000_Chem_Pert_Consensus_Sigs", "LINCS_L1000_CRISPR_KO_Consensus_Sigs"]
allres = {}
for lib in libs:
    e = requests.get(BASE + "/enrich",
                     params={"userListId": ulid, "backgroundType": lib},
                     proxies=PROX, timeout=120)
    data = e.json()[lib]
    df = pd.DataFrame(data, columns=["rank", "term", "pvalue", "zscore",
                                     "combined_score", "overlap_genes",
                                     "adjusted_pvalue", "old_pvalue", "old_adj"])
    df = df[["rank", "term", "pvalue", "adjusted_pvalue", "combined_score", "overlap_genes"]]
    allres[lib] = df
    df.to_csv(f"Enrichr_{lib}.csv", index=False, encoding="utf-8-sig")
    print(f"\n=== {lib}: {len(df)} terms, significant(adj<0.05):", (df.adjusted_pvalue < 0.05).sum())
    print(df.head(8).to_string(index=False))

# 3. 检索我们的化合物
print("\n\n=== TARGET COMPOUND SEARCH ===")
for kw in ["luteolin", "quercetin", "kaempferol", "paeoniflorin", "cyper", "flavone"]:
    for lib, df in allres.items():
        m = df[df.term.str.lower().str.contains(kw)]
        if len(m):
            print(f"\n[{kw}] in {lib}")
            print(m.head(6).to_string(index=False))
