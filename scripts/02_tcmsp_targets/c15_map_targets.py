# -*- coding: utf-8 -*-
import os, re, pandas as pd
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
u = pd.read_csv(os.path.join(OUT,"uniprot_human_reviewed.tsv"), sep="\t", dtype=str)
T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv"))

def norm(s):
    s = s.lower()
    s = re.sub(r"\(ec [\d\.\- ]+\)", " ", s)      # drop EC numbers
    s = re.sub(r"[^a-z0-9]+", " ", s).strip()
    return re.sub(r"\s+", " ", s)

def recommended(s):
    # first chunk before first parenthesis, also strip EC tail
    s = re.sub(r"\(ec [\d\.\- ]+\)", " ", s)
    return s.split("(")[0].strip()

# build lookups
exact, rec = {}, {}
for _, r in u.iterrows():
    gene = r["Gene Names (primary)"]
    if not isinstance(gene, str) or not gene.strip(): continue
    pn = r["Protein names"]
    exact.setdefault(norm(pn), gene)
    rec.setdefault(norm(recommended(pn)), gene)
    # also index every parenthesized alternative name
    for alt in re.findall(r"\(([^()]+)\)", pn):
        if not re.search(r"ec ", alt.lower()):
            rec.setdefault(norm(alt), gene)

targets = sorted(T["target_name"].dropna().unique())
mapping = {}
for tname in targets:
    n = norm(tname)
    g = exact.get(n) or rec.get(n) or rec.get(norm(recommended(tname)))
    mapping[tname] = g

mapped = {k: v for k, v in mapping.items() if v}
unmapped = [k for k, v in mapping.items() if not v]
print(f"total targets: {len(targets)} | mapped to human gene: {len(mapped)} | unmapped: {len(unmapped)}")
print("\n--- sample unmapped (first 30) ---")
for x in unmapped[:30]: print("  ", x)
# save mapping
pd.DataFrame([{"target_name": k, "gene": v} for k, v in mapping.items()]) \
  .to_csv(os.path.join(OUT,"target_gene_mapping.csv"), index=False, encoding="utf-8-sig")
# verify pain targets
for k in ["Transient receptor potential cation channel subfamily V member 1",
          "Prostaglandin G/H synthase 2","Tumor necrosis factor","Interleukin-1 beta",
          "Interleukin-6","Sodium channel protein type 9 subunit alpha"]:
    print("CHECK", k, "->", mapping.get(k))
