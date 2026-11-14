# -*- coding: utf-8 -*-
import os, re, pandas as pd
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
u = pd.read_csv(os.path.join(OUT,"uniprot_human_reviewed.tsv"), sep="\t", dtype=str)
T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv"))

def toks(s):
    s = re.sub(r"\(ec [\d\.\- ]+\)", " ", s.lower())
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return [w for w in s.split() if w]

def variants(pn):
    out = [pn.split("(")[0]]
    for alt in re.findall(r"\(([^()]+)\)", pn):
        if not re.search(r"ec ", alt.lower()): out.append(alt)
    return out

# index all human name token sets (frozenset and multiset tuple)
name_index = []   # (token_tuple, gene)
for _, r in u.iterrows():
    gene = r["Gene Names (primary)"]
    if not isinstance(gene, str) or not gene.strip(): continue
    for v in variants(r["Protein names"]):
        tt = tuple(sorted(toks(v)))
        if tt: name_index.append((tt, gene))

def match(tname):
    q = tuple(sorted(toks(tname)))
    qset = set(q)
    # 1) exact multiset
    for tt, gene in name_index:
        if tt == q: return gene, "exact_token"
    # 2) same token set ignoring order/duplicates, high overlap
    best = None; bestj = 0.0
    for tt, gene in name_index:
        ts = set(tt)
        if not ts: continue
        j = len(qset & ts) / len(qset | ts)
        # require all query tokens covered and strong jaccard
        if j > bestj and qset.issubset(ts) and len(qset) >= 3:
            bestj, best = j, gene
    if best and bestj >= 0.85:
        return best, f"fuzzy_j={bestj:.2f}"
    return None, "unmapped"

targets = sorted(T["target_name"].dropna().unique())
res = {}
for tname in targets:
    res[tname] = match(tname)
m = {k: v[0] for k, v in res.items() if v[0]}
conf = {k: v[1] for k, v in res.items()}
print(f"total {len(targets)} | mapped {len(m)} | unmapped {len(targets)-len(m)}")
print("fuzzy matches:", sum(1 for c in conf.values() if c.startswith("fuzzy")))
for k in ["Sodium channel protein type 9 subunit alpha",
          "5-hydroxytryptamine 2A receptor","5-hydroxytryptamine 2C receptor",
          "Transient receptor potential cation channel subfamily V member 1",
          "Prostaglandin G/H synthase 2"]:
    print("CHECK", k, "->", res.get(k))
print("\nremaining unmapped (likely non-human), sample:")
um = [k for k, v in res.items() if not v[0]]
for x in um[:25]: print("  ", x)
pd.DataFrame([{"target_name": k, "gene": v[0], "match_type": v[1]} for k, v in res.items()]) \
  .to_csv(os.path.join(OUT,"target_gene_mapping.csv"), index=False, encoding="utf-8-sig")
