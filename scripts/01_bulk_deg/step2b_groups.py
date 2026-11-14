# -*- coding: utf-8 -*-
import gzip
P = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\GSE89632_series_matrix.txt.gz"
with gzip.open(P, "rt", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()
def fields(prefix):
    out=[]
    for ln in lines:
        if ln.startswith(prefix):
            out.append([x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
    return out
titles = fields("!Sample_title")[0]
gsm = fields("!Sample_geo_accession")[0]
chars = fields("!Sample_characteristics_ch1")
print("n samples:", len(titles))
print("\n=== first sample title ===", titles[0])
print("\n=== characteristic labels (from sample 0) ===")
for c in chars:
    print("  ", c[0])
print("\n=== group-like characteristic across all samples ===")
for c in chars:
    joined = " | ".join(c)
    if any(k in joined.lower() for k in ["group", "diagnosis", "condition", "nash", "status", "class"]):
        from collections import Counter
        print(Counter(c))
