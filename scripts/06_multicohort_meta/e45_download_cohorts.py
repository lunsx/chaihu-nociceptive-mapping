# -*- coding: utf-8 -*-
"""下载 5 个独立队列数据到 multicohort 目录并校验。"""
import os, gzip, hashlib
import requests

PROXIES = {"http": "http://127.0.0.1:8902", "https": "http://127.0.0.1:8902"}
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
OUT = "multicohort"

FILES = {
    "GSE281797": [
        "GSE281nnn/GSE281797/suppl/GSE281797_GeneCount.txt.gz",
        "GSE281nnn/GSE281797/suppl/GSE281797_tpm.txt.gz",
    ],
    "GSE260666": [
        "GSE260nnn/GSE260666/suppl/GSE260666_raw_counts.txt.gz",
    ],
    "GSE147304": [
        "GSE147nnn/GSE147304/suppl/GSE147304_gene_expression.txt.gz",
    ],
    "GSE193066": [
        "GSE193nnn/GSE193066/suppl/GSE193066_NAFLD.HUn164.gct.gz",
    ],
    "GSE192959": [
        "GSE192nnn/GSE192959/suppl/GSE192959_NAFLD.TUn48.gct.gz",
    ],
}

os.makedirs(OUT, exist_ok=True)
for gse, paths in FILES.items():
    d = os.path.join(OUT, gse)
    os.makedirs(d, exist_ok=True)
    for p in paths:
        url = f"{BASE}/{p}"
        name = p.split("/")[-1]
        dst = os.path.join(d, name)
        r = requests.get(url, proxies=PROXIES, timeout=300)
        r.raise_for_status()
        with open(dst, "wb") as f:
            f.write(r.content)
        # 验证可 gunzip
        try:
            with gzip.open(dst, "rb") as g:
                head = g.read(200)
            ok = "gunzip OK"
        except Exception as e:
            ok = f"GUNZIP FAIL {e}"
        print(f"{gse}/{name}  bytes={len(r.content)}  {ok}")
print("\nDownloaded tree:")
for root, _, fs in os.walk(OUT):
    for f in fs:
        fp = os.path.join(root, f)
        print(" ", fp, os.path.getsize(fp))
