# -*- coding: utf-8 -*-
import os, time, csv, requests
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
U = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951&targ=self&form=text&view=full"
RAW = os.path.join(OUT, "GPL14951_platform_table.txt")
MAP = os.path.join(OUT, "probe_gene_map.csv")
HDR = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
NO_PROXY = {"http": None, "https": None}

def stream_download(tries=5):
    for t in range(1, tries+1):
        try:
            with requests.get(U, headers=HDR, timeout=60, proxies=NO_PROXY, stream=True) as r:
                with open(RAW, "wb") as f:
                    for chunk in r.iter_content(65536):
                        if chunk: f.write(chunk)
            if os.path.getsize(RAW) > 1_000_000:
                return os.path.getsize(RAW)
        except Exception as e:
            print("retry", t, e); time.sleep(2*t)
    return -1

sz = stream_download()
print("platform raw bytes:", sz)
if sz < 0:
    raise SystemExit("download failed")

# parse table
in_table = False
header = None
n_rows = 0
with open(RAW, "r", encoding="utf-8", errors="replace") as f, \
     open(MAP, "w", newline="", encoding="utf-8") as out:
    w = csv.writer(out)
    w.writerow(["probe_id", "symbol", "entrez_id", "ilmn_gene"])
    for line in f:
        if line.startswith("!platform_table_begin"):
            in_table = True; continue
        if line.startswith("!platform_table_end"):
            in_table = False; break
        if not in_table: continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = {name: i for i, name in enumerate(parts)}
            continue
        def g(col):
            i = header.get(col)
            return parts[i] if i is not None and i < len(parts) else ""
        probe = g("ID"); sym = g("Symbol"); ent = g("Entrez_Gene_ID"); ilmn = g("ILMN_Gene")
        w.writerow([probe, sym, ent, ilmn])
        n_rows += 1
print("mapped probes:", n_rows)
