# -*- coding: utf-8 -*-
import os, requests, time
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
DST = os.path.join(OUT, "uniprot_human_reviewed.tsv")
PX = {"http": None, "https": None}
H = {"User-Agent": "Mozilla/5.0"}
URL = ("https://rest.uniprot.org/uniprotkb/stream?format=tsv"
       "&query=organism_id%3A9606+AND+reviewed%3Atrue"
       "&fields=accession%2Cprotein_name%2Cgene_names%2Cgene_primary")
def go():
    for t in range(1, 5):
        try:
            with requests.get(URL, headers=H, proxies=PX, timeout=120, stream=True) as r:
                print("status", r.status_code)
                with open(DST, "wb") as f:
                    for ch in r.iter_content(65536):
                        if ch: f.write(ch)
            if os.path.getsize(DST) > 500000:
                return os.path.getsize(DST)
        except Exception as e:
            print("retry", t, repr(e)[:100]); time.sleep(2*t)
    return -1
print("size:", go())
