# -*- coding: utf-8 -*-
"""Resumable, retrying downloader for the GSE344087 RAW tar."""
import time
import urllib.request
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE344nnn/GSE344087/suppl/GSE344087_RAW.tar"
DST = Path(__file__).resolve().parent / "GSE344087_RAW.tar"
TOTAL = 170096640

def have():
    return DST.stat().st_size if DST.exists() else 0

for attempt in range(60):
    pos = have()
    if pos >= TOTAL:
        break
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
    if pos:
        req.add_header("Range", f"bytes={pos}-")
    try:
        with urllib.request.urlopen(req, timeout=60) as r, open(DST, "ab") as fh:
            while True:
                chunk = r.read(1 << 16)
                if not chunk:
                    break
                fh.write(chunk)
    except Exception as e:
        print("attempt", attempt, "at", have(), "retry:", type(e).__name__)
        time.sleep(2)

print("final size", have(), "complete:", have() == TOTAL)
