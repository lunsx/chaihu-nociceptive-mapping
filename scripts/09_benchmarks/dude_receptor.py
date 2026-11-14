# -*- coding: utf-8 -*-
import http.client
import time
import urllib.request
from pathlib import Path

BASE = "https://dude.docking.org/targets/pgh2/"
OUT = Path(__file__).resolve().parent / "dude_pgh2"

def robust(url, tries=6):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except (http.client.IncompleteRead, TimeoutError, OSError) as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise last

for f in ["receptor.pdb", "crystal_ligand.mol2"]:
    dst = OUT / f
    if dst.exists():
        print("exists", f)
        continue
    dst.write_bytes(robust(BASE + f))
    print("downloaded", f, dst.stat().st_size)

print((OUT / "crystal_ligand.mol2").read_text()[:500])
