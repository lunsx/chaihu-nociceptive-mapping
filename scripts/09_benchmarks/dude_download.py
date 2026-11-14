# -*- coding: utf-8 -*-
"""Download DUD-E pgh2 (PTGS2) actives/decoys and benchmark receptor."""
import gzip
import http.client
import time
import urllib.request
from pathlib import Path

BASE = "https://dude.docking.org/targets/pgh2/"
OUT = Path(__file__).resolve().parent / "dude_pgh2"
OUT.mkdir(exist_ok=True)

FILES = [
    "actives_final.ism", "decoys_final.ism",
    "actives_final.sdf.gz", "decoys_final.sdf.gz",
    "receptor.pdb", "crystal_ligand.mol2",
    "actives_clustered.ism",
]

def get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except (http.client.IncompleteRead, TimeoutError, OSError) as e:
            last = e
            time.sleep(2 * (i + 1))
    raise last

for f in FILES:
    dst = OUT / f
    if dst.exists():
        print("exists", f)
        continue
    data = get(BASE + f)
    dst.write_bytes(data)
    print("downloaded", f, len(data))

# count SMILES lines
for f in ["actives_final.ism", "decoys_final.ism", "actives_clustered.ism"]:
    lines = (OUT / f).read_text().splitlines()
    print(f, "n =", len(lines), "| first:", lines[0][:90] if lines else "")
