"""Narrow CELLxGENE search: only dorsal root / trigeminal / peripheral ganglion."""
import re

import requests

s = requests.Session()
r = s.get("https://api.cellxgene.cziscience.com/dp/v1/datasets/index", timeout=90)
d = r.json()

pat = re.compile(r"dorsal root|trigeminal|peripheral ganglion|sensory ganglion", re.I)
hits = []
for x in d:
    labels = [str(t.get("label", "")) for t in (x.get("tissue") or [])]
    if any(pat.search(lab) for lab in labels):
        hits.append((x, labels))

print("dorsal-root/trigeminal datasets:", len(hits))
for x, labels in hits:
    print("-" * 70)
    print("id:", x["id"])
    print("name:", x.get("name"))
    print("tissue labels:", [lab for lab in labels if pat.search(lab)])
    print("organism:", [o.get("label") for o in x.get("organism", [])],
          "cells:", x.get("cell_count"))
    for a in x.get("dataset_assets", []):
        print("  asset:", a.get("files"), a.get("filetype"), a.get("url", "")[:80])
