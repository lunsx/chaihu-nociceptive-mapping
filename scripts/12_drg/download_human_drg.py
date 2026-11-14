"""Download human DRG local.h5ad (27 MB) from the CELLxGENE datasets host."""
import os

import requests

WORKDIR = os.environ.get("CHAIHU_DRG_WORKDIR", os.path.dirname(__file__))
DSID = "6a6e7f1e-26b0-4f0a-ab34-66e9dec3bed2"
url = f"https://datasets.cellxgene.cziscience.com/{DSID}.h5ad"
s = requests.Session()
with s.get(url, stream=True, timeout=120) as resp:
    resp.raise_for_status()
    total = 0
    out = os.path.join(WORKDIR, "human_drg.h5ad")
    with open(out, "wb") as fh:
        for chunk in resp.iter_content(1 << 20):
            fh.write(chunk)
            total += len(chunk)
            print(f"\r{total/1e6:.1f} MB", end="")
print("\ndone")
