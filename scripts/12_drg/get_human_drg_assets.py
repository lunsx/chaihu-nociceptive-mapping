"""Inspect the human DRG (Nguyen eLife 2021) H5AD asset."""
import json

import requests

DSID = "6a6e7f1e-26b0-4f0a-ab34-66e9dec3bed2"
s = requests.Session()
r = s.get("https://api.cellxgene.cziscience.com/dp/v1/datasets/index", timeout=90)
d = r.json()
rec = next(x for x in d if x["id"] == DSID)
print(json.dumps(rec.get("dataset_assets", []), indent=1)[:3000])
