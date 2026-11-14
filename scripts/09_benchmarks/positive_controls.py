# -*- coding: utf-8 -*-
"""Fetch canonical isomeric SMILES for named PTGS2 clinical inhibitors (positive controls)."""
import json
from pathlib import Path
import urllib.request
import urllib.parse

def pubchem(name):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name) + "/property/IsomericSMILES,CanonicalSMILES,MolecularFormula/JSON")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    props = json.loads(urllib.request.urlopen(req, timeout=40).read())["PropertyTable"]["Properties"][0]
    return props

controls = {}
for name in ["celecoxib", "diclofenac", "ibuprofen", "indomethacin"]:
    controls[name] = pubchem(name)
    print(name, controls[name])

json.dump(controls, open(Path(__file__).resolve().parent / "positive_controls.json", "w"),
          ensure_ascii=False, indent=2)
