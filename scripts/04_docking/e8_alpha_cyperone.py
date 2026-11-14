# -*- coding: utf-8 -*-
import os, time, requests, json
PX = {"http": None, "https": None}
work = r"C:\dockwork\alpha"
os.makedirs(work, exist_ok=True)
CID = 6452086

def retry_get(url, tries=8, want="text"):
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, proxies=PX, timeout=45)
            if r.status_code == 200:
                return r
            last = "HTTP %s" % r.status_code
        except Exception as e:
            last = str(e)
        time.sleep(1.2 * (i + 1))
    raise RuntimeError("fail after %d tries: %s url=%s" % (tries, last, url))

# 属 性
r = retry_get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/MolecularFormula,ConnectivitySMILES/JSON")
props = r.json()["PropertyTable"]["Properties"][0]
print("formula:", props.get("MolecularFormula"))
print("smiles:", props.get("ConnectivitySMILES"))

# 3D SDF
r = retry_get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/SDF?record_type=3d")
sdf = r.text
sdf_path = os.path.join(work, "alpha_cyperone.sdf")
with open(sdf_path, "w", encoding="utf-8", newline="\n") as f:
    f.write(sdf)
print("3D SDF chars:", len(sdf), "->", sdf_path)
with open(os.path.join(work, "alpha_info.txt"), "w", encoding="utf-8") as f:
    f.write(json.dumps(props, ensure_ascii=False))
