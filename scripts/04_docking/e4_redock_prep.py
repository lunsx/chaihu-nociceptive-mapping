# -*- coding: utf-8 -*-
import os, requests
PX = {"http": None, "https": None}
raw = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking\receptors\5KIR.pdb"
work = r"C:\dockwork\redock"
os.makedirs(work, exist_ok=True)

# 提取 RCX chain A resseq 601
lines = []
with open(raw, encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("HETATM") and line[17:20].strip() == "RCX" and line[21].strip() == "A" and line[22:26].strip() == "601":
            lines.append(line)
out_pdb = os.path.join(work, "RCX_crystal.pdb")
with open(out_pdb, "w", encoding="utf-8", newline="\n") as f:
    f.write("".join(lines) + "END\n")
print("crystal RCX atoms:", len(lines), "->", out_pdb)

# PubChem 取 rofecoxib
def get(url):
    for i in range(4):
        try:
            r = requests.get(url, timeout=40, proxies=PX)
            if r.status_code == 200:
                return r
        except Exception:
            pass
    return None

j = get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/rofecoxib/cids/JSON")
cid = j.json()["IdentifierList"]["CID"][0]
print("rofecoxib CID:", cid)
p = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES/JSON")
props = p.json()["PropertyTable"]["Properties"][0]
print("SMILES:", props["CanonicalSMILES"])
with open(os.path.join(work, "rcx_template.txt"), "w", encoding="utf-8") as f:
    f.write(str(cid) + "\n" + props["CanonicalSMILES"] + "\n")
