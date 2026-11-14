# -*- coding: utf-8 -*-
import os, requests, numpy as np
PX = {"http": None, "https": None}
work = r"C:\dockwork\ensemble"
pdb_path = os.path.join(work, "6COX.pdb")

# S58 身份
r = requests.get("https://data.rcsb.org/rest/v1/core/chemcomp/S58", proxies=PX, timeout=40)
if r.status_code == 200:
    print("S58 name:", r.json().get("chem_comp", {}).get("name"),
          "formula:", r.json().get("chem_comp", {}).get("formula"))

# 计算 S58 chain A 中心，并清洗受体（仅留链 A 蛋白 ATOM）
coords = []
clean = []
with open(pdb_path, encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("HETATM") and line[17:20].strip() == "S58" and line[21].strip() == "A":
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if line.startswith("ATOM") and line[21].strip() == "A":
            clean.append(line)
center = np.array(coords).mean(axis=0)
print("S58 centroid:", [round(x, 3) for x in center], "from", len(coords), "atoms")

clean_path = os.path.join(work, "6COX_clean.pdb")
with open(clean_path, "w", encoding="utf-8", newline="\n") as f:
    f.write("".join(clean) + "END\n")
print("cleaned receptor ATOM lines:", len(clean), "->", clean_path)

import json
with open(os.path.join(work, "6COX_center.json"), "w") as f:
    json.dump([round(x, 3) for x in center], f)
