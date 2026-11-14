# -*- coding: utf-8 -*-
"""Dock named clinical inhibitors and luteolin; rank against the DUD-E benchmark."""
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

W = Path(__file__).resolve().parent
PC = json.load(open("C:/Users/Administrator/Desktop/草药/新选题_醋香附胁痛/计算分析/补强分析_20261114/positive_controls.json"))
LUTEOLIN = "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
items = {**{k: v["SMILES"] for k, v in PC.items()}, "luteolin": LUTEOLIN}

def build(smi, name):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    p = AllChem.ETKDGv3(); p.randomSeed = 42
    AllChem.EmbedMolecule(m, p)
    AllChem.MMFFOptimizeMolecule(m, maxIters=300)
    m.SetProp("_Name", name)
    return m

sdf = W / "controls.sdf"
with open(sdf, "w", encoding="utf-8", newline="") as fh:
    w = Chem.SDWriter(fh)
    for name, smi in items.items():
        w.write(build(smi, name))
    w.close()
subprocess.run(["mk_prepare_ligand", "-i", str(sdf), "--multimol_outdir",
                str(W / "control_pdbqt"), "--multimol_prefix", "c"], check=True)

bench = pd.read_csv(W / "docking_results.csv")
all_scores = bench["affinity"].dropna().values
res = {}
for f in sorted((W / "control_pdbqt").glob("*.pdbqt")):
    out = W / "control_docked" / f.name
    out.parent.mkdir(exist_ok=True)
    argv = ["C:/dockwork/vina.exe", "--receptor", str(W / "receptor_rigid.pdbqt"),
            "--ligand", str(f), "--center_x", "30.801", "--center_y", "-22.631",
            "--center_z", "-16.008", "--size_x", "22", "--size_y", "22", "--size_z", "22",
            "--exhaustiveness", "8", "--num_modes", "8", "--seed", "42",
            "--out", str(out), "--verbosity", "1"]
    r = subprocess.run(argv, capture_output=True, timeout=120)
    txt = r.stdout.decode(errors="replace")
    aff = float(re.search(r"^\s*1\s+(-?\d+\.\d+)", txt, re.M).group(1))
    # percentile: fraction of benchmark ligands with affinity >= this (more negative)
    pct = 100 * (all_scores <= aff).mean()
    res[f.stem] = {"affinity": aff, "benchmark_percentile": round(float(pct), 1)}

# map meeko index back to names
order = list(items.keys())
named = {}
for i, name in enumerate(order, start=1):
    key = f"c-{i}"
    if key in res:
        named[name] = res[key]
json.dump(named, open(W / "controls_luteolin.json", "w"), indent=2)
print(json.dumps(named, indent=2))
