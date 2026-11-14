# -*- coding: utf-8 -*-
"""Dock all benchmark ligands and record best affinity. Resumable."""
import json
import re
import subprocess
from pathlib import Path

from rdkit import Chem

W = Path(__file__).resolve().parent
PDBQT = W / "pdbqt"
OUT = W / "docked"
LOGS = W / "logs"
OUT.mkdir(exist_ok=True)
LOGS.mkdir(exist_ok=True)

# index -> class/id from SDF order (matches meeko lig-N numbering)
suppl = Chem.SDMolSupplier(str(W / "benchmark_ligands.sdf"), removeHs=False)
idmap = {}
for i, mol in enumerate(suppl, start=1):
    idmap[i] = {"name": mol.GetProp("_Name"),
                "cls": mol.GetProp("lig_class"), "cid": mol.GetProp("lig_id")}
json.dump(idmap, open(W / "idmap.json", "w"), indent=2)

center = json.load(open(Path("C:/Users/Administrator/Desktop/草药/新选题_醋香附胁痛/计算分析/补强分析_20261114/dude_pgh2/pocket.json")))
cx, cy, cz = center["center"]
results_csv = W / "docking_results.csv"
done = {}
if results_csv.exists():
    for line in results_csv.read_text(encoding="utf-8").splitlines()[1:]:
        p = line.split(",")
        done[int(p[0])] = True

def run_one(i):
    lig = PDBQT / f"lig-{i}.pdbqt"
    out = OUT / f"lig-{i}.pdbqt"
    argv = ["C:/dockwork/vina.exe", "--receptor", str(W / "receptor_rigid.pdbqt"),
            "--ligand", str(lig), "--center_x", str(cx), "--center_y", str(cy),
            "--center_z", str(cz), "--size_x", "22", "--size_y", "22", "--size_z", "22",
            "--exhaustiveness", "8", "--num_modes", "8", "--seed", "42",
            "--out", str(out), "--verbosity", "1"]
    r = subprocess.run(argv, capture_output=True, timeout=180)
    txt = r.stdout.decode(errors="replace")
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", txt, re.M)
    return float(m.group(1)) if m else None

with open(results_csv, "a", encoding="utf-8") as fh:
    if not results_csv.exists() or results_csv.stat().st_size == 0:
        fh.write("index,cls,cid,affinity\n")
    for i in range(1, len(idmap) + 1):
        if i in done:
            continue
        aff = run_one(i)
        info = idmap[i]
        fh.write(f"{i},{info['cls']},{info['cid']},{aff}\n")
        fh.flush()
        if i % 20 == 0:
            print("docked", i, "/", len(idmap))
print("ALL DONE", len(idmap))
