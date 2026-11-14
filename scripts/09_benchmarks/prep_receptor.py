# -*- coding: utf-8 -*-
"""Prepare DUD-E pgh2 receptor: pocket center from crystal ligand and PDBQT."""
import json
import subprocess
from pathlib import Path

import numpy as np

D = Path(__file__).resolve().parent / "dude_pgh2"

# --- pocket center from crystal ligand mol2 ATOM records ---
coords = []
in_atom = False
for line in (D / "crystal_ligand.mol2").read_text().splitlines():
    if line.startswith("@<TRIPOS>ATOM"):
        in_atom = True
        continue
    if line.startswith("@<TRIPOS>"):
        in_atom = False
    if in_atom and line.strip():
        p = line.split()
        coords.append([float(p[2]), float(p[3]), float(p[4])])
coords = np.array(coords)
center = coords.mean(axis=0)
print("pocket center", np.round(center, 3), "from", len(coords), "atoms")

# --- receptor pdbqt via obabel ---
src = D / "receptor.pdb"
dst = D / "receptor.pdbqt"
if not dst.exists():
    r = subprocess.run(["obabel", "-ipdb", str(src), "-opdbqt", "-O", str(dst)],
                       capture_output=True, timeout=120)
    print("obabel rc", r.returncode, r.stdout[:120], r.stderr[:120])
print("receptor pdbqt size", dst.stat().st_size)

json.dump({"center": [float(x) for x in center], "box_size": [22.0, 22.0, 22.0]},
          open(D / "pocket.json", "w"), indent=2)
