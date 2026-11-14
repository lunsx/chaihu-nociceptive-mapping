# -*- coding: utf-8 -*-
"""Build 3D benchmark ligands from DUD-E SMILES and convert to PDBQT with Meeko."""
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

D = Path(__file__).resolve().parent / "dude_pgh2"
PDBQT_DIR = D / "ligands_pdbqt"
PDBQT_DIR.mkdir(exist_ok=True)

N_ACTIVE, N_DECOY, SEED = 50, 250, 42

def parse_ism(f):
    out = []
    for line in (D / f).read_text().splitlines():
        p = line.split()
        if len(p) >= 2:
            out.append((p[0], p[1]))
    return out

def select(rows, n):
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(rows), size=min(n, len(rows)), replace=False)
    return [rows[i] for i in sorted(idx)]

def build(smiles, cid, label):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    except Exception:
        pass
    mol.SetProp("_Name", f"{label}_{cid}")
    mol.SetProp("lig_id", cid)
    mol.SetProp("lig_class", label)
    return mol

act = select(parse_ism("actives_final.ism"), N_ACTIVE)
dec = select(parse_ism("decoys_final.ism"), N_DECOY)

mols, failed = [], []
for label, rows in (("active", act), ("decoy", dec)):
    for smi, cid in rows:
        m = build(smi, cid, label)
        (mols if m is not None else failed).append(m if m is not None else (label, cid))

sdf = D / "benchmark_ligands.sdf"
print("built", len(mols), "ligands; failed", len([f for f in failed if f]))
with open(sdf, "w", encoding="utf-8", newline="") as fh:
    w = Chem.SDWriter(fh)
    for m in mols:
        w.write(m)
    w.close()
print("sdf written", sdf.stat().st_size)

# Meeko: one pdbqt per molecule
r = subprocess.run(["mk_prepare_ligand", "-i", str(sdf), "--multimol_outdir", str(PDBQT_DIR),
                    "--multimol_prefix", "lig"], capture_output=True, timeout=300)
made = list(PDBQT_DIR.glob("*.pdbqt"))
print("meeko pdbqt files:", len(made))
print(r.stderr[:200])
