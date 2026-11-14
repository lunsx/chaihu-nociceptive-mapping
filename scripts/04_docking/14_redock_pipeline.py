# -*- coding: utf-8 -*-
"""Redock co-crystal ligands into their own grids and compute heavy-atom RMSD.

Structures:
  1CX2  SC-558 (S58)  selective COX-2 inhibitor
  6COX  SC-558 (S58)  same ligand, I222 crystal form
  3PGH  flurbiprofen (FLP) non-selective inhibitor
Same rigid receptors, grid centers, box (22 A) as run_ensemble.py.
"""
import json
import os
import pathlib
import subprocess
import time

import numpy as np
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

WORK = pathlib.Path(r"C:\dockwork\redock2")
WORK.mkdir(exist_ok=True)
MK = r"C:\Users\Administrator\AppData\Local\Doubao\User Data\sandbox_runtime\bases\772984f6fd3a8ae0ae7d2422cd7c2bfa\python\Scripts\mk_prepare_ligand.exe"
VINA = r"C:\dockwork\vina.exe"

JOBS = [
    ("1CX2", r"C:\dockwork\ens2\1CX2.pdb", "S58",
     r"C:\dockwork\ens2\1CX2_rigid.pdbqt", (24.26, 21.53, 16.50)),
    ("6COX", r"C:\dockwork\ensemble\6COX.pdb", "S58",
     r"C:\dockwork\ensemble\6COX.pdbqt", None),
    ("3PGH", r"C:\dockwork\ens2\3PGH.pdb", "FLP",
     r"C:\dockwork\ens2\3PGH_rigid.pdbqt", (25.71, 22.23, 15.05)),
]
PUBCHEM_NAME = {"S58": "SC-558", "FLP": "flurbiprofen"}


def pubchem_smiles(name):
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    candidates = [name]
    if name == "SC-558":
        candidates = ["SC-558", "SC558"]
    data = None
    for q in candidates:
        r = requests.get(f"{base}/compound/name/{q}/property/IsomericSMILES/JSON", timeout=40)
        if r.status_code == 200:
            data = r.json()
            break
    if data is None:
        raise RuntimeError(f"PubChem name lookup failed for {name}")
    props = data["PropertyTable"]["Properties"]
    return props[0].get("IsomericSMILES") or props[0]["SMILES"]


def extract_crystal(pdb, resn):
    lines = []
    with open(pdb, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if (line.startswith("HETATM") and line[17:20].strip() == resn
                    and line[21].strip() == "A" and line[22:26].strip() == "701"):
                lines.append(line)
    out = WORK / f"{resn}_{pathlib.Path(pdb).stem}_crystal.pdb"
    out.write_text("".join(lines) + "END\n", encoding="utf-8")
    return out, len(lines)


def prepare_ligand(crystal_pdb, smiles, tag):
    pdbmol = Chem.MolFromPDBBlock(crystal_pdb.read_text(), removeHs=False, sanitize=False)
    template = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, pdbmol)
    Chem.SanitizeMol(mol)
    mol_h = Chem.AddHs(mol, addCoords=True)
    mol_file = WORK / f"{tag}.mol"
    mol_file.write_text(Chem.MolToMolBlock(mol_h), encoding="utf-8")
    pdbqt = WORK / f"{tag}.pdbqt"
    p = subprocess.run([MK, "-i", str(mol_file), "-o", str(pdbqt)],
                       capture_output=True, text=True, timeout=120)
    if p.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(f"mk_prepare failed: {p.stderr}")
    return mol_h, pdbqt


def dock(receptor, center, lig_pdbqt, tag):
    out_pdbqt = WORK / f"{tag}_redock_out.pdbqt"
    cmd = [VINA, "--receptor", receptor, "--ligand", str(lig_pdbqt),
           "--center_x", str(center[0]), "--center_y", str(center[1]),
           "--center_z", str(center[2]),
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", "16", "--num_modes", "5",
           "--seed", "42", "--cpu", "4",
           "--out", str(out_pdbqt)]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=900)
    best = None
    for line in p.stdout.splitlines():
        if line.strip().startswith("1 ") and best is None:
            parts = line.split()
            best = float(parts[1])
    if best is None:
        raise RuntimeError(f"vina failed: {p.stderr[:200]} {p.stdout[:200]}")
    return best, out_pdbqt, round(time.time() - t0, 1)


def pose_to_mol(out_pdbqt, smiles):
    """Take MODEL 1 ATOM lines and rebuild a molecule with template bond orders.

    Vina PDBQT stores AutoDock atom types (A, OA, ...) in the element columns,
    which RDKit cannot parse. We rebuild clean PDB lines with the element
    inferred from the atom name, keep heavy atoms only, then assign bond orders
    from the PubChem SMILES template.
    """
    keep = []
    in_model1 = True
    for line in out_pdbqt.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL") and "1" not in line.split():
            in_model1 = False
        if in_model1 and line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            elem = "".join(ch for ch in name if ch.isalpha())
            if elem not in ("Br", "Cl"):
                elem = elem[0]
            if elem == "H":
                continue
            clean = line[:66] + " " * 10 + elem.rjust(2)
            keep.append(clean)
    pdb_like = "\n".join(keep) + "\nEND\n"
    pdbmol = Chem.MolFromPDBBlock(pdb_like, removeHs=True, sanitize=False)
    if pdbmol is None:
        raise RuntimeError("failed to parse cleaned redock PDB")
    template = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, pdbmol)
    Chem.SanitizeMol(mol)
    return mol


def rmsds(crystal_mol, docked_mol):
    ref = Chem.RemoveHs(crystal_mol)
    prb = docked_mol
    assert ref.GetNumAtoms() == prb.GetNumAtoms()
    # free alignment: atom correspondence is searched, not assumed by index
    best_rmsd = float(AllChem.GetBestRMS(ref, prb, maxMatches=100))
    from rdkit.Chem import rdMolAlign
    _, _, pairs = rdMolAlign.GetBestAlignmentTransform(prb, ref, maxMatches=100)
    c1 = ref.GetConformer().GetPositions()
    c2 = prb.GetConformer().GetPositions()
    d2 = [np.sum((c1[j] - c2[i]) ** 2) for i, j in pairs]
    direct = float(np.sqrt(np.mean(d2)))
    return best_rmsd, direct


def main():
    center_6cox = json.loads(pathlib.Path(r"C:\dockwork\ensemble\6COX_center.json").read_text())
    rows = []
    smiles_cache = {}
    for pid, pdb, resn, receptor, center in JOBS:
        if pid == "6COX":
            center = tuple(center_6cox)
        print(f"==== {pid} / {resn} ====", flush=True)
        if resn not in smiles_cache:
            smiles_cache[resn] = pubchem_smiles(PUBCHEM_NAME[resn])
        smiles = smiles_cache[resn]
        print("smiles:", smiles, flush=True)
        crystal_pdb, n = extract_crystal(pdb, resn)
        print("crystal atoms:", n, flush=True)
        crystal_mol, lig_pdbqt = prepare_ligand(crystal_pdb, smiles, f"{pid}_{resn}")
        print("prepared pdbqt:", lig_pdbqt.name, flush=True)
        best, out_pdbqt, secs = dock(receptor, center, lig_pdbqt, f"{pid}_{resn}")
        print("docked best:", best, f"({secs}s)", flush=True)
        docked_mol = pose_to_mol(out_pdbqt, smiles)
        best_rmsd, direct = rmsds(crystal_mol, docked_mol)
        print(f"bestRMS {best_rmsd:.3f} A | direct RMSD {direct:.3f} A", flush=True)
        rows.append((pid, resn, best, best_rmsd, direct, secs))
    csv = "structure,crystal_ligand,best_affinity_kcal_mol,GetBestRMS_aligned_A,direct_RMSD_A,seconds,pass_2A\n"
    for pid, resn, best, br, dr, secs in rows:
        csv += f"{pid},{resn},{best},{br:.3f},{dr:.3f},{secs},{'yes' if dr < 2 else 'no'}\n"
    (WORK / "redock2_rmsd.csv").write_text(csv, encoding="utf-8")
    print(csv)


if __name__ == "__main__":
    main()
