# -*- coding: utf-8 -*-
"""Recompute redock RMSDs from existing output PDBQTs without re-docking."""
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from redock_pipeline import JOBS, PUBCHEM_NAME, extract_crystal, pose_to_mol, pubchem_smiles, prepare_ligand

WORK = Path(r"C:\dockwork\redock2")
rows = []
for pid, pdb, resn, receptor, center in JOBS:
    smiles = pubchem_smiles(PUBCHEM_NAME[resn])
    crystal_pdb, _ = extract_crystal(pdb, resn)
    crystal_mol, _ = prepare_ligand(crystal_pdb, smiles, f"{pid}_{resn}")
    docked_mol = pose_to_mol(WORK / f"{pid}_{resn}_redock_out.pdbqt", smiles)
    ref, prb = Chem.RemoveHs(crystal_mol), docked_mol
    best = float(AllChem.GetBestRMS(ref, prb, maxMatches=100))
    from rdkit.Chem import rdMolAlign
    _, _, pairs = rdMolAlign.GetBestAlignmentTransform(prb, ref, maxMatches=100)
    c1 = ref.GetConformer().GetPositions()
    c2 = prb.GetConformer().GetPositions()
    direct = float(np.sqrt(np.mean([np.sum((c1[j] - c2[i]) ** 2) for i, j in pairs])))
    print(f"{pid}/{resn}: aligned {best:.3f} A | direct {direct:.3f} A")
    # affinity from the existing CSV row / re-read output remark
    aff = None
    for line in (WORK / f"{pid}_{resn}_redock_out.pdbqt").read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            aff = float(line.split()[3])
            break
    rows.append((pid, resn, aff, best, direct))

csv = "structure,crystal_ligand,best_affinity_kcal_mol,GetBestRMS_aligned_A,direct_RMSD_A,pass_2A\n"
for pid, resn, aff, best, direct in rows:
    csv += f"{pid},{resn},{aff},{best:.3f},{direct:.3f},{'yes' if direct < 2 else 'no'}\n"
(WORK / "redock2_rmsd.csv").write_text(csv, encoding="utf-8")
print(csv)
