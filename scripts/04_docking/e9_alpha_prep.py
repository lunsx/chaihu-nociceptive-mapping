# -*- coding: utf-8 -*-
import os
from rdkit import Chem
from meeko import MoleculePreparation

work = r"C:\dockwork\alpha"
raw = open(os.path.join(work, "alpha_cyperone.sdf"), encoding="utf-8").read()
block = raw[:raw.find("M  END") + 6]
mol = Chem.MolFromMolBlock(block, removeHs=False)
print("parsed heavy atoms:", mol.GetNumAtoms())
mol_h = Chem.AddHs(mol, addCoords=True)
print("with H:", mol_h.GetNumAtoms())

with open(os.path.join(work, "alpha_cyperone_H.mol"), "w", encoding="utf-8", newline="\n") as f:
    f.write(Chem.MolToMolBlock(mol_h))

prep = MoleculePreparation()
prep.prepare(mol_h)
pdbqt = prep.write_pdbqt_string()
with open(os.path.join(work, "alpha_cyperone.pdbqt"), "w", encoding="utf-8", newline="\n") as f:
    f.write(pdbqt)
print("pdbqt chars:", len(pdbqt))
print("\n".join(pdbqt.splitlines()[:4]))
