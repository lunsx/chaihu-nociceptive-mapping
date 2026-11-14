# -*- coding: utf-8 -*-
import os
from rdkit import Chem
from meeko import MoleculePreparation

work = r"C:\dockwork\redock"
mol = Chem.MolFromMolFile(os.path.join(work, "RCX_H.mol"), removeHs=False)
print("loaded mol atoms:", mol.GetNumAtoms())

prep = MoleculePreparation()
setups = prep.prepare(mol)
print("n setups:", len(setups))
pdbqt = prep.write_pdbqt_string()
with open(os.path.join(work, "RCX.pdbqt"), "w", encoding="utf-8", newline="\n") as f:
    f.write(pdbqt)
print("RCX.pdbqt written, chars:", len(pdbqt))
print("\n".join(pdbqt.splitlines()[:6]))
