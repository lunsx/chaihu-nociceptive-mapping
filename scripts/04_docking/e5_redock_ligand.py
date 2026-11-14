# -*- coding: utf-8 -*-
import os
from rdkit import Chem
from rdkit.Chem import AllChem

work = r"C:\dockwork\redock"
smiles = "CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3"

# 晶体配体（PDB 块，rdkit 按距离推键，默认单键）
pdb_block = open(os.path.join(work, "RCX_crystal.pdb"), encoding="utf-8").read()
pdbmol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
print("pdb mol atoms:", pdbmol.GetNumAtoms())
print("pdb elements:", sorted(a.GetSymbol() for a in pdbmol.GetAtoms()))

template = Chem.MolFromSmiles(smiles)
print("template heavy atoms:", template.GetNumAtoms())
print("template elements:", sorted(a.GetSymbol() for a in template.GetAtoms()))

# 用模板分配键级（保留晶体坐标）
mol = AllChem.AssignBondOrdersFromTemplate(template, pdbmol)
Chem.SanitizeMol(mol)
print("bond-order assigned, atoms:", mol.GetNumAtoms())

# 加氢：基于晶体几何放置氢坐标，保留重原子晶体位置不动
mol_h = Chem.AddHs(mol, addCoords=True)

mol_block = Chem.MolToMolBlock(mol_h)
with open(os.path.join(work, "RCX_H.mol"), "w", encoding="utf-8", newline="\n") as f:
    f.write(mol_block)
print("written RCX_H.mol, total atoms:", mol_h.GetNumAtoms())
