# -*- coding: utf-8 -*-
import os, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import linear_sum_assignment

work = r"C:\dockwork\redock"

# 晶体（重原子，键级正确）
cry_h = Chem.MolFromMolFile(os.path.join(work, "RCX_H.mol"), removeHs=False)
cry = Chem.RemoveHs(cry_h)

# 对接输出：读完整文件（保留头部 REMARK SMILES），只取第 1 姿态
raw = open(os.path.join(work, "RCX_redock_out.pdbqt"), encoding="utf-8").read()
from meeko import PDBQTMolecule, RDKitMolCreate
pm = PDBQTMolecule(raw, poses_to_read=1, skip_typing=True)
docked_h = RDKitMolCreate.from_pdbqt_mol(pm)[0]
docked = Chem.RemoveHs(docked_h)
print("crystal heavy:", cry.GetNumAtoms(), " docked heavy:", docked.GetNumAtoms())

# 口径1：对称校正、对齐后的最佳 RMSD
best_rms = AllChem.GetBestRMS(docked, cry)
print("GetBestRMS (symmetry, aligned): %.3f A" % best_rms)

# 口径2：同一坐标框、不做整体对齐，按元素匈牙利匹配
def heavy_coords(mol):
    out = []
    conf = mol.GetConformer()
    for a in mol.GetAtoms():
        out.append((a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx()))))
    return out

c, d = heavy_coords(cry), heavy_coords(docked)
elems = sorted(set(x[0] for x in c))
sq, n = 0.0, 0
for el in elems:
    pc = np.array([p for s, p in c if s == el])
    pd_ = np.array([p for s, p in d if s == el])
    cost = np.linalg.norm(pc[:, None, :] - pd_[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(cost)
    sq += float((cost[ri, ci] ** 2).sum())
    n += len(ri)
direct = np.sqrt(sq / n)
print("direct frame RMSD (no alignment, symmetry): %.3f A" % direct)

import csv
with open(os.path.join(work, "redock_rmsd.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["receptor", "ligand", "best_affinity_kcal_mol", "GetBestRMS_aligned_A", "direct_RMSD_A", "pass_2A"])
    w.writerow(["5KIR", "RCX(rofecoxib)", -6.4, round(best_rms, 3), round(direct, 3),
                "yes" if min(best_rms, direct) < 2.0 else "no"])
print("saved redock_rmsd.csv")
