# -*- coding: utf-8 -*-
import os, pickle
from openmm import *
from openmm.app import *
from openmm.unit import nanometer, molar, amu
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem

WORK = r"C:\dockwork\md"
PROT_PDB = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking\receptors\5KIR_clean.pdb"
LIG_SDF = os.path.join(WORK, "luteolin_docked.sdf")
SAGE = os.path.join(WORK, "openff-2.1.0.offxml")

# ---------- 1. 蛋白：PDBFixer 补缺失原子，再去现有氢 ----------
from pdbfixer import PDBFixer
fixer = PDBFixer(filename=PROT_PDB)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
mod = Modeller(fixer.topology, fixer.positions)
mod.deleteWater()
h_atoms = [a for a in mod.topology.atoms() if a.element == Element.getBySymbol("H")]
mod.delete(h_atoms)
print("protein heavy atoms:", mod.topology.getNumAtoms())

# ---------- 2. 配体（带对接坐标，MMFF94 电荷） ----------
offmol = Molecule.from_rdkit(Chem.MolFromMolFile(LIG_SDF, removeHs=False))
offmol.assign_partial_charges(partial_charge_method="mmff94")
lig_top = offmol.to_topology().to_openmm()
lig_pos = offmol.conformers[0].to_openmm()
print("ligand atoms:", lig_top.getNumAtoms())

# ---------- 3. 力场：ff14SB + TIP3P，注册 Sage(SMIRNOFF) ----------
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
smirnoff = SMIRNOFFTemplateGenerator(molecules=[offmol], forcefield=SAGE)
ff.registerTemplateGenerator(smirnoff.generator)

# 蛋白补氢（ff14SB）
mod.addHydrogens(ff)
# 记录蛋白原子数（补氢后），用于定位配体索引
n_prot = mod.topology.getNumAtoms()
# 加入配体（对接姿态）
mod.add(lig_top, lig_pos)
n_lig = lig_top.getNumAtoms()
lig_idx = list(range(n_prot, n_prot + n_lig))
# 水盒 + 0.15M 离子
mod.addSolvent(ff, model="tip3p", padding=1.1 * nanometer,
               ionicStrength=0.15 * molar, neutralize=True)
print("total atoms after solvation:", mod.topology.getNumAtoms())
print("ligand atom indices:", lig_idx[0], "..", lig_idx[-1])

# ---------- 4. 体系：PME + HMR(4fs) ----------
system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                        nonbondedCutoff=1.0 * nanometer,
                        constraints=HBonds, rigidWater=True,
                        hydrogenMass=3.0 * amu)

with open(os.path.join(WORK, "system.xml"), "w") as f:
    f.write(XmlSerializer.serialize(system))
PDBFile.writeFile(mod.topology, mod.positions, open(os.path.join(WORK, "complex_solvated.pdb"), "w"))
with open(os.path.join(WORK, "modeller.pkl"), "wb") as f:
    pickle.dump({"topology": mod.topology, "positions": mod.positions,
                 "lig_idx": lig_idx, "n_prot": n_prot}, f)
print("saved system.xml, complex_solvated.pdb, modeller.pkl")
