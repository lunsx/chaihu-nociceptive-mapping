# -*- coding: utf-8 -*-
"""Build one solvated OpenMM system per distinct pose (multi-pose validation).

Every arm uses the SAME OpenFF molecule (fixed graph and atom order, taken from
the canonical SDF). Only the 3D coordinates differ. For blind poses the pose
PDB is mapped onto the heavy-atom template and its coordinates are remapped
into the canonical atom order through an explicit isomorphism. This keeps the
SMIRNOFF template generation identical across arms.
"""
import os, pickle, sys
from pathlib import Path

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import nanometer, molar, amu
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

PROT_PDB = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking\receptors\5KIR_clean.pdb"
CANON_SDF = r"C:\dockwork\md\luteolin_docked.sdf"
SAGE = r"C:\dockwork\md\openff-2.1.0.offxml"
WP10 = PATHS["dock_root"] / "md" / "wp10_multipose"
N_BLIND = 9

CANON_RD = Chem.MolFromMolFile(CANON_SDF, removeHs=False)
TEMPLATE = Chem.MolFromMolFile(CANON_SDF, removeHs=True)


def canonical_coords():
    c = CANON_RD.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(CANON_RD.GetNumAtoms())])


def blind_coords(pdb):
    raw = Chem.MolFromPDBFile(pdb, removeHs=True, sanitize=False)
    raw = Chem.RemoveHs(raw, sanitize=False)
    mapped = Chem.AddHs(AllChem.AssignBondOrdersFromTemplate(TEMPLATE, raw), addCoords=True)
    # explicit isomorphism: match[i] = mapped atom index for canonical atom i
    match = mapped.GetSubstructMatch(CANON_RD)
    assert len(match) == CANON_RD.GetNumAtoms(), "incomplete atom mapping"
    conf = mapped.GetConformer()
    return np.array([list(conf.GetAtomPosition(match[i])) for i in range(CANON_RD.GetNumAtoms())])


def prepare_protein():
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=PROT_PDB)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    base = Modeller(fixer.topology, fixer.positions)
    base.deleteWater()
    base.delete([a for a in base.topology.atoms()
                 if a.element == Element.getBySymbol("H")])
    return base


def build_arm(arm_id, coords):
    out_dir = WP10 / f"arm{arm_id}"
    out_dir.mkdir(parents=True, exist_ok=True)
    mod_pkl = out_dir / "modeller.pkl"
    sys_xml = out_dir / "system.xml"
    if mod_pkl.exists() and sys_xml.exists():
        print(f"arm{arm_id}: existing, skip")
        return

    offmol = Molecule.from_rdkit(CANON_RD)
    offmol.assign_partial_charges(partial_charge_method="mmff94")
    # overwrite conformer with this arm's coordinates
    offmol.conformers[0] = coords * offmol.conformers[0].units
    lig_top = offmol.to_topology().to_openmm()
    lig_pos = offmol.conformers[0].to_openmm()

    base = prepare_protein()
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[offmol], forcefield=SAGE)
    ff.registerTemplateGenerator(smirnoff.generator)
    base.addHydrogens(ff)
    n_prot = base.topology.getNumAtoms()
    base.add(lig_top, lig_pos)
    n_lig = lig_top.getNumAtoms()
    lig_idx = list(range(n_prot, n_prot + n_lig))
    base.addSolvent(ff, model="tip3p", padding=1.1 * nanometer,
                    ionicStrength=0.15 * molar, neutralize=True)
    system = ff.createSystem(base.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0 * nanometer,
                             constraints=HBonds, rigidWater=True,
                             hydrogenMass=3.0 * amu)
    with open(sys_xml, "w") as fh:
        fh.write(XmlSerializer.serialize(system))
    with open(mod_pkl, "wb") as fh:
        pickle.dump({"topology": base.topology, "positions": base.positions,
                     "lig_idx": lig_idx, "n_prot": n_prot}, fh)
    print(f"arm{arm_id}: atoms={base.topology.getNumAtoms()} saved")


def main():
    poses = [(0, canonical_coords())]
    for i in range(1, N_BLIND + 1):
        poses.append((i, blind_coords(str(WP10 / f"pose_{i}.pdb"))))
    for arm_id, coords in poses:
        print(f"== building arm{arm_id} ==")
        build_arm(arm_id, coords)
    print("ALL ARM SYSTEMS BUILT")


if __name__ == "__main__":
    main()
