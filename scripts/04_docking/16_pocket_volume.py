# -*- coding: utf-8 -*-
"""Grid-based binding-pocket volume for the four COX-2 structures.

fpocket is not installable in this environment (no compiler/conda), so we use a
transparent geometric equivalent:

  1. parse protein heavy atoms and the co-crystal ligand from the PDBQT/PDB;
  2. lay a cubic grid over a box around the ligand and flag every voxel inside a
     protein atom's van der Waals radius as occupied;
  3. seed the free voxels adjacent to the ligand and flood-fill the connected
     free region (6-connectivity, so the fill cannot leak through narrow gaps);
  4. pocket volume = reached voxels x voxel volume.

The fill is bounded by the protein vdW surface and the grid box, so it measures
the ligand-connected cavity rather than bulk solvent. Two grid spacings are run
to show the estimate is stable. This is a geometric pocket descriptor, not a
binding free energy.

Path configuration: set CHAIHU_DOCK_DIR to the local docking work root
(default C:\\dockwork); output is written next to this module's data via
CHAIHU_POCKET_OUT if set, otherwise the dockwork root.
"""
import os
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.Chem import GetPeriodicTable

DOCK = Path(os.environ.get("CHAIHU_DOCK_DIR", r"C:\dockwork"))
OUT = Path(os.environ.get("CHAIHU_POCKET_OUT", str(DOCK)))
PT = GetPeriodicTable()

STRUCTURES = {
    "5KIR": {
        "receptor": DOCK / "receptors" / "5KIR.pdbqt",
        "ligand": DOCK / "redock" / "RCX_crystal.pdb",
    },
    "6COX": {
        "receptor": DOCK / "ensemble" / "6COX.pdbqt",
        "ligand": DOCK / "redock2" / "S58_6COX_crystal.pdb",
    },
    "1CX2": {
        "receptor": DOCK / "ens2" / "1CX2_rigid.pdbqt",
        "ligand": DOCK / "redock2" / "S58_1CX2_crystal.pdb",
    },
    "3PGH": {
        "receptor": DOCK / "ens2" / "3PGH_rigid.pdbqt",
        "ligand": DOCK / "redock2" / "FLP_3PGH_crystal.pdb",
    },
}


def infer_element(atom_name):
    name = atom_name.strip().lstrip("0123456789")
    if not name:
        return None
    two = name[:2]
    if two in ("Cl", "Br"):
        return two
    return name[0]


def parse_atoms(path):
    """Return list of (element, x, y, z) heavy atoms from PDB/PDBQT."""
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        elem = infer_element(line[12:16])
        if elem is None or elem == "H":
            continue
        atoms.append((elem, x, y, z))
    return atoms


def rvdw(elem):
    try:
        return float(PT.GetRvdw(PT.GetAtomicNumber(elem)))
    except Exception:
        return 1.70


def pocket_volume(protein, ligand, spacing, margin=8.0):
    lig_xyz = np.array([(x, y, z) for _, x, y, z in ligand])
    lo = lig_xyz.min(axis=0) - margin
    hi = lig_xyz.max(axis=0) + margin
    grid_axes = [np.arange(lo[a], hi[a], spacing) for a in range(3)]
    nx, ny, nz = (len(a) for a in grid_axes)
    occupied = np.zeros((nx, ny, nz), dtype=bool)

    # flag voxels inside each protein atom's vdW radius (local box per atom)
    for elem, ax, ay, az in protein:
        r = rvdw(elem)
        ix0 = max(int(np.searchsorted(grid_axes[0], ax - r) - 1), 0)
        ix1 = min(int(np.searchsorted(grid_axes[0], ax + r) + 1), nx)
        iy0 = max(int(np.searchsorted(grid_axes[1], ay - r) - 1), 0)
        iy1 = min(int(np.searchsorted(grid_axes[1], ay + r) + 1), ny)
        iz0 = max(int(np.searchsorted(grid_axes[2], az - r) - 1), 0)
        iz1 = min(int(np.searchsorted(grid_axes[2], az + r) + 1), nz)
        if ix1 <= ix0 or iy1 <= iy0 or iz1 <= iz0:
            continue
        gx, gy, gz = np.meshgrid(grid_axes[0][ix0:ix1],
                                 grid_axes[1][iy0:iy1],
                                 grid_axes[2][iz0:iz1], indexing="ij")
        d2 = (gx - ax) ** 2 + (gy - ay) ** 2 + (gz - az) ** 2
        occupied[ix0:ix1, iy0:iy1, iz0:iz1] |= d2 <= r ** 2

    # seed free voxels within 1.6 A of any ligand atom
    seed_set = set()
    for _, lx, ly, lz in ligand:
        ix0 = max(int(np.searchsorted(grid_axes[0], lx - 1.6) - 1), 0)
        ix1 = min(int(np.searchsorted(grid_axes[0], lx + 1.6) + 1), nx)
        iy0 = max(int(np.searchsorted(grid_axes[1], ly - 1.6) - 1), 0)
        iy1 = min(int(np.searchsorted(grid_axes[1], ly + 1.6) + 1), ny)
        iz0 = max(int(np.searchsorted(grid_axes[2], lz - 1.6) - 1), 0)
        iz1 = min(int(np.searchsorted(grid_axes[2], lz + 1.6) + 1), nz)
        gx, gy, gz = np.meshgrid(grid_axes[0][ix0:ix1],
                                 grid_axes[1][iy0:iy1],
                                 grid_axes[2][iz0:iz1], indexing="ij")
        d2 = (gx - lx) ** 2 + (gy - ly) ** 2 + (gz - lz) ** 2
        near = np.argwhere(d2 <= 1.6 ** 2)
        for ix, iy, iz in near:
            v = (ix + ix0, iy + iy0, iz + iz0)
            if not occupied[v]:
                seed_set.add(v)

    # flood fill, 6-connectivity
    reached = np.zeros_like(occupied)
    queue = deque(seed_set)
    n_fill = 0
    while queue:
        ix, iy, iz = queue.popleft()
        if reached[ix, iy, iz] or occupied[ix, iy, iz]:
            continue
        reached[ix, iy, iz] = True
        n_fill += 1
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                           (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            jx, jy, jz = ix + dx, iy + dy, iz + dz
            if 0 <= jx < nx and 0 <= jy < ny and 0 <= jz < nz:
                if not reached[jx, jy, jz] and not occupied[jx, jy, jz]:
                    queue.append((jx, jy, jz))

    volume = n_fill * spacing ** 3
    return {
        "spacing": spacing, "grid_voxels": nx * ny * nz,
        "n_filled_voxels": n_fill, "volume_A3": round(volume, 2),
    }


def main():
    rows = []
    for name, paths in STRUCTURES.items():
        protein = parse_atoms(paths["receptor"])
        ligand = parse_atoms(paths["ligand"])
        print(name, "protein atoms", len(protein), "ligand atoms", len(ligand))
        for spacing in (0.5, 0.75):
            res = pocket_volume(protein, ligand, spacing)
            res["structure"] = name
            rows.append(res)
            print(" ", res)
    df = pd.DataFrame(rows)[["structure", "spacing", "grid_voxels",
                             "n_filled_voxels", "volume_A3"]]
    df.to_csv(OUT / "pocket_volume.csv", index=False, encoding="utf-8-sig")
    piv = df.pivot(index="structure", columns="spacing", values="volume_A3")
    piv["ratio_0.75_to_0.5"] = (piv[0.75] / piv[0.5]).round(3)
    print("\n", piv)
    print("\nNote: grid flood-fill geometric estimate, not fpocket; volumes in A^3.")


if __name__ == "__main__":
    main()
