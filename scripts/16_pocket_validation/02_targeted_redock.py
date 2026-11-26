# -*- coding: utf-8 -*-
"""Targeted re-docking of luteolin into the canonical pocket on the monomer.

A small box is centred on the canonical ligand centroid (32.58, 20.05, 25.07),
using the same single-chain receptor and ligand as whole-protein blind docking
but with a local, pocket-sized search box. This isolates whether the canonical
site is dockable when the search is not spread over the whole protein.
"""
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
RECEPTOR = PATHS["dock_root"] / "receptors" / "5KIR.pdbqt"
LIGAND = PATHS["dock_root"] / "ligands" / "luteolin.pdbqt"
OUT_PDBQT = OUT_DIR / "targeted_redock.pdbqt"
LOG = OUT_DIR / "targeted_redock.console.txt"

CENTER = (32.58, 20.05, 25.07)
BOX = 25.0


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PATHS["vina"]),
        "--receptor", str(RECEPTOR),
        "--ligand", str(LIGAND),
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX),
        "--size_y", str(BOX),
        "--size_z", str(BOX),
        "--exhaustiveness", "32",
        "--num_modes", "20",
        "--energy_range", "8",
        "--min_rmsd", "1.5",
        "--cpu", "8",
        "--out", str(OUT_PDBQT),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    LOG.write_text(proc.stdout + proc.stderr, encoding="utf-8")
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"vina failed rc={proc.returncode}")


if __name__ == "__main__":
    main()
