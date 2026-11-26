# -*- coding: utf-8 -*-
"""Targeted re-docking with the Vinardo scoring function (same pocket box).

Runs the identical pocket-sized search as the Vina targeted job but switches
scoring function, giving a local, reproducible multi-scoring consensus on
whether the canonical site is dockable without a GPU or large download.
"""
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
OUT_PDBQT = OUT_DIR / "targeted_vinardo.pdbqt"
LOG = OUT_DIR / "targeted_vinardo.console.txt"

CENTER = (32.58, 20.05, 25.07)
BOX = 25.0


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PATHS["vina"]),
        "--receptor", str(PATHS["dock_root"] / "receptors" / "5KIR.pdbqt"),
        "--ligand", str(PATHS["dock_root"] / "ligands" / "luteolin.pdbqt"),
        "--scoring", "vinardo",
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
        "--exhaustiveness", "32", "--num_modes", "20",
        "--energy_range", "8", "--min_rmsd", "1.5", "--cpu", "8",
        "--out", str(OUT_PDBQT),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    LOG.write_text(proc.stdout + proc.stderr, encoding="utf-8")
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"vina failed rc={proc.returncode}")


if __name__ == "__main__":
    main()
