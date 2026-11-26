# -*- coding: utf-8 -*-
"""Convert the blind-docking receptor PDBQT into a standard PDB for P2Rank.

Only ATOM/HETATM protein records are kept. Coordinates, chain and residue
numbering are preserved verbatim; the element is derived from the atom name
and written into the standard PDB element columns. The output therefore
represents the exact single-chain receptor used in whole-protein blind docking
(5400 atoms, chain A), rather than the dimer reference structure.
"""
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

SRC = PATHS["dock_root"] / "receptors" / "5KIR.pdbqt"
OUT_DIR = PATHS["dock_root"] / "pocket_validation"
OUT = OUT_DIR / "5KIR_monomer.pdb"


def element_of(atom_name: str) -> str:
    token = atom_name.strip().lstrip("0123456789")
    return token[:1].upper() if token else "C"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = SRC.read_text(encoding="utf-8", errors="ignore").splitlines()
    serial = 0
    out_lines = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        record = line[0:6]
        name = line[12:16]
        altloc = line[16:17]
        resname = line[17:20]
        chain = line[21:22]
        resseq = line[22:26]
        icode = line[26:27]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occ = 1.00
        bfac = 0.00
        elem = element_of(name)
        serial += 1
        out_lines.append(
            f"{record:<6}{serial:5d} {name:<4}{altloc:1}{resname:>3} "
            f"{chain:1}{resseq:4}{icode:1}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}"
            f"{'':10}{elem:>2}"
        )
    out_lines.append("END")
    OUT.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} atoms={serial}")


if __name__ == "__main__":
    main()
