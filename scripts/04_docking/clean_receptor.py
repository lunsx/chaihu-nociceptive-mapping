# -*- coding: utf-8 -*-
"""Strip the obabel torsion tree to produce a rigid receptor (ATOM records only)."""
from pathlib import Path

W = Path(__file__).resolve().parent
src = W / "receptor.pdbqt"
keep = []
for line in src.read_text(errors="replace").splitlines():
    if line.startswith(("ATOM", "HETATM")):
        keep.append(line)
(W / "receptor_rigid.pdbqt").write_text("\n".join(keep) + "\n")
print("rigid receptor atoms:", len(keep))
