# -*- coding: utf-8 -*-
"""Cross-version pocket-detection robustness with P2Rank 2.5.1.

The 2.3.1 analysis in 04_pocket_match.py is repeated with the latest P2Rank
release (v3 model) to verify that the canonical pocket ranking is not an
artefact of one model version.

P2Rank 2.5.1 requires Java 17 or newer. With a current JDK on the workstation
the detector was run as:

    set JAVA_HOME=C:\\Program Files\\Java\\jdk-25.0.5.1
    prank.bat predict -f 5KIR_monomer.pdb -o p2rank251_out

This script reads the resulting predictions and applies the same centre-distance
and contact-residue matching criteria as 04_pocket_match.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
P2RANK251_CSV = OUT_DIR / "p2rank251_out" / "5KIR_monomer.pdb_predictions.csv"
CANON_CENTROID = np.array([32.576, 20.047, 25.066])
CANON_RESIDUES = {142, 226, 374, 375, 376}


def parse_residue_ids(field):
    ids = set()
    for token in field.split():
        if token.startswith("A_"):
            num = token[2:].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            if num.isdigit():
                ids.add(int(num))
    return ids


def main():
    rows = []
    header = None
    for line in P2RANK251_CSV.read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split(",")]
        if header is None:
            header = parts
            continue
        if not parts[0].startswith("pocket"):
            break
        rec = dict(zip(header, parts))
        center = np.array([float(rec["center_x"]), float(rec["center_y"]),
                           float(rec["center_z"])])
        residues = parse_residue_ids(rec["residue_ids"])
        overlap = sorted(residues & CANON_RESIDUES)
        rows.append({
            "pocket": rec["name"],
            "rank": int(rec["rank"]),
            "probability": float(rec["probability"]),
            "center_distance_to_canonical": round(
                float(np.linalg.norm(center - CANON_CENTROID)), 3),
            "n_canonical_residues_hit": len(overlap),
            "canonical_residues": " ".join(str(r) for r in overlap),
        })

    frame = pd.DataFrame(rows)
    out_csv = OUT_DIR / "pocket_match_251.csv"
    frame.to_csv(out_csv, index=False)
    print(frame.to_string(index=False))
    matched = frame[frame.n_canonical_residues_hit >= 3]
    print("\nP2Rank 2.5.1 pocket matching the canonical site:")
    print(matched.to_string(index=False))


if __name__ == "__main__":
    main()
