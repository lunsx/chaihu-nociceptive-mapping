# -*- coding: utf-8 -*-
"""Match P2Rank pockets to the canonical luteolin site.

For each predicted pocket we report the centre distance to the canonical
ligand centroid and the overlap with the canonical contact residues reported
in WP7 (PHE142, HIS226, GLN374, ASN375, ARG376). This establishes whether an
independent pocket detector recognises the canonical site and how it ranks it.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
P2RANK_CSV = OUT_DIR / "p2rank_out" / "5KIR_monomer.pdb_predictions.csv"
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
    for line in P2RANK_CSV.read_text(encoding="utf-8").splitlines():
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
    frame.to_csv(OUT_DIR / "pocket_match_summary.csv", index=False)
    print(frame.to_string(index=False))
    matched = frame[frame.n_canonical_residues_hit >= 3]
    print("\nPockets matching the canonical site (>=3 contact residues):")
    print(matched.to_string(index=False))


if __name__ == "__main__":
    main()
