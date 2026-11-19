# -*- coding: utf-8 -*-
"""Greedily select distinct binding-site poses for multi-pose MD.

Ranks blind poses by affinity and repeatedly keeps the best pose while
excluding every pose within a fixed-frame separation, giving a set of
representatives from distinct sites. Writes each selected pose as a PDB.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

BLIND = PATHS["dock_root"] / "blind"
OUT = PATHS["processed_small"]
WP10 = PATHS["dock_root"] / "md" / "wp10_multipose"
MIN_SEP = 8.0


def parse_file(path):
    text = path.read_text(encoding="utf-8", errors="ignore").replace("\x00", "")
    poses, cur, affinity = [], [], None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if cur:
                poses.append((affinity, cur))
            cur, affinity = [], None
        elif line.startswith("REMARK VINA RESULT"):
            affinity = float(line.split()[3])
        elif line.startswith("ATOM"):
            cur.append(line)
    if cur:
        poses.append((affinity, cur))
    return poses


def fixed_rmsd(p, q):
    return float(np.sqrt(np.mean(np.sum((p - q) ** 2, axis=1))))


def heavy_coords(lines):
    coords = []
    for l in lines:
        elem = l[76:79].strip() if len(l) >= 77 else l[12:16].strip()
        if elem[:1].upper() == "H":
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(coords)


def std_element(type_token):
    t = type_token.strip()
    if t.startswith("H"):
        return "H"
    if t.startswith("O"):
        return "O"
    return "C"


def write_pdb(lines, path):
    with open(path, "w") as fh:
        for i, l in enumerate(lines, 1):
            type_token = l[76:79] if len(l) >= 77 else " "
            elem = std_element(type_token)
            head = l[:76]
            body = f"{elem:>2} "
            fh.write(head + body + "\n")
        fh.write("END\n")


def main():
    all_poses = []
    for f in sorted(BLIND.glob("blind_*.pdbqt")):
        tag = f.stem.replace("blind_", "")
        for i, (aff, raw) in enumerate(parse_file(f), 1):
            all_poses.append({
                "label": f"{tag}_m{i}", "affinity": aff,
                "raw": raw, "coords": heavy_coords(raw),
            })
    all_poses.sort(key=lambda x: x["affinity"])

    selected, chosen = [], []
    for cand in all_poses:
        if any(fixed_rmsd(cand["coords"], s["coords"]) < MIN_SEP
               for s in selected):
            continue
        selected.append(cand)

    WP10.mkdir(parents=True, exist_ok=True)
    rows = []
    for idx, s in enumerate(selected, 1):
        pdb = WP10 / f"pose_{idx}.pdb"
        write_pdb(s["raw"], pdb)
        c = s["coords"].mean(0)
        rows.append({
            "site_pose": idx, "source_pose": s["label"],
            "affinity": s["affinity"],
            "center_x": round(float(c[0]), 2),
            "center_y": round(float(c[1]), 2),
            "center_z": round(float(c[2]), 2),
            "pdb": str(pdb),
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "selected_distinct_sites.csv", index=False)
    print(f"min_separation={MIN_SEP} A -> {len(selected)} distinct sites")
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
