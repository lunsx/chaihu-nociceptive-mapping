# -*- coding: utf-8 -*-
"""Fixed-frame RMSD of docked poses against the canonical luteolin pose.

Heavy atoms are matched optimally with a Hungarian assignment restricted to
the same element, so atom ordering differences between SDF and PDBQT do not
matter. No structural alignment is performed: every pose shares the receptor
frame, so this RMSD measures binding-site position as well as internal pose.
Both the targeted re-dock and the whole-protein blind poses are evaluated.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

OUT_DIR = PATHS["dock_root"] / "pocket_validation"
BLIND = PATHS["dock_root"] / "blind"
CANON_SDF = PATHS["dock_root"] / "md" / "luteolin_docked.sdf"
OUT_CSV = OUT_DIR / "pose_rmsd_summary.csv"


def load_canonical():
    elems, coords = [], []
    lines = CANON_SDF.read_text(encoding="utf-8").splitlines()
    n_atoms = 0
    atom_start = 0
    for i, line in enumerate(lines):
        if line.endswith("V2000"):
            n_atoms = int(line[0:3])
            atom_start = i + 1
            break
    for line in lines[atom_start:atom_start + n_atoms]:
        x = float(line[0:10]); y = float(line[10:20]); z = float(line[20:30])
        elem = line[31:34].strip()
        if elem == "H":
            continue
        elems.append(elem); coords.append((x, y, z))
    return elems, np.array(coords)


def parse_pdbqt(path):
    text = path.read_text(encoding="utf-8", errors="ignore").replace("\x00", "")
    models, cur, aff = [], [], None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if cur:
                models.append((aff, cur))
            cur, aff = [], None
        elif line.startswith("REMARK VINA RESULT"):
            aff = float(line.split()[3])
        elif line.startswith("ATOM"):
            name = line[12:16].strip()
            elem = name.lstrip("0123456789")[:1].upper()
            if elem == "H":
                continue
            cur.append((elem, [float(line[30:38]), float(line[38:46]),
                               float(line[46:54])]))
    if cur:
        models.append((aff, cur))
    return models


def fixed_rmsd(ref_elems, ref_xyz, pose):
    p_elems = [e for e, _ in pose]
    p_xyz = np.array([c for _, c in pose])
    assert len(p_elems) == len(ref_elems), (len(p_elems), len(ref_elems))
    cost = np.zeros((len(ref_elems), len(p_elems)))
    for i, re in enumerate(ref_elems):
        for j, pe in enumerate(p_elems):
            d = np.linalg.norm(ref_xyz[i] - p_xyz[j])
            cost[i, j] = d if re == pe else 1e6
    ri, ci = linear_sum_assignment(cost)
    assert np.all(cost[ri, ci] < 1e6), "no element-consistent bijection"
    d2 = np.sum((ref_xyz[ri] - p_xyz[ci]) ** 2)
    return float(np.sqrt(d2 / len(ri)))


def main():
    ref_elems, ref_xyz = load_canonical()
    ref_centroid = ref_xyz.mean(axis=0)
    rows = []
    targets = [
        ("targeted_vina", OUT_DIR / "targeted_redock.pdbqt"),
        ("targeted_vinardo", OUT_DIR / "targeted_vinardo.pdbqt"),
    ]
    for f in sorted(BLIND.glob("blind_*.pdbqt")):
        targets.append((f"blind:{f.stem}", f))
    for source, path in targets:
        if not path.exists():
            continue
        for mid, (aff, pose) in enumerate(parse_pdbqt(path), 1):
            rmsd = fixed_rmsd(ref_elems, ref_xyz, pose)
            p_xyz = np.array([c for _, c in pose])
            centroid_dist = float(np.linalg.norm(
                p_xyz.mean(axis=0) - ref_centroid))
            rows.append({
                "source": source, "model": mid, "affinity": aff,
                "rmsd_to_canonical": round(rmsd, 3),
                "centroid_distance": round(centroid_dist, 3),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_CSV, index=False)
    print(f"reference heavy atoms={len(ref_elems)} centroid="
          f"{np.round(ref_centroid, 3).tolist()}")
    print(frame.to_string(index=False))
    for src in ("targeted_vina", "targeted_vinardo"):
        t = frame[frame.source == src]
        if len(t):
            best = t.loc[t.affinity.idxmin()]
            print(f"\nbest {src}: model={best.model} aff={best.affinity} "
                  f"RMSD={best.rmsd_to_canonical} centroid_dist="
                  f"{best.centroid_distance}")


if __name__ == "__main__":
    main()
