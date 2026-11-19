# -*- coding: utf-8 -*-
"""Cluster blind-docking poses and identify distinct binding sites.

Parses all blind PDBQT outputs (null bytes stripped), extracts heavy-atom
coordinates and each pose's Vina affinity, clusters at 2.0 A, and reports
each site with its best affinity and a representative centroid.
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
THRESH = 5.0


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
            elem = line[76:79].strip() if len(line) >= 77 else line[12:16].strip()
            if elem[:1].upper() == "H":
                continue
            cur.append([float(line[30:38]), float(line[38:46]),
                        float(line[46:54])])
    if cur:
        poses.append((affinity, cur))
    return poses


def fixed_rmsd(p, q):
    # No alignment: all poses share the receptor frame, so this retains
    # binding-site position (matches Vina's upper-bound RMSD).
    return float(np.sqrt(np.mean(np.sum((p - q) ** 2, axis=1))))


def main():
    labels, coords, affinities = [], [], []
    for f in sorted(BLIND.glob("blind_*.pdbqt")):
        tag = f.stem.replace("blind_", "")
        for i, (aff, c) in enumerate(parse_file(f), 1):
            labels.append(f"{tag}_m{i}")
            coords.append(np.array(c))
            affinities.append(aff)

    n = len(labels)
    n_heavy = len(coords[0])
    assert all(len(c) == n_heavy for c in coords)
    rmsd = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = fixed_rmsd(coords[i], coords[j])
            rmsd[i, j] = rmsd[j, i] = r

    clusters = [[i] for i in range(n)]
    while len(clusters) > 1:
        best = None
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                dmin = min(rmsd[i, j] for i in clusters[a] for j in clusters[b])
                if best is None or dmin < best[0]:
                    best = (dmin, a, b)
        if best[0] > THRESH:
            break
        _, a, b = best
        clusters[a] = clusters[a] + clusters[b]
        clusters.pop(b)

    rows = []
    for cid, members in enumerate(clusters, 1):
        affs = [(affinities[i], i) for i in members]
        best_aff, best_idx = min(affs)
        site_center = np.mean([coords[i] for i in members], axis=(0, 1))
        for aff, i in affs:
            rows.append({
                "pose": labels[i], "cluster": cid,
                "affinity": aff,
                "is_best_in_cluster": i == best_idx,
                "site_center_x": round(float(site_center[0]), 2),
                "site_center_y": round(float(site_center[1]), 2),
                "site_center_z": round(float(site_center[2]), 2),
            })
    frame = pd.DataFrame(rows).sort_values(["cluster", "affinity"])
    frame.to_csv(OUT / "blind_pose_clusters.csv", index=False)
    pd.DataFrame(rmsd, index=labels, columns=labels).round(3).to_csv(
        OUT / "blind_pose_rmsd_matrix.csv")

    print(f"blind poses={n} heavy_atoms={n_heavy} clusters={len(clusters)}")
    summary = frame.groupby("cluster").agg(
        size=("pose", "count"),
        best_affinity=("affinity", "min"),
        cx=("site_center_x", "first"),
        cy=("site_center_y", "first"),
        cz=("site_center_z", "first")).reset_index()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
