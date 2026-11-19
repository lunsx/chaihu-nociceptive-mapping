# -*- coding: utf-8 -*-
"""Cluster all available 5KIR-luteolin docking poses by heavy-atom RMSD.

Poses are ligand-only PDBQT files in the same receptor docking frame, so
pairwise best-fit (Kabsch) RMSD is computed directly on heavy atoms, then
hierarchical clustering groups them at 2.0 A.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

POSE_DIR = PATHS["dock_root"] / "ens2" / "poses"
OUT = PATHS["processed_small"]
THRESH = 2.0


def parse_models(path):
    poses, cur = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("MODEL"):
            if cur:
                poses.append(cur)
            cur = []
        elif line.startswith("ATOM"):
            cur.append(line)
    if cur:
        poses.append(cur)
    heavy = []
    for block in poses:
        coords, names = [], []
        for l in block:
            elem = l[76:79].strip() if len(l) >= 77 else l[12:16].strip()
            if elem[:1].upper() == "H":
                continue
            coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            names.append(l[12:16].strip())
        heavy.append((np.array(coords), names))
    return heavy


def kabsch_rmsd(p, q):
    pc = p - p.mean(0)
    qc = q - q.mean(0)
    h = pc.T @ qc
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    diag = np.diag([1, 1, d])
    rot = vt.T @ diag @ u.T
    fitted = pc @ rot.T
    return float(np.sqrt(np.mean(np.sum((fitted - qc) ** 2, axis=1))))


def main():
    files = sorted(POSE_DIR.glob("5KIR_luteolin_*.pdbqt"))
    labels, coords = [], []
    for f in files:
        for i, (c, names) in enumerate(parse_models(f), 1):
            labels.append(f"{f.stem}_m{i}")
            coords.append(c)
    extra = PATHS["dock_root"] / "out" / "5KIR_luteolin.pdbqt"
    if extra.exists():
        for i, (c, _) in enumerate(parse_models(extra), 1):
            labels.append(f"out_luteolin_m{i}")
            coords.append(c)

    n = len(labels)
    n_heavy = len(coords[0])
    assert all(len(c) == n_heavy for c in coords), "heavy atom count mismatch"
    rmsd = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = kabsch_rmsd(coords[i], coords[j])
            rmsd[i, j] = rmsd[j, i] = r

    # Agglomerative single-link clustering by repeatedly merging nearest
    # clusters until the nearest distance exceeds the threshold.
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

    assign = []
    for cid, members in enumerate(clusters, 1):
        for idx in members:
            # centroid = member with smallest mean RMSD to the cluster
            d = {i: rmsd[i, members].mean() for i in members}
            centroid = min(d, key=d.get)
            assign.append({
                "pose": labels[idx], "cluster": cid,
                "is_centroid": idx == centroid,
                "mean_rmsd_to_cluster": round(float(d[idx]), 3),
            })
    frame = pd.DataFrame(assign).sort_values(["cluster", "pose"])
    frame.to_csv(OUT / "pose_cluster_rmsd.csv", index=False)
    pd.DataFrame(rmsd, index=labels, columns=labels).round(3).to_csv(
        OUT / "pose_pairwise_rmsd_matrix.csv")

    print(f"poses={n} heavy_atoms={n_heavy} threshold={THRESH} A")
    print(f"clusters={len(clusters)}")
    for cid, members in enumerate(clusters, 1):
        print(f"  cluster {cid}: {len(members)} poses")
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
