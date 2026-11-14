# -*- coding: utf-8 -*-
"""Multi-conformation x multi-scoring x multi-seed docking robustness scan.

Receptors : 5KIR (rofecoxib), 1CX2 and 6COX (SC-558), 3PGH (flurbiprofen);
            all four are ligand-bound COX-2 co-crystal structures
Ligands   : five formula compounds + celecoxib (positive control)
Scoring   : vina, vinardo
Seeds     : 42, 2026, 7
Each run is appended to results.csv, so an interrupted scan keeps what it finished.

Run:  python run_ensemble.py
"""
import csv
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(r"C:\dockwork\ens2")
VINA = r"C:\dockwork\vina.exe"
OUT = ROOT / "poses"
OUT.mkdir(exist_ok=True)

RECEPTORS = {
    "5KIR": (r"C:\dockwork\receptors\5KIR.pdbqt", (23.459, 20.752, 19.058)),
    "6COX": (r"C:\dockwork\ensemble\6COX.pdbqt", None),
    "1CX2": (str(ROOT / "1CX2_rigid.pdbqt"), (24.26, 21.53, 16.50)),
    "3PGH": (str(ROOT / "3PGH_rigid.pdbqt"), (25.71, 22.23, 15.05)),
}
LIGANDS = {
    "luteolin": r"C:\dockwork\ligands\luteolin.pdbqt",
    "quercetin": r"C:\dockwork\ligands\quercetin.pdbqt",
    "kaempferol": r"C:\dockwork\ligands\kaempferol.pdbqt",
    "cyperotundone": r"C:\dockwork\ligands\cyperotundone.pdbqt",
    "alpha_cyperone": r"C:\dockwork\alpha\alpha_cyperone.pdbqt",
    "celecoxib": r"C:\dockwork\dude_bench\control_pdbqt\c-1.pdbqt",
}
SCORINGS = ["vina", "vinardo"]
SEEDS = [42, 2026, 7]
BOX = 22
EXH = 8

centre_6cox = json.loads((pathlib.Path(r"C:\dockwork\ensemble") / "6COX_center.json")
                         .read_text())
if isinstance(centre_6cox, list):
    centre_6cox = tuple(centre_6cox)
RECEPTORS["6COX"] = (RECEPTORS["6COX"][0], tuple(centre_6cox))

CSV = ROOT / "results.csv"
FIELDS = ["receptor", "ligand", "scoring", "seed", "affinity", "mode2", "n_modes", "seconds"]


def first_affinity(stdout):
    """Parse the vina stdout table: '   1        -7.0      0.000  0.000'."""
    rows = []
    for line in stdout.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+", line)
        if m:
            rows.append((int(m.group(1)), float(m.group(2))))
    rows.sort()
    return (rows[0][1] if rows else None, rows[1][1] if len(rows) > 1 else None, len(rows))


def main():
    done = set()
    if CSV.exists():
        with open(CSV, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                done.add((r["receptor"], r["ligand"], r["scoring"], r["seed"]))
    new_file = not CSV.exists()
    fh = open(CSV, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(fh, fieldnames=FIELDS)
    if new_file:
        writer.writeheader()

    total = len(RECEPTORS) * len(LIGANDS) * len(SCORINGS) * len(SEEDS)
    n = 0
    import time
    for rec, (pdbqt, centre) in RECEPTORS.items():
        for lig, ligpath in LIGANDS.items():
            for sc in SCORINGS:
                for seed in SEEDS:
                    n += 1
                    key = (rec, lig, sc, str(seed))
                    if key in done:
                        continue
                    pose = OUT / f"{rec}_{lig}_{sc}_{seed}.pdbqt"
                    cmd = [VINA, "--receptor", pdbqt, "--ligand", ligpath,
                           "--center_x", str(centre[0]), "--center_y", str(centre[1]),
                           "--center_z", str(centre[2]),
                           "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
                           "--exhaustiveness", str(EXH), "--num_modes", "5",
                           "--seed", str(seed), "--cpu", "4", "--scoring", sc,
                           "--out", str(pose)]
                    t0 = time.time()
                    r = subprocess.run(cmd, capture_output=True, text=True,
                                       encoding="utf-8", errors="replace")
                    dt = round(time.time() - t0, 1)
                    best, second, nm = first_affinity(r.stdout)
                    if best is None:
                        print(f"FAIL {rec}/{lig}/{sc}/{seed}: {r.stderr.strip()[:200]}",
                              file=sys.stderr)
                        continue
                    writer.writerow(dict(receptor=rec, ligand=lig, scoring=sc, seed=seed,
                                         affinity=best, mode2=second, n_modes=nm,
                                         seconds=dt))
                    fh.flush()
                    print(f"[{n}/{total}] {rec} {lig} {sc} seed{seed} -> {best} ({dt}s)")
    fh.close()
    print("done")


if __name__ == "__main__":
    main()
