# -*- coding: utf-8 -*-
"""B1: consolidate threshold sensitivity into one redistributable table.

Mechanically merges the GSE130970 TPM threshold table and the GSE89632 /
formula background permutation summary; no values are hand-entered. Inputs and
output are resolved from config. Writes
``data/processed_small/threshold_sensitivity_summary.csv``.

Run:
  python scripts/09_benchmarks/threshold_sensitivity_summary.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

SMALL = PATHS["processed_small"]
G130970 = SMALL / "neutral_threshold_sensitivity.csv"
G89632 = SMALL / "second_background_permutation.json"
OUT = SMALL / "threshold_sensitivity_summary.csv"


def main():
    rows = []

    a = pd.read_csv(G130970)
    for _, r in a.iterrows():
        rows.append({
            "dataset": "GSE130970",
            "threshold_rule": f"TPM >= {r['TPM_threshold']}",
            "background_size": int(r["universe_size"]),
            "precision_pct": round(float(r["precision"]) * 100, 2),
            "empirical_p": float(r["empirical_p"]),
        })

    b = json.load(open(G89632, encoding="utf-8"))
    key2rule = {
        "GSE89632_log2_8": "log2 intensity >= 8",
        "GSE89632_log2_9": "log2 intensity >= 9",
        "GSE89632_log2_10": "log2 intensity >= 10",
        "conditional_664": "formula-derived 664 background",
    }
    for key, rule in key2rule.items():
        rows.append({
            "dataset": "GSE89632" if key.startswith("GSE") else "formula",
            "threshold_rule": rule,
            "background_size": int(b["backgrounds"][key]),
            "precision_pct": round(float(b["precision_by_background"][key]) * 100, 2),
            "empirical_p": float(b["empirical_p_by_background"][key]),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print("saved", OUT)
    print(out.to_string(index=False))
    print("precision range:", out["precision_pct"].min(), "to", out["precision_pct"].max())


if __name__ == "__main__":
    main()
