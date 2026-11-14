# -*- coding: utf-8 -*-
"""Enrichment metrics from benchmark docking: ROC AUC, EF, BEDROC."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

W = Path(__file__).resolve().parent
df = pd.read_csv(W / "docking_results.csv").dropna(subset=["affinity"])
df["is_active"] = (df["cls"] == "active").astype(int)
# score: higher = more likely active; more negative affinity is better
df["score"] = -df["affinity"]

n_act = int(df["is_active"].sum())
n_dec = int((1 - df["is_active"]).sum())

# --- ROC AUC (rank based, Mann-Whitney) ---
def auc_manual(scores, labels):
    pos = scores[labels == 1]
    ranks = pd.Series(scores).rank(method="average")
    rpos = ranks[labels == 1].sum()
    return (rpos - n_act * (n_act + 1) / 2) / (n_act * n_dec)

auc = auc_manual(df["score"].values, df["is_active"].values)

# --- enrichment factor at fraction x ---
ordered = df.sort_values("score", ascending=False).reset_index(drop=True)
def ef_at(frac):
    k = max(1, int(round(frac * len(ordered))))
    top = ordered.head(k)
    a = int(top["is_active"].sum())
    return (a / k) / (n_act / len(ordered)), a, k

efs = {f"EF_{int(f*100)}pct": ef_at(f) for f in [0.01, 0.05, 0.10, 0.20]}

# --- BEDROC (alpha=20) ---
def bedroc(alpha=20.0):
    s = ordered["is_active"].values
    n = len(s)
    w = np.exp(-alpha * np.arange(1, n + 1) / n)
    ra = (s * w).sum() / w.sum()
    rrand = n_act / n
    return ra / rrand if rrand else np.nan

summary = {
    "n_active": n_act, "n_decoy": n_dec,
    "ROC_AUC": float(auc),
    "mean_affinity_active": float(df[df.is_active == 1]["affinity"].mean()),
    "mean_affinity_decoy": float(df[df.is_active == 0]["affinity"].mean()),
    "median_affinity_active": float(df[df.is_active == 1]["affinity"].median()),
    "median_affinity_decoy": float(df[df.is_active == 0]["affinity"].median()),
    "BEDROC_alpha20": float(bedroc()),
}
for k, v in efs.items():
    summary[k + "_EF"] = float(v[0])
    summary[k + "_actives_in_top"] = int(v[1])

json.dump(summary, open(W / "enrichment_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
