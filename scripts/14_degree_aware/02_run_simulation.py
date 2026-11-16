# -*- coding: utf-8 -*-
"""Run the degree-aware simulation benchmark and save the rejection-rate table.

Writes data/processed_small/simulation_validation.csv. The empirical rejection
rate is the type I error under S0/S1 and power under S2.
"""
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

sys.path.insert(0, str(PATHS["repo_root"] / "src"))
from decimap.simulation import run_simulation  # noqa: E402

result = run_simulation()
out = PATHS["processed_small"] / "simulation_validation.csv"
result.to_csv(out, index=False, encoding="utf-8")
print(result.to_string(index=False))
print("saved:", out)
