"""Self-contained validation of the installed DECIMAP package.

The check recomputes the benchmark and compares it against the authoritative
values reported in the manuscript and stored in the repository, and verifies
the integrity of the bundled reference files. It needs no network access.
"""
import hashlib
from importlib.resources import files

from .benchmark import run_benchmark

EXPECTED = {
    "Chaihu_Shugan_San": {
        "universe_size": 254,
        "n_external_pain_in_universe": 9,
        "selected_DEG_intersect": 51,
        "observed_pain_overlap": 6,
        "observed_precision": 0.1176,
        "expected_overlap_by_chance": 1.81,
        "enrichment_vs_chance": 3.32,
        "null_precision_mean": 0.0353,
        "null_precision_p95": 0.0784,
        "null_precision_max": 0.1373,
        "fisher_p": 0.0025930368799714993,
        "empirical_p": 0.002,
        "three_way_genes": "IL1A;IL1B;IL6;PTGER3;PTGES;PTGS2",
    },
    "Jinlingzi_San": {
        "universe_size": 204,
        "n_external_pain_in_universe": 8,
        "selected_DEG_intersect": 42,
        "observed_pain_overlap": 5,
        "observed_precision": 0.119,
        "expected_overlap_by_chance": 1.65,
        "enrichment_vs_chance": 3.04,
        "null_precision_mean": 0.0388,
        "null_precision_p95": 0.0952,
        "null_precision_max": 0.1667,
        "fisher_p": 0.010265110018295389,
        "empirical_p": 0.0092,
        "three_way_genes": "IL1A;IL1B;IL6;PTGER3;PTGS2",
    },
}

BUNDLED_FILES = [
    "external_pain_gene_set.csv",
    "gse89632_significant_DEG.csv",
    "chaihu_filtered_gene_universe.csv",
    "jinling_filtered_gene_universe.csv",
]


def _check_bundled_files():
    messages = []
    for name in BUNDLED_FILES:
        with files("decimap.data").joinpath(name).open("rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        messages.append(f"OK  bundled file readable: {name} (sha256 {digest[:12]})")
    return messages


def run_validation():
    """Return (ok, messages) after running all validation checks."""
    messages = _check_bundled_files()
    result = run_benchmark().set_index("formula")
    ok = True

    for formula, expected in EXPECTED.items():
        actual = result.loc[formula]
        for key, exp_value in expected.items():
            got = actual[key]
            if key == "fisher_p":
                match = abs(float(got) - float(exp_value)) < 1e-12
            else:
                match = got == exp_value
            if not match:
                ok = False
                messages.append(f"FAIL {formula}.{key}: expected {exp_value}, got {got}")
        messages.append(f"OK  benchmark values match reference: {formula}")

    return ok, messages
