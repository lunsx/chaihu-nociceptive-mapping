"""Lightweight sanity checks for the curated repository.

Run:  python tests/test_sanity.py
These checks verify structure and that the bundled summary files parse; they do
not rerun the analyses or require the large raw data.
"""
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODULES = [
    "01_bulk_deg", "02_tcmsp_targets", "03_ppi_enrichment", "04_docking",
    "05_md_mmgbsa", "06_multicohort_meta", "07_singlecell_spatial",
    "08_perturbation", "09_benchmarks", "10_manuscript",
    "11_mr", "12_drg", "13_cross_formula",
]


def test_module_dirs_and_scripts():
    scripts = ROOT / "scripts"
    for mod in MODULES:
        d = scripts / mod
        assert d.is_dir(), f"missing module {mod}"
        assert any(d.glob("*.py")), f"no scripts in {mod}"


def test_bundled_json_parses():
    small = ROOT / "data" / "processed_small"
    for name in ["dilution_stats.json", "enrichment_summary.json",
                 "controls_luteolin.json", "masld_summary.json",
                 "drg_summary.json", "smr_heidi_results.json",
                 "docking_ensemble_stats.json", "md_robustness.json",
                 "framework_permutation.json", "masld_full_summary.json",
                 "structure_alert_summary.json",
                 "neutral_universe_permutation.json",
                 "heldout_prediction_summary.json",
                 "pseudobulk_titration_summary.json",
                 "independent_correlation.json",
                 "disease_stratified_matrix.json",
                 "lr_curation_summary.json", "mr_power_summary.json",
                 "md_convergence_diagnostic.json",
                 "second_background_permutation.json",
                 "psc_pseudobulk_titration_summary.json",
                 "mouse_pseudobulk_titration_summary.json"]:
        with open(small / name, encoding="utf-8") as fh:
            json.load(fh)


def test_checksums_validate():
    import hashlib
    small = ROOT / "data" / "processed_small"
    lines = (small / "CHECKSUMS.sha256").read_text(encoding="utf-8").strip().splitlines()
    for line in lines:
        digest, fname = line.split("  ", 1)
        h = hashlib.sha256()
        with open(small / fname, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        assert h.hexdigest() == digest, f"checksum mismatch for {fname}"


def test_required_docs_exist():
    for rel in ["README.md", "LICENSE", "requirements.txt", "pyproject.toml",
                "docs/workflow.md", "docs/code_map.md", "docs/data_sources.md"]:
        assert (ROOT / rel).is_file(), f"missing {rel}"


def test_packaged_data_matches_repo():
    import hashlib
    small = ROOT / "data" / "processed_small"
    bundled = ROOT / "src" / "decimap" / "data"
    names = [
        "external_pain_gene_set.csv",
        "gse89632_significant_DEG.csv",
        "chaihu_filtered_gene_universe.csv",
        "jinling_filtered_gene_universe.csv",
    ]
    for name in names:
        def digest(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    h.update(chunk)
            return h.hexdigest()
        assert digest(bundled / name) == digest(small / name), f"drift in {name}"


if __name__ == "__main__":
    failures = 0
    for fn in [test_module_dirs_and_scripts, test_bundled_json_parses,
               test_checksums_validate, test_required_docs_exist,
               test_packaged_data_matches_repo]:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {fn.__name__}: {exc}")
    print("ALL CHECKS PASSED" if failures == 0 else f"{failures} CHECK(S) FAILED")
    raise SystemExit(1 if failures else 0)
