# DECIMAP: de-circular, evidence-graded computational mapping

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24283746.svg)](https://doi.org/10.5281/zenodo.24283746)

Reproducible code, data, and a command line tool for the manuscript:

> **DECIMAP: a de-circularized, evidence-graded computational framework for mapping Chinese medicine formulas onto a disease axis**
> Fengkai Liu, Tianhui Gao*
> Target journal: *Computer Methods and Programs in Biomedicine* (Elsevier)
> Affiliation: School of Pharmacy, Qilu Medical University, Zibo, Shandong, China

DECIMAP is both the name of the framework and an installable Python package. It
packages the reusable bias-control procedures used in the paper, an external-only
de-circularized enrichment benchmark and a controlled pseudobulk titration, and
it applies them to two formulas on the same liver inflammation to pain-sensitization
axis. Version 2 adds a degree-aware estimator: a fixed-size null that samples
targets in proportion to their network degree (a Pareto key model), replacing the
uniform permutation, with a bias-corrected enrichment ratio, negative-control
empirical calibration, a conditional portability metric, and a synthetic-graph
simulation that verifies type I error control and power. The study is entirely
computational (public de-identified data, no new human subjects, no animals, no
cell experiments).

## Install the tool

```bash
python -m venv .venv
.venv\Scripts\activate                  # Windows
python -m pip install .
```

The package depends only on numpy, pandas, and scipy and bundles the small
reference files it needs, so the benchmark runs offline in a clean install.

```bash
decimap benchmark              # run the symmetric de-circularized benchmark
decimap benchmark --output metrics.csv
decimap degree                 # run the degree-aware de-circularized enrichment
decimap calibrate              # negative-control empirical calibration
decimap portability            # conditional cross-formula portability
decimap simulate               # synthetic-graph type I error and power benchmark
decimap validate               # verify bundled data and reproduce reference values
decimap --version
```

`decimap validate` recomputes the benchmark and checks it against the values
reported in the manuscript. It also verifies the integrity of the bundled data.
The degree-aware commands read the bundled compound-gene edges and degree tables
and run offline in a clean install.

## What the analysis does

The framework maps compounds of a formula onto human targets, then asks whether
those targets are (1) differentially expressed in inflammatory and metabolic
liver disease, (2) localized to specific cell types, (3) supported by perturbation
and human genetic evidence, and (4) physically plausible by molecular docking and
molecular dynamics. The IL1B -> PTGS2 -> PGE2 route is established by earlier
experiments; the contribution here is a reusable, de-circularized way to show how
formula compounds map onto that route, and to quantify how strongly the mapping
enriches against external references. Direct binding, causal effects, and analgesic
efficacy are not demonstrated and are left to later experimental work.

Two formulas are studied on the same disease axis. Chaihu Shugan San is the first
case and Jinlingzi San is the portability test; the two share no herbs but their
filtered genes overlap at Jaccard 0.703, so the second case is presented as
portability and convergence, not an independent gene set. An exploratory third
formula, Xiaoyao San, is kept as a control.

## Repository layout

```
chaihu-nociceptive-mapping/
├── README.md                 # this file
├── LICENSE                   # MIT
├── pyproject.toml            # DECIMAP package build and console entry point
├── requirements.txt          # pinned packages for the full analysis environment
├── environment.yml           # equivalent conda environment
├── src/
│   └── decimap/              # installable package (cli, benchmark, degree_aware,
│                              # simulation, validate, data)
├── scripts/                  # curated, final analysis scripts by module
│   ├── 01_bulk_deg/          # GSE89632 QC, annotation, differential expression
│   ├── 02_tcmsp_targets/     # TCMSP parsing, UniProt mapping, intersections
│   ├── 03_ppi_enrichment/    # PPI, enrichment, ssGSEA, ROC
│   ├── 04_docking/           # receptor/ligand prep, AutoDock Vina, DUD-E, redock
│   ├── 05_md_mmgbsa/         # OpenMM setup, production, analysis, MM-GBSA
│   ├── 06_multicohort_meta/  # five-cohort voting, random-effects meta-analysis
│   ├── 07_singlecell_spatial/# single-cell and spatial localization
│   ├── 08_perturbation/      # L1000/Enrichr and in-silico perturbation
│   ├── 09_benchmarks/        # dilution, CI/TOST, framework benchmark, MASLD
│   ├── 10_manuscript/        # builds the Word manuscript, letter, highlights
│   ├── 11_mr/                # eQTLGen instruments, MR, colocalization, SMR/HEIDI
│   ├── 12_drg/               # human DRG receptor/channel expression bridge
│   ├── 13_cross_formula/     # self-contained cross-formula fetch, parse, benchmark
│   └── 14_degree_aware/      # build compound-gene edges and run the simulation
├── data/
│   ├── processed_small/      # small result tables versioned in git
│   └── raw_external_links/   # manifests pointing to raw/GEO data (large files)
├── figures/                  # versioned publication figures
├── docs/                     # workflow, code map, data sources, methodology reports
├── notebooks/                # optional exploratory notebooks (not authoritative)
└── tests/                    # lightweight sanity checks
```

## Reproduce the full pipeline

```bash
# Full environment (the package alone needs only numpy/pandas/scipy)
python -m pip install -r requirements.txt

# Run modules in numeric order 01 -> 13. Each script writes its own outputs.
# Module 11 needs an OpenGWAS token: set OPENGWAS_JWT (or OPENGWAS_JWT_FILE).
python tests/test_sanity.py
```

Molecular docking requires the standalone AutoDock Vina executable and molecular
dynamics uses OpenMM; ligand preparation uses Meeko and Open Babel. See
`docs/workflow.md` and `docs/data_sources.md`.

## Reproducibility and honest boundaries

- Random seeds are fixed in every stochastic script (permutation seed 20261116;
docking seed 42). Multiple testing uses Benjamini-Hochberg FDR.
- The symmetric benchmark returns near-identical precision for the two formulas,
11.76% and 11.90%, each about three-fold above the permutation null, with Fisher
p = 0.0026 and 0.0103. The two cases are not independent confirmations because
the gene universes overlap.
- The degree-aware estimator is validated on synthetic graphs with known truth.
When both the selected and reference sets independently favour high-degree hubs,
the uniform permutation rejects 97.6% of the time (false positives), while the
degree-aware null rejects 1.8% and controls type I error; it retains power for
genuine signal planted on low-degree genes (0.29 to 0.73 as the number of planted
genes rises from 2 to 8). On the real data the degree correction roughly halves
the apparent enrichment (Chaihu 3.32 to 1.76; Jinling 3.04 to 1.56), which stays
significant.
- Co-crystal redock passes on all four PTGS2 structures, with RMSDs of 1.073,
0.718, 0.105, and 1.478 A, and the DUD-E benchmark confirms the docking protocol
enriches active ligands. Docking scores are not binding free energies, and the
MM-GBSA values do not converge across technical replicates, so they are used only
relatively.
- The single fully held-out cohort does not support prospective direction
prediction (42.3% agreement, below chance). This is reported as a failed end-to-end
generalization test, distinct from the component-level controls that pass.
- MR, colocalization, and SMR do not support a causal claim, and the instruments
derive from whole-blood rather than liver eQTL. The null is not evidence of absence.
- The MASLD single-cell series has one healthy control, so it supports localization
but not disease-contrast inference. Spatial sections are technical replicates from
one donor per group and are treated as descriptive.

## Citation and archived data

A frozen snapshot of each release is archived on Zenodo. Version 2.0.0 is citable
as 10.5281/zenodo.24283746 (https://doi.org/10.5281/zenodo.24283746). The all-versions
concept DOI 10.5281/zenodo.24283745 always resolves to the latest release, including
version 2.1.0, which adds the independent-replicate molecular dynamics robustness
analysis. All source data come from public repositories named in `docs/data_sources.md`.

## License

Code is released under the MIT License (`LICENSE`). Third-party databases retain
their own licenses; see `docs/data_sources.md`.

## Contact

- Corresponding author: Tianhui Gao, gaotianhui@qlmu.edu.cn
- First author: Fengkai Liu
