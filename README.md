# Computational mapping of inflammatory nociceptive signaling in Chaihu-Shugan-San

Reproducible code and data guide for the manuscript:

> **A multi-omics and molecular-modeling framework maps inflammatory nociceptive signaling in Chaihu-Shugan-San**
> Fengkai Liu, Tianhui Gao*
> Target journal: *Computers in Biology and Medicine* (Elsevier)
> Affiliation: School of Pharmacy, Qilu Medical University, Zibo, Shandong, China

This repository is the computational companion to the paper. It documents every
public data source, the analysis pipeline from raw data to figures and tables,
fixed random seeds, and the software environment. The study is entirely
computational (public de-identified data, no new human subjects, no animals).

## What the analysis does

The framework links compounds of the seven-herb formula Chaihu-Shugan-San (which
contains vinegar-processed *Cyperi Rhizoma*) to human targets, then asks whether
those targets are (1) differentially expressed in inflammatory/metabolic liver
disease, (2) localized to specific cell types, (3) supported by perturbation and
causal evidence, and (4) physically plausible by molecular docking, molecular
dynamics, and an external benchmark. The IL1B -> PTGS2 -> PGE2 signaling route
itself is already established by earlier experiments; the question this study
addresses is how the formula's compounds map onto that route and how robust the
mapping is. Direct binding and analgesic efficacy are not demonstrated here and
are explicitly left to later experimental work.

## Repository layout

```
chaihu-nociceptive-mapping/
├── README.md                 # this file
├── LICENSE                   # MIT
├── requirements.txt          # pinned Python packages (pip)
├── environment.yml           # equivalent conda environment
├── .gitignore
├── scripts/                  # curated, final analysis scripts by module
│   ├── 01_bulk_deg/          # GSE89632 QC, annotation, differential expression
│   ├── 02_tcmsp_targets/     # TCMSP parsing, UniProt mapping, intersections
│   ├── 03_ppi_enrichment/    # PPI, enrichment, ssGSEA, ROC
│   ├── 04_docking/           # receptor/ligand prep, AutoDock Vina, DUD-E, redock
│   ├── 05_md_mmgbsa/         # OpenMM setup, production, analysis, MM-GBSA
│   ├── 06_multicohort_meta/  # five-cohort voting, random-effects meta-analysis
│   ├── 07_singlecell_spatial/# single-cell and spatial localization
│   ├── 08_perturbation/      # L1000/Enrichr, CellOracle-style perturbation
│   ├── 09_benchmarks/        # dilution, CI/TOST, framework benchmark, MASLD
│   ├── 10_manuscript/        # builds the Word manuscript, letter, highlights
│   ├── 11_mr/                # eQTLGen instruments, MR estimators, ABF colocalization, SMR/HEIDI
│   └── 12_drg/               # human DRG receptor/channel expression bridge
├── data/
│   ├── processed_small/      # small result tables that can be versioned in git
│   └── raw_external_links/   # manifests pointing to raw/GEO data (large files)
├── figures/                  # versioned publication figures (e.g. DRG dot plot)
├── docs/
│   ├── workflow.md           # step-by-step pipeline and figure/table mapping
│   ├── code_map.md           # curated script -> original script -> output mapping
│   ├── data_sources.md       # every database, URL, accession, and version
│   └── 第3/4/6轮 reports   # round-by-round analysis reports (第6轮_防守项/)
├── notebooks/                # optional exploratory notebooks (not authoritative)
└── tests/                    # lightweight sanity checks
```

## Quick start

```bash
# 1. Create the environment
python -m venv .venv
.venv\Scripts\activate            # Windows
python -m pip install -r requirements.txt

# 2. Download/prepare raw data following docs/data_sources.md
#    (large files are not stored in git; manifests record exact URLs and checksums)

# 3. Run modules in numeric order 01 -> 12
#    Each script writes its outputs under a module-specific results directory.
#    Module 11 needs an OpenGWAS token: set OPENGWAS_JWT (or OPENGWAS_JWT_FILE).
```

Molecular docking additionally requires the standalone **AutoDock Vina 1.2.7**
executable (path configured per module), and molecular dynamics uses **OpenMM**.
Ligand protonation uses **Meeko** and **Open Babel**. See `docs/workflow.md`.

## Reproducibility notes

- Random seeds are fixed inside each stochastic script (docking seed `42`;
  numeric and resampling seeds recorded where used).
- Multiple testing uses Benjamini-Hochberg FDR; the meta-analysis uses Hedges' g
  with a DerSimonian-Laird random-effects model.
- Null and inconclusive results are reported as such. The meta-analysis is
  underpowered and inconclusive (confidence intervals cross zero); it is not
  interpreted as evidence of a zero effect.
- MR finds no robust causal effect of the seven genes' expression on abdominal
  pain, and colocalization gives PP.H4 all <= 0.003. SMR/HEIDI agrees: one
  nominal SMR p=0.0246 on a fallback instrument, none survives Bonferroni or BH
  across 28 tests. The human DRG analysis shows neurons carry inflammatory-mediator
  receptors (notably PTGER3) without producing the mediators; both are reported as
  bridging, not causal, evidence.
- Docking conclusions are reported at the level the grid supports. All four PTGS2
  structures are co-crystal complexes (5KIR rofecoxib; 1CX2/6COX SC-558; 3PGH
  flurbiprofen), and co-crystal redock passes on all four, with RMSDs of 1.073,
  0.718, 0.105, and 1.478 A. Seed variation is negligible (mean spread 0.04
  kcal/mol), Vinardo and Vina rank ligands similarly but differ by 2.35 kcal/mol
  in absolute value, cross-structure agreement ranges from 1.00 to below 0.26, and
  the positive control celecoxib outscores the formula compounds in the SC-558
  pockets. Luteolin is placed among known inhibitors, not above them.
- The selection benchmark is reported with its residual circularity: against an
  eight-gene external reference the framework keeps 8/8 with 14.6% precision and a
  permutation p<1e-4 (13.55x enrichment over the 664-gene universe), while the
  degree-based and strict-intersection heuristics recover 1/8 and 0/8. Evidence
  tiers do not predict external confirmation, and the paper says so.
- The MASLD single-cell series is complete (19 case + 1 control, 48,594 cells) but
  the single control is composition-mismatched, so it is used for donor-level
  localization only, never for a disease-versus-normal comparison.
- The spatial analysis is replicated across nine Visium files. SPP1 spatial
  clustering holds on all four PSC011 sections (Moran I 0.31-0.49) and separates
  completely from the four C73 normal sections (0.03-0.12); the pooled four-donor
  normal reference and the mouse MASLD case-control table point the same way.
  Each donor group has one donor, so the sections establish technical
  reproducibility rather than population-level inference, and the paper states
  this. Sparse cytokines keep near-zero Moran values because each spot mixes
  multiple cells, which is treated as a detection limit rather than absence.

## Citation and archived data

A frozen snapshot of this repository is archived on **Zenodo** and assigned a DOI;
the DOI is inserted here once the repository is made public. All source data come
from public repositories named in `docs/data_sources.md`.

## License

Code is released under the MIT License (`LICENSE`). Third-party databases retain
their own licenses; see `docs/data_sources.md`.

## Contact

- Corresponding author: Tianhui Gao, gaotianhui@qlmu.edu.cn
- First author: Fengkai Liu
