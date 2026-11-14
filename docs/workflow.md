# Analysis workflow and output mapping

Modules run in numeric order under `scripts/`. Each module is self-contained and
writes to its own results directory. Paths to raw data and external executables
are configured at the top of each script.

## Pipeline order

1. **Bulk differential expression (01_bulk_deg).** Download GSE89632, inspect
   sample groups, annotate probes, run Welch t-tests with Benjamini-Hochberg FDR,
   verify and align the DEG table. GSE130970 provides severity correlations.
2. **Compounds and targets (02_tcmsp_targets).** Parse TCMSP for the seven herbs,
   apply OB>=30 and DL>=0.18, map targets through UniProt (exact/fuzzy/rescue),
   and build the strict three-way and broad intersections.
3. **PPI and enrichment (03_ppi_enrichment).** Build the STRING network, run
   GO/KEGG over-representation, ssGSEA immune fractions, and ROC checks.
4. **Docking (04_docking).** Prepare rigid receptors and ligands, run AutoDock
   Vina (seed 42), validate with a co-crystal redock, run the DUD-E benchmark and
   positive controls, run a four-structure x two-scoring x three-seed grid
   (144 runs) with a celecoxib positive control in every panel, and estimate each
   co-crystal pocket volume with a grid flood-fill at two voxel spacings.
5. **Molecular dynamics (05_md_mmgbsa).** Build the solvated 5KIR-luteolin
   complex, minimize/equilibrate, run production, analyze RMSD/H-bonds, compute
   MM-GBSA, and run within-trajectory robustness checks (blocks, half split,
   drift regression, pose stability, frame bootstrap).
6. **Multi-cohort meta-analysis (06_multicohort_meta).** Harmonize five cohorts,
   run direction voting, Hedges' g random-effects pooling, sensitivity and BH.
7. **Single-cell and spatial (07_singlecell_spatial).** Localize targets across
   cell types and spatial spots, including disease-vs-normal percentage shifts,
   retrieve all 20 GSE344087 samples individually for the full-series map, and
   run the nine-file multi-section Visium analysis (4 PSC011 vs 4 C73 sections
   plus the pooled e84 normal reference) with memory-safe backed reads and a
   mouse MASLD donor-level case-control table. Sections within a donor are
   technical replicates, so these spatial contrasts stay descriptive. The 20-sample
   MASLD series is then integrated with Harmony (sample as batch) and checked with
   a same-sample neighbor mixing metric, marker heatmap, and composition stack.
8. **Perturbation (08_perturbation).** Enrichr/L1000 reversal signatures and
   CellOracle-style in silico knockout/overexpression.
9. **Benchmarks (09_benchmarks).** Bulk-dilution analysis, confidence intervals
   and TOST, framework precision/recall, de-circularized external-reference
   benchmark with an exact permutation null, and the human MASLD localization.
10. **Manuscript build (10_manuscript).** Assemble the Word manuscript, cover
    letter, highlights, and full target list.
11. **Mendelian randomization (11_mr).** Pull LD-clumped eQTLGen instruments,
    classify cis/trans with GRCh37 coordinates, harmonize against four UKB
    abdominal-pain outcomes, run Wald/IVW/Egger/weighted-median and Cochran Q,
    run Wakefield ABF colocalization with prior sensitivity, and run SMR with
    the analytic HEIDI test at the top cis SNP (fallback cis instrument when the
    top SNP is absent), with Bonferroni/BH adjustment across 28 tests.
12. **Human DRG bridge (12_drg).** Download the Nguyen 2021 human DRG h5ad,
    quantify per-cluster receptor/channel expression, build annotation and
    mediator-receptor tables, and render the receptor dot plot.

## Figure mapping

| Figure | Content | Producing module |
|---|---|---|
| Fig 1 | Overall technical workflow | 10_manuscript / e77 |
| Fig 2 | Direction-voting matrix | 06_multicohort_meta |
| Fig 3 | Meta-analysis forest panels | 06_multicohort_meta |
| Fig 4 | PPI network and enrichment | 03_ppi_enrichment |
| Fig 5 | Single-cell localization | 07_singlecell_spatial |
| Fig 6 | Spatial expression | 07_singlecell_spatial |
| Fig 7 | Docking poses/energies | 04_docking |
| Fig 8 | Molecular dynamics | 05_md_mmgbsa |
| Fig 9 | Causal perturbation | 08_perturbation |
| DRG dot plot | Receptor/channel expression across 15 DRG clusters | 12_drg (figures/drg_receptor_dotplot.png) |
| Harmony UMAP | MASLD pre/post-correction UMAP by sample and condition | 07_singlecell_spatial (figures/masld_harmony_umap.png) |
| Marker heatmap | 9 MASLD cell types by marker-gene z-score | 07_singlecell_spatial (figures/masld_marker_heatmap.png) |
| Composition stack | 20 donor cell-type composition | 07_singlecell_spatial (figures/masld_composition_stacked.png) |

## Table mapping

| Table | Content | Producing module |
|---|---|---|
| Table 1 | Cohort characteristics | 06_multicohort_meta |
| Table 2 | Direction voting | 06_multicohort_meta |
| Table 3 | Target evidence grading | 02_tcmsp_targets / 09_benchmarks |
| Table 4 | Docking binding energies | 04_docking |
| Table 5 | Meta effects and BH | 06_multicohort_meta |
| Ensemble docking | docking_ensemble_summary, _ranks, structure_concordance, _stats | 04_docking |
| Pocket volume | pocket_volume (4 structures x 2 spacings) | 04_docking |
| Harmony integration | masld_harmony_summary, masld_harmony_mixing, masld_marker_heatmap | 07_singlecell_spatial |
| MD robustness | md_block_average, md_robustness | 05_md_mmgbsa |
| De-circularized benchmark | framework_decirculation, framework_permutation | 09_benchmarks |
| UniProt adjudication | uniprot_round2_check, _trembl_canonical, _trembl_adjudicated, _name_resolved | 02_tcmsp_targets |
| MASLD 20-sample | masld_full_localization, _composition, _contrast, _summary | 07_singlecell_spatial |
| Multi-section spatial | 空间_全切片_逐切片, _Moran, _组间汇总, 空间_PSC四切片_中期核验 | 07_singlecell_spatial |
| Mouse case-control | 空间_小鼠病例对照_供体级 | 07_singlecell_spatial |
| MR tables | mr_results_table, coloc_table, instrument_table | 11_mr |
| SMR table | smr_heidi_table, smr_heidi_results | 11_mr |
| DRG tables | drg_cluster_annotation, mediator_receptor_mapping, drg_summary | 12_drg |

## Key statistical and reproducibility settings

- All four PTGS2 structures are co-crystal complexes: 5KIR with rofecoxib, 1CX2
  and 6COX with SC-558, 3PGH with flurbiprofen. None is an apo structure.
- Co-crystal redock RMSDs on all four structures: 5KIR 1.073 A, 1CX2 0.718 A,
  6COX 0.105 A, 3PGH 1.478 A (all pass the < 2 A threshold).
- Ensemble grid: 4 structures (5KIR, 6COX, 1CX2, 3PGH) x 2 scorings (vina,
  vinardo) x 3 seeds (42, 2026, 7), exhaustiveness 8, box 22 A, 5 modes; mean
  seed spread 0.04 kcal/mol (max 0.65); Vinardo systematically 2.35 kcal/mol less
  negative than Vina; cross-structure Spearman is 1.00 between the two SC-558
  structures and 0.886 between 1CX2 and 3PGH, dropping to 0.257 or below when
  5KIR enters the comparison; pocket shape differs across the co-crystal forms.
- Celecoxib preparation note: Meeko and Open Babel outputs crash Vina with an
  internal tree error (`src/lib/tree.h:101`); the working ligand file is the
  DUD-E control `c-1.pdbqt`, and this is the only celecoxib file used in the grid.
- Differential expression: Welch t, BH FDR.
- Meta-analysis: Hedges' g, DerSimonian-Laird random effects; equivalence tested
  with TOST (no gene passes at bounds 0.2/0.3).
- MM-GBSA: 125 frames, OBC2, single trajectory, entropy omitted; treated as a
  stability/relative metric only because absolute values are biased negative.
  Five 25-frame blocks span minus 29.76 to minus 28.35 kcal/mol (range 1.40);
  drift after discarding 12 frames is plus 0.22 kcal/mol per ns (p 0.009); the
  ligand stays within 0.3 nm in 89.2% of frames; frame bootstrap 95% interval
  minus 29.49 to minus 28.67. Pose stable, energy not fully converged, no
  independent replicate.
- Selection benchmark: 8 external pain genes (after removing the 6 seeded ones),
  framework precision 14.6% (55 targets) and 16.3% (49), permutation p<1e-4,
  enrichment 13.55x over the 664-gene pool; degree top-60 1.7% (p 0.53); strict
  intersection 0 hits; tier ordering does not survive the external reference
  (core 1/8, strong 1/8).
- MASLD single-cell: 19 MASLD + 1 control, 48,594 cells; donor-level medians
  IL1B 80.0%, TNF 52.6%, PTGS2 44.0% in macrophages/Kupffer; control composition
  not matched, so localization only.
- Multi-section Visium: 9 files (4 PSC011, 4 C73, pooled e84); self-made lognorm
  on in-tissue spots, Delaunay Moran with 200 perms; backed reads with 2,000-gene
  chunks keep peak memory at about 3.2 GB. SPP1 Moran I is 0.31-0.49 on PSC versus
  0.03-0.12 on C73 (no overlap), positive spots 29.3-40.1% versus 1.7-2.5%; one
  donor per group, so technical reproducibility only, no cross-section p-values.
  The pooled e84 normal has no coordinates and cannot return Moran.
- Mouse MASLD Visium: 1 SD vs 2 WD donors, descriptive only; Spp1 4.52x
  (0.4254->1.9241), Ccl2 4.90x, Il1b 1.61x; Il6/Calcr/Trpv1/Tacr1 are all zero.
- Harmony (harmonypy 0.0.6, pure Python): sample as batch, random_state 0;
  same-sample neighbor fraction 0.359 to 0.237 overall (MASLD 0.319 to 0.191,
  HC 0.593 to 0.509). A descriptive mixing check only; one healthy donor means
  no inferential disease comparison. The embeddings .npz and per-cell .csv are
  regenerable and not committed.
- Pocket volume (grid flood-fill, 16_pocket_volume.py): 0.5 A voxel volumes
  5KIR 6967, 6COX 7490, 1CX2 7578, 3PGH 5659 A^3; 0.5/0.75 A spacing agreement
  within 2%. A geometric cavity descriptor, not a binding energy.
- DRG and other external tissue resources are analyzed separately and never
  pooled with liver single-cell data.
- MR: instruments p<=5e-8, LD clump r2=0.001, 5000 kb, EUR; cis within 1 Mb;
  all cis F>=10; no robust causal effect, coloc PP.H4 all <=0.003.
- SMR/HEIDI: single top cis SNP, delta-method SE, region gene body +/-100 kb,
  HEIDI eQTL p<0.05, proxies r2>0.9 removed; one nominal SMR p=0.0246 on a
  fallback instrument, none survives Bonferroni or BH across 28 tests.
- DRG: log-normalized X, 15 clusters; percentages and cluster means computed
  per gene; pipeline output verified byte-identical to the original tables.
