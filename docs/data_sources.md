# Data sources

All inputs are public, de-identified data. No new human-subject or animal data
were collected. Local raw downloads and large derived files are intentionally
excluded from git; this file records where each one comes from.

## 1. Bulk transcriptomics (liver)

| Dataset | Role | Key facts | Source |
|---|---|---|---|
| GSE89632 | Primary bulk DEG | 63 samples; NAFLD/NASH vs control; Welch t + BH, 2,982 DEG | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89632 |
| GSE130970 | Severity correlation | 78 cases; graded severity; TPM; IL1B r=+0.55, TNF r=+0.46 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130970 |

## 2. Multi-cohort validation (five cohorts, 332 samples)

Case-control cohorts used in the pooled meta-analysis:

| Dataset | n (case/control) |
|---|---|
| GSE281797 | 64 / 30 |
| GSE260666 | 4 / 6 |
| GSE147304 | 5 / 5 |

Two additional graded-severity cohorts are used for correlation only (not pooled
as case-control). Meta-analysis: Hedges' g, DerSimonian-Laird random effects.
All accessions are resolved through the GEO portal: https://www.ncbi.nlm.nih.gov/geo/

## 3. Single-cell and spatial

| Resource | Role | Key facts / identifier |
|---|---|---|
| Human PSC/PBC single-cell | Cell-type localization | CELLxGENE dataset `ff713455-e6b7-45b3-b2dd-1c9dd26251f6`; 89,637 cells, 16 donors; PSC 47,156, PBC 18,240, normal 24,241 |
| Human MASLD GSE344087 | Correct-disease localization | 20 samples (19 MASLD + 1 HC); all 20 retrieved per sample; 48,594 cells, 9 types after QC; donor-level medians reported |
| Human PSC Visium | Spatial autocorrelation | PSC011 four sections `085648af-…`, `1d929cac-…`, `975d4575-…`, `c8c7681f-…`; SPP1 Moran I 0.402/0.308/0.390/0.485, reproducible across all four sections |
| Human normal Visium | Spatial comparison | C73 four sections `36122151-…`, `878d13cd-…`, `adbebdbd-…`, `c12dd8cc-…`; SPP1 Moran I 0.116/0.053/0.025/0.106, no overlap with the PSC range |
| Multi-donor normal Visium | Independent normal reference | pooled `e84dbcac-…` (4,333 spots, 4 donors H35-H38); no image coordinates retained, so fractions/expression only (Moran not estimable) |
| Mouse MASLD Visium | Cross-species spatial | pooled accession `4a34293d-929f-41cb-bdc1-ff707e3318b6` (6,671 spots, 1 SD + 2 WD donors); Spp1 0.4254->1.9241 (4.52x), descriptive only |

Harmony batch integration (derived from the GSE344087 MASLD series, not a new
dataset): sample treated as the batch, random_state 0. The same-sample neighbor
fraction fell from 0.359 to 0.237 overall (MASLD 0.319 to 0.191, HC 0.593 to
0.509), a descriptive mixing check; with one healthy donor it is not a disease
comparison. Outputs: `masld_harmony_summary.json`, `masld_harmony_mixing.csv`,
`masld_marker_heatmap.csv`, and the three figures. The large
`masld_harmony_embeddings.npz` and per-cell `masld_harmony_obs.csv` are
regenerable intermediates and are not committed.

GSE344087 per-sample matrices: `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9970nnn/GSM<GSM>/suppl/<GSM>_<MASLD_S0xx|HC_S020>_LIVER_filtered_feature_bc_matrix.h5`
(GSM9970276-9970294 = MASLD_S001-S019, GSM9970295 = HC_S020). The series-level FTP is
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE344nnn/GSE344087/suppl/
CELLxGENE: https://cellxgene.cziscience.com/

The full list of 28 candidate Visium records is retained in the project file
`spatial/visium_liver_records.json`; the mouse MASLD accession is confirmed before
public release.

## 4. Compounds and targets

| Source | Role | Key facts |
|---|---|---|
| TCMSP | Formula compounds/targets | 7 herbs; 938 unique compounds; 139 pass OB>=30 and DL>=0.18; 6,712 compound-target pairs; 946 targets |
| UniProt | Human target mapping | reviewed human reference; 693 mapped rows, 253 unmapped, 664 unique human genes (two TrEMBL replacements written back) |
| HERB | Reproduction check | 234/237 reproduced |
| TTD | Target verification | 196 verified |
| STITCH | Interaction network | luteolin 84 edges; no PTGS2 edge |
| SwissTargetPrediction | Target probability | luteolin PTGS2 probability 0.941 |
| GeneCards | Pathway annotation | IL1B -> PTGS2 -> pain-sensitization chain |

- TCMSP: https://www.tcmsp-e.com/
- UniProt: https://rest.uniprot.org/
- HERB: http://herb.ac.cn/
- TTD: https://db.idrblab.net/ttd/
- STITCH: https://stitch-db.org/
- SwissTargetPrediction: http://www.swisstargetprediction.ch/
- GeneCards: https://www.genecards.org/

## 5. Protein network and enrichment

- STRING (PPI): https://string-db.org/  (239 nodes, 5,047 edges)
- GO/KEGG enrichment via clusterProfiler-style over-representation
- Enrichr / L1000 perturbation signatures: https://maayanlab.cloud/Enrichr
- CellOracle gene-perturbation modeling (Python, separate environment)

## 6. Structural data and docking benchmark

| Item | Detail | Source |
|---|---|---|
| TNF structure | 2AZ5 | https://www.rcsb.org/structure/2AZ5 |
| PTGS2 structure | 5KIR, co-crystallized with rofecoxib | https://www.rcsb.org/structure/5KIR |
| PTGS2 structure | 6COX, co-crystallized with SC-558 (I222 form) | https://www.rcsb.org/structure/6COX |
| PTGS2 structure | 1CX2, co-crystallized with SC-558 | https://www.rcsb.org/structure/1CX2 |
| PTGS2 structure | 3PGH, co-crystallized with flurbiprofen | https://www.rcsb.org/structure/3PGH |
| Redock RMSD (all < 2 A) | 5KIR 1.073, 1CX2 0.718, 6COX 0.105, 3PGH 1.478 | data/processed_small/redock2_rmsd.csv |
| SC-558 / flurbiprofen | PubChem CID 1396 / CID 3394 | https://pubchem.ncbi.nlm.nih.gov/ |
| Pocket volume (grid flood-fill) | 0.5 A voxel volumes (A^3): 5KIR 6967, 6COX 7490, 1CX2 7578, 3PGH 5659; 0.5/0.75 A spacing agreement within 2%; geometric cavity descriptor, not an energy | data/processed_small/pocket_volume.csv |
| MD pose clustering (average linkage) | one backbone basin for all 250 frames (threshold 0.16 nm); six ligand pose clusters, dominant 72.8% and secondary 18% (threshold 0.15 nm) | data/processed_small/md_pose_clusters.csv, md_pose_cluster_summary.json |
| Independent MD replicate (seed 2027) | second 10 ns trajectory through the same pipeline; backbone RMSD 0.150 vs 0.148 nm, HB 4.96 vs 4.62, ARG376/ASN375 anchors reproduced, MM-GBSA minus 31.60 vs minus 29.08 (2.53 apart); pose reproducible, absolute energy still uncertain | data/processed_small/md2_replicate_comparison.json |
| DUD-E target | pgh2 (PTGS2), 3LN1 receptor, celecoxib co-crystal | https://dude.docking.org/targets/pgh2 |
| DUD-E ligands | 50 actives + 250 decoys (benchmark subset) | https://dude.docking.org/ |
| Positive controls | celecoxib, diclofenac, ibuprofen, indomethacin | PubChem CIDs 2662, 3033, 3673, 3715 |

PubChem: https://pubchem.ncbi.nlm.nih.gov/
RCSB: https://www.rcsb.org/

## 7. Mendelian randomization and colocalization

| Resource | Role | Key facts |
|---|---|---|
| OpenGWAS | MR API | POST-first at https://api.opengwas.io/api/; token from `OPENGWAS_JWT`, not committed |
| eQTLGen cis-eQTL | Exposure | 7 `eqtl-a-<ENSG>` blood eQTL datasets, European, hg19; N 14,263-31,684; every gene has cis instruments with F>=10 |
| UKB abdominal pain | Outcomes | ukb-b-11413 (primary), 19097, 6608, 6223; European binary log-odds |
| Ensembl GRCh37 | Coordinates | hg19 gene bodies; cis rule = same chromosome within 1 Mb |

Result: no robust causal effect (PTGS2 IVW fixed p=0.041 is outlier/heterogeneity
driven and not reproducible); colocalization PP.H4 all <= 0.003 across the seven
cis regions. SMR/HEIDI gives the same direction: one nominal SMR p=0.0246 on a
fallback instrument, none survives Bonferroni or BH across 28 tests, and HEIDI
shows no heterogeneity where estimable. No liver-tissue eQTL is hosted on
OpenGWAS; the only available eQTL is eQTLGen whole blood. See
`docs/第3轮_MR与共定位正式分析报告.md`.

## 8. Human dorsal root ganglion

| Resource | Role | Key facts |
|---|---|---|
| Nguyen et al. eLife 2021 | Human DRG snRNA-seq | CELLxGENE `6a6e7f1e-26b0-4f0a-ab34-66e9dec3bed2`; 1,837 nuclei, 15 clusters, 6 donors; log-normalized CSR |
| Tavares-Ferreira 2022 | Independent cross-check | spatial transcriptomics, 8 donors; PTGER3 enriched in PENK+ nociceptors |

Result: DRG neurons do not produce inflammatory mediators (PTGS2 0%, IL1B 0%,
IL6 0.05%) but carry their receptors; PTGER3 (EP3) is the main PGE2 receptor
(mean 30%, max 83.9%). No human DRG dataset with a pain/inflammation condition
exists on CELLxGENE or GEO; the only disease-condition DRG is the mouse iPain
Atlas (8-21 GB, not downloadable here), so a pain-state receptor differential
cannot be computed in this environment. See `docs/第4轮_人DRG伤害感受器桥接分析.md`.

## 9. Software (external executables)

- AutoDock Vina 1.2.7 (standalone executable)
- Open Babel (`obabel`)
- Meeko (`mk_prepare_ligand`)
- OpenMM 8.6.1, OpenFF Sage 2.1.0 force field, TIP3P water
- harmonypy 0.0.6 (pure-Python batch correction; the 2.1.0 release needs a C++ toolchain and is not used)
- MDTraj, SciPy (hierarchical average-linkage pose clustering)

## 10. Review-response resources (data rounds 1-8)

| Resource | Role | Key facts / identifier |
|---|---|---|
| PubChem PUG REST | Structure retrieval | 118/139 retained compounds mapped with explicit match status; parent-name and alias recovery; https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest |
| RDKit FilterCatalog | PAINS/BRENK alerts | PAINS_A/B/C and BRENK; zero formal PAINS hits across 118 structures |
| CellPhoneDB static data | Ligand-receptor curation | github.com/ventolab/cellphonedb-data; interaction/protein/complex/gene inputs; 11/11 protein pairs curated |
| Guide to Pharmacology | PGE2 GPCR pharmacology | familyId 58; PGE2 ligandId 1883; EP1-EP4 = PTGER1-PTGER4 with G-protein coupling; https://www.guidetopharmacology.org |
| eQTLGen / GTEx | Power and tissue check | eQTLGen N=31,684 for instrument R2; no liver eQTL available; multi-tissue GTEx feasibility |

The 23 review-response products in `data/processed_small/` are covered by the
regenerated SHA-256 checksum file and the automated sanity tests. Per-round
narratives are in `docs/数据轮1-8_审稿响应/`.

## 11. Second-round reuse demonstration (reviewer response)

Two additions demonstrate that the bias-control procedures transfer rather than
merely asserting it.

| Product | Role | Key facts |
|---|---|---|
| GSE89632 second neutral background | Independent de-circularized background | Illumina GPL14951; probes collapsed to genes; log2 detection at 8/9/10 in >=50% of samples gives 19,670/15,003/11,503-gene background definitions (nested, not independent); precision 24.1%/19.6%/15.8%, all empirical p = 0.0001 (resolution limit) |
| PSC/PBC independent titration | Reuse of the controlled titration | CELLxGENE `ff713455`; 89,637 cells, 27 types; 60 genes, 58 producer-specific; net signal retained 1.0/0.592/0.326/0.14/0.068/0.017/0 across the same fractions; monotone, near-linear |
| Inference ledger | Evidence hierarchy | `inference_ledger.csv` (manuscript Table 6): each inference step with evidence type, status, figure/table, and the closing experiment |

Across the two datasets the five background definitions are the 664 conditional
set, the GSE130970 neutral set, and the three nested GSE89632 log2 thresholds;
the two datasets are the independent units and the three thresholds are nested
variants of one background. Precision stays within 16 to 24% and does not
collapse. Because the platforms and detection rules differ, only precision and
its significance are compared, not absolute background sizes. Outputs:
`second_background_gene_list.csv`, `second_background_permutation.json`,
`psc_pseudobulk_titration.csv`, `psc_pseudobulk_titration_summary.json`,
`severity_correlation_verified.csv`, and `inference_ledger.csv`.

## 12. Third-round additions (cross-species, external healthy reference, threshold summary)

| Product | Role | Key facts |
|---|---|---|
| Mouse pseudobulk titration | Cross-species reuse of the controlled titration | mouse MASH LoxP/SIRT2, 21,715 cells, 14 types; 55 of 60 genes mapped (CXCL11 no ortholog), 53 producer-specific; net signal retained 1.0/0.595/0.340/0.148/0.078/0.014/0; monotone, near-linear, matching the two human collections at 0.5 |
| External healthy human liver reference | Independent normal reference for localization | Aizarani integrated atlas `135bdd30`; 73,295 cells, 4 donors, 12 author types mapped to comparable MASLD labels (the healthy atlas has no mast or plasma cells); macrophage/Kupffer IL1B 2.95% vs 80.0%, TNF 1.07% vs 52.6%, PTGS2 1.33% vs 44.0%; descriptive, no cross-dataset p-values |
| Threshold sensitivity summary | Mechanically merged threshold table | seven rows; precision 15.79 to 24.07% across the two datasets; no hand-entered values |

Outputs: `mouse_pseudobulk_titration.csv`,
`mouse_pseudobulk_titration_summary.json`,
`external_healthy_reference_comparison.csv`,
`threshold_sensitivity_summary.csv`, and the figure
`figures/fig_external_healthy_reference.png` (manuscript Figure 10). The large
source h5ad files are not redistributed; they are documented above and in the
manuscript.

## 13. Fourth-round additions (suspension sensitivity, mapping files, large-file registry)

| Product | Role | Key facts |
|---|---|---|
| Suspension-type sensitivity | Tests whether the healthy-vs-disease contrast depends on whole-cell vs nuclear suspension | nine rows; the highest whole-cell macrophage IL1B is 7.52%, still about an order of magnitude below the MASLD median of 80.0%; IL1B, TNF, and PTGS2 stay low in both suspensions; FOS is strongly suspension-dependent and is reported as such |
| Healthy cell-type mapping | Author-to-comparable label registry | 12 author types with cell/nucleus counts and notes; documents the two labels absent from the healthy atlas (mast, plasma) |
| Healthy per-type gene percentages | Per-type nonzero percentages for 17 genes | 12 rows covering all, cell, and nucleus percentages; supports the suspension sensitivity analysis |
| English 60-gene list | Machine-readable target list used by the cross-species titration | 60 genes with English tier and direction; zero Chinese |
| Large-file registry | Integrity and provenance for h5ad files that are not redistributed | source URL, size, and SHA-256 for the three large h5ad files; a config-based download script is provided |

Outputs: `external_healthy_suspension_sensitivity.csv`,
`healthy_celltype_mapping.csv`,
`healthy_gene_pct_by_author_type.csv`, `target_gene_list_60.csv`,
`figures/fig_suspension_sensitivity.png` (manuscript Supplementary Figure S7),
`data/raw_external_links/large_files_manifest.csv`, and the config-based
scripts in `scripts/07_singlecell_spatial/` and
`scripts/09_benchmarks/threshold_sensitivity_summary.py`. The large h5ad
source files remain outside the repository and are registered with checksums.

## Third-party licenses

Databases and tools are used under their respective public licenses. TCMSP, STRING,
RCSB, PubChem, UniProt, and DUD-E terms apply to the redistributed summaries; only
small derived tables are included in this repository.
