# Code map: curated script -> original -> output

The repository keeps the final, canonical scripts. The full project additionally
contains many exploratory probes (one-off `c*`, `d6*`, `e*` investigation files);
those are not needed to reproduce reported results and stay in the project archive.

| Curated module | Original script (project) | Purpose / output |
|---|---|---|
| 01_bulk_deg | step1_download.py, step2b_groups.py, step3_annotation.py, step4_deg.py, step4b_verify.py, step4c_align.py | GSE89632 download, groups, annotation, DEG table, verification |
| 02_tcmsp_targets | c14_uniprot.py, c15_map_targets.py, c16_map_v2.py, c17_intersect.py, c18_master.py | UniProt mapping, strict/broad intersections, master target table |
| 02_tcmsp_targets | uniprot_round2_check.py, uniprot_trembl_canonical.py, uniprot_trembl_adjudicate.py, uniprot_name_resolve.py | NO66/DHRS4 resolution, 18 TrEMBL adjudications (REPLACE/AMBIGUOUS/REJECT/KEEP), protein-name resolution for symbol-less targets |
| 02_tcmsp_targets | 09_writeback_trembl.py | writes the two REPLACE accessions back into the master mapping table (693 mapped/253 unmapped/664 unique genes), emits uniprot_TrEMBL_final_verdict.csv |
| 03_ppi_enrichment | d2_ppi.py, d3_enrich.py, d4_roc.py, d5_ssgsea.py | PPI, GO/KEGG enrichment, ROC, immune ssGSEA |
| 04_docking | d8e_receptors.py, d8f_ligands.py, d8h2_clean_rec.py | receptor/ligand preparation, rigid receptor |
| 04_docking | e4_redock_prep.py, e5_redock_ligand.py, e6_meeko.py, e7_rmsd.py | co-crystal redock and RMSD validation |
| 04_docking | e8_alpha_cyperone.py, e9_alpha_prep.py | alpha-cyperone docking |
| 04_docking | e11_ensemble_get.py, e12_ensemble_clean.py, e13_ensemble_compare.py | multi-structure ensemble docking |
| 04_docking | run_ensemble.py, analyze_ensemble.py (C:\dockwork\ens2) | four-structure two-scoring three-seed grid (144 runs), seed/scoring/structural concordance |
| 04_docking | 14_redock_pipeline.py, 15_recompute_rmsd.py (C:\dockwork\redock2) | three-structure co-crystal redock (1CX2/6COX/3PGH), element/atom-order repair, GetBestRMS alignment, redock2_rmsd.csv |
| 04_docking | 16_pocket_volume.py (C:\dockwork) | grid-based flood-fill pocket volume around each co-crystal ligand at two voxel spacings; geometric cavity descriptor, not an energy; outputs pocket_volume.csv |
| 04_docking | run_docking.py, enrichment.py, controls_luteolin.py, clean_receptor.py (C:\dockwork\dude_bench) | DUD-E screen, enrichment metrics, controls |
| 05_md_mmgbsa | e33_system_build.py, e35_min_eq.py, e38_production.py, e39_md_analysis.py, e40_mmgbsa.py | OpenMM system, production, MD analysis, MM-GBSA |
| 05_md_mmgbsa | md_robustness.py | block means, half-trajectory test, drift regression, pose stability, frame bootstrap (no new trajectory) |
| 05_md_mmgbsa | 07_pose_clustering.py (C:\dockwork\md) | average-linkage clustering of backbone-superposed production frames on ligand heavy-atom RMSD (0.15 nm) and backbone RMSD (0.16 nm); outputs md_pose_clusters.csv, md_pose_cluster_summary.json, figures/md_pose_clustering.png |
| 05_md_mmgbsa | 08_md2_analysis.py (C:\dockwork\md2) | independent replicate analysis with the same RMSD/HB/distance and single-trajectory MM-GBSA (OBC2, NoCutoff, entropy omitted) pipeline, plus two-trajectory frame-by-frame comparison; outputs md2_* CSVs and md2_replicate_comparison.json |
| 06_multicohort_meta | e45_download_cohorts.py, e47_build_meta.py, e48_multicohort.py, e49_meta.py, e50_sensitivity.py, e71_meta_bh.py | cohort harmonization, voting, meta, sensitivity, BH |
| 07_singlecell_spatial | d9l_localize.py, d9m_figs.py | target localization and figures |
| 07_singlecell_spatial | download_gsm.py, analyze_masld_full.py (C:\dockwork\masld) | per-sample GSE344087 retrieval (19 MASLD + 1 control), full-series localization |
| 07_singlecell_spatial | all_sections_analysis.py (cloud spatial) | nine-file multi-section Visium analysis (4 PSC011 + 4 C73 + pooled e84), memory-safe backed reads, per-section fractions/Moran and group contrast; outputs 空间_全切片_*.csv |
| 07_singlecell_spatial | mouse_case_control.py (cloud spatial) | mouse MASLD donor-level descriptive case-control (1 SD vs 2 WD); outputs 空间_小鼠病例对照_供体级.csv |
| 07_singlecell_spatial | 10_harmony_integration.py (C:\dockwork\masld) | Harmony batch integration on the 20-sample MASLD series (sample as batch), rebuilt neighbors/UMAP, same-sample neighbor mixing metric, marker heatmap and composition stack; outputs masld_harmony_*.json/csv and figures |
| 08_perturbation | e63_enrichr_lincs.py, e64_summarize.py (lincs/) | Enrichr/L1000 signatures and summary |
| 09_benchmarks | dilution_analysis.py, meta_ci_analysis.py, framework_benchmark.py | dilution, CI/TOST, framework benchmark |
| 09_benchmarks | framework_decirculation.py | external-only reference sets, leave-one-out labels, exact permutation null (10,000 draws) |
| 09_benchmarks | 05_round6_extra_stats.py | gene-set level permutation matrix (8 references x 7 strategies, 10,000 draws), direction 3x3 concordance, docking panel z-scores and per-scoring ranks |
| 09_benchmarks | dude_download.py, dude_receptor.py, prep_receptor.py, prep_ligands.py, positive_controls.py | DUD-E prep and positive controls |
| 09_benchmarks | analyze_masld.py, download_tar.py (C:\dockwork\masld) | early seven-sample human MASLD localization |
| 10_manuscript | build_manuscript.py, build_letter_highlights.py, build_target_list.py | Word manuscript, letter, highlights, target list; build_manuscript inserts Figures 1-9 and Supplementary Figures S1-S6 and all five end tables |
| 10_manuscript | plot_send_receive_panels.py | three-panel send/receive figure (numbers pulled from tables, not hard-coded) |
| 11_mr | step1_instruments.py, step2_classify_cis.py, step3_harmonize_mr.py, step4_coloc.py, step5_export_tables.py (cloud mr3) | eQTLGen instruments, cis/trans, MR estimators, ABF colocalization, CSV export |
| 11_mr | og_client.py, mr_methods.py, coloc.py | OpenGWAS client (token from env), pure-Python MR, Wakefield ABF |
| 11_mr | smr_heidi.py, smr_multiplicity.py (cloud mr3) | pure-Python SMR at the top cis SNP with delta-method SE, analytic HEIDI (Zhu 2016), fallback cis instrument, Bonferroni/BH adjustment |
| 12_drg | find_cellxgene_drg.py, get_human_drg_assets.py, download_human_drg.py, explore_drg.py, check_x.py | CELLxGENE discovery, asset metadata, h5ad download, structure, normalization check |
| 12_drg | analyze_drg.py, build_final_tables.py, plot_drg.py, check_source_genes.py | per-cluster expression, annotation/mediator tables, dot plot, ligand-source check |

## Path and environment caveat

The originals were written during the analysis and use absolute paths
(e.g. `C:\dockwork`, the project work root). When run from this repository, update
the path constants at the top of each script or set the documented environment
variables. The curated modules read all paths from environment variables with
sensible in-module defaults: `CHAIHU_MR_WORKDIR` / `CHAIHU_MR_OUTDIR`,
`CHAIHU_DRG_WORKDIR` / `CHAIHU_DRG_OUTDIR` / `CHAIHU_FIGDIR`,
`CHAIHU_SPATIAL_DIR` for the human multi-section and mouse case-control Visium
scripts, `CHAIHU_MASLD_DIR` for the Harmony integration,
`CHAIHU_DOCK_DIR` / `CHAIHU_POCKET_OUT` for the grid-based pocket-volume
script, `CHAIHU_MD_DIR` for the pose-clustering script, and `CHAIHU_MD2_DIR`
for the independent-replicate analysis. The OpenGWAS token
is read from `OPENGWAS_JWT` or `OPENGWAS_JWT_FILE` and is never committed.
Refactoring every older script into a single shared `config.py` is tracked as a
follow-up cleanup; the current curation preserves the exact code that produced
the reported numbers, which is the priority for auditability.
