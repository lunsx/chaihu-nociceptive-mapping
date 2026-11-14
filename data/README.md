# Data directory

## processed_small/

Small derived tables and summary JSON files that are safe to version in git.
They let a reviewer inspect the headline numbers without rerunning the pipeline.

| File | Module | Contents |
|---|---|---|
| dilution_stats.json | 09 | dilution correlations, median expression percentage, attenuation |
| dilution_gene_table.csv | 09 | per-gene dilution statistics (54 genes) |
| meta_ci_table.csv | 06/09 | Hedges' g, 95% CI, TOST outcomes |
| framework_benchmark.csv | 09 | precision/recall/OR across selection stages |
| enrichment_summary.json | 04 | DUD-E ROC AUC, BEDROC, enrichment factors |
| controls_luteolin.json | 04 | positive-control and luteolin docking energies |
| docking_results.csv | 04 | DUD-E benchmark ligand docking scores |
| masld_summary.json | 09 | human MASLD localization summary |
| Table1..Table5 | manuscript | the five submission tables |

## raw_external_links/

Large raw files are not stored in git. `manifest.csv` records each required input,
its public URL, expected size, and where the module expects it. Download the raw
data following the manifest, then run the matching module.

## Regenerating the processed tables

Each curated script under `scripts/` reads the raw data and writes the file of the
same name. The mapping between script, raw input, and output is in
`docs/code_map.md`.
