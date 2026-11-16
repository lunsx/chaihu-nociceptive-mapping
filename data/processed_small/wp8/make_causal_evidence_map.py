# -*- coding: utf-8 -*-
"""Generate the graded causal evidence map (WP8).

The grades are design-based, not p-value based. Each claim receives the highest
tier with positive support; replication at the same tier is recorded as the
number of independent layers but never raises the tier; a null at a higher tier
caps the causal claim; unclosed links get no tier and are never filled by
assumption. Standard library only.
"""
import csv
import collections

OUTDIR = "/home/user/compute_work/wp8"

# Tier definitions: tier, name, study_design, causal_weight, supports, cannot_show
TIERS = [
    ("T0", "Asserted database association",
     "Curated or literature-derived edges imported from secondary databases without an intervention or a direction measured in this study",
     "Minimal; existence of a reported relationship",
     "That a relationship has been asserted in a source",
     "Direction, effect size, or that the relationship occurs in the studied system"),
    ("T1", "Replicated observational co-change",
     "Concordant direction across independent human cohorts, single-cell and spatial localization, and disease-stratified comparison",
     "Low to moderate; correlated change that remains confounded",
     "That a target changes with disease and the cell type or tissue where it is expressed",
     "That the change causes the outcome or that confounding is removed"),
    ("T2", "Curated physical or biochemical interaction",
     "Validated docking with redock and DUD-E, and independently curated ligand-receptor pairs with secretion and receptor annotations",
     "Moderate; a physical mechanism or contact",
     "A plausible physical binding or signaling interaction",
     "An in-system measured effect or receptor activation"),
    ("T3", "Interventional or perturbation evidence",
     "Directed in silico knockout and overexpression with signal propagation, converging with independent LINCS L1000 signatures and high-throughput experimental targets",
     "Highest available in this study; intervention-like and directional",
     "A directed perturbation-response relationship",
     "That a model-based or external-context perturbation reproduces in vivo"),
    ("T4", "Human genetic causal inference",
     "Two-sample MR, Bayesian colocalization, and SMR with HEIDI using germline instruments",
     "Strongest design for a human causal effect; the result here is null",
     "A germline causal effect when positive; here no positive support",
     "Tissue-specific or local non-germline effects; a null is not evidence of absence"),
]

# Claim map: id, claim, source, target, highest_positive_tier, independent_layers,
# direction, status, reference, upgrade_experiment
CLAIMS = [
    ("C1", "Catalog compounds are linked to targets through source databases",
     "Catalog compounds", "Mapped targets", "T0",
     "TCMSP; HERB aggregation of TCMSP, SymMap, TCMID, and HIT",
     "No direction", "Positive existence", "Methods 2.3",
     "Direct binding assay"),
    ("C2", "Retained compounds physically engage inflammatory targets",
     "Retained compounds", "Inflammatory targets", "T2",
     "Docking ensemble; redock; DUD-E benchmark",
     "Binding without a sign", "Positive", "Methods 2.8; Fig 7",
     "Binding assay and co-crystal structure"),
    ("C3", "Inflammatory cells produce and upregulate the mediators in disease",
     "Macrophages, monocytes, endothelial cells", "IL1B, TNF, PTGS2, PGE2", "T1",
     "Single-cell; spatial; multi-cohort differential expression",
     "Up in disease", "Positive", "Methods 2.4 and 2.7",
     "Protein-level measurement in tissue"),
    ("C4", "Upstream transcription factors drive the inflammatory module",
     "FOS, RELA, NFKB1", "IL1B, TNF, PTGS2, NFKBIA", "T3",
     "CellOracle knockout and overexpression; LINCS L1000",
     "Knockout reduces; overexpression increases", "Positive, model-based",
     "Methods 2.10; Fig 9", "Physical CRISPR perturbation"),
    ("C5", "High-throughput experiments confirm herb target engagement",
     "HERB experiment channel", "Experimental targets including pain genes", "T3",
     "High-throughput assays with p and FDR; 209 targets for Chaihu and 134 for Jinling",
     "Empirical", "Positive, external context", "Methods 2.19",
     "Confirmatory targeted assay"),
    ("C6", "Mediators bind their receptors",
     "PGE2, IL6, TNF, NGF", "PTGER3 and PTGER4, IL6ST, TNFRSF1A, NTRK1", "T2",
     "Curated ligand-receptor; GPCR pharmacology",
     "Agonism", "Positive", "Methods 2.12",
     "Functional receptor assay"),
    ("C7", "Mediators are secreted",
     "Inflammatory cells", "Mediators", "T2",
     "Curated secretion annotation plus expression",
     "Release", "Positive in literature", "Methods 2.12",
     "Secretion assay"),
    ("C8", "Human nociceptors express the matching receptors",
     "Human DRG neurons", "PTGER3, NTRK1, IL6ST, OSMR", "T1",
     "1,837 nuclei across six donors",
     "Presence", "Positive", "Methods 2.7",
     "Functional receptor assay"),
    ("C9", "Mediators are transported to the ganglion",
     "Liver", "Ganglion", "No tier",
     "None",
     "Unknown", "Unclosed", "Methods 2.12; Table 6",
     "In vivo pharmacokinetics and nerve concentration"),
    ("C10", "Mediators reach the ganglion at an active concentration",
     "Circulating mediators", "Ganglion", "No tier",
     "None",
     "Unknown", "Unclosed", "Methods 2.12; Table 6",
     "Ganglion microdialysis"),
    ("C11", "Receptors activate and nociceptors fire",
     "Mediators", "Nociceptor firing", "No tier",
     "None",
     "Unknown", "Unclosed", "Methods 2.12; Table 6",
     "Electrophysiology"),
    ("C12", "Genetically determined expression causes abdominal pain",
     "Germline expression", "Abdominal pain outcomes", "T4",
     "MR; colocalization; SMR with HEIDI",
     "No supported direction", "Null", "Methods 2.11; Fig 9",
     "Liver eQTL and a larger abdominal pain GWAS"),
    ("C13", "The core six-gene identity is retained across data sources",
     "TCMSP and HERB", "Core six genes", "T1",
     "Cross-source identity; high-throughput channel",
     "Identity retained; enrichment is source-dependent",
     "Positive identity", "Methods 2.19",
     "An independent primary-source graph"),
]


def write_csv(name, header, rows):
    path = f"{OUTDIR}/{name}"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print("wrote", name, len(rows), "rows")


def main():
    write_csv(
        "causal_evidence_tiers.csv",
        ["tier", "name", "study_design", "causal_weight", "supports", "cannot_show"],
        TIERS,
    )
    write_csv(
        "causal_evidence_map.csv",
        ["claim_id", "mechanistic_claim", "source_node", "target_node",
         "highest_positive_tier", "independent_layers", "direction", "status",
         "reference", "upgrade_experiment"],
        CLAIMS,
    )

    counts = collections.Counter()
    for c in CLAIMS:
        counts[c[4]] += 1
    summary = [(tier, name, counts.get(tier, 0)) for tier, name, *_ in TIERS]
    summary.append(("No tier", "Unclosed links", counts.get("No tier", 0)))
    write_csv(
        "causal_grade_summary.csv",
        ["tier", "name", "n_claims"],
        summary,
    )

    print("\nTier counts:")
    for tier, name, n in summary:
        print(f"  {tier:7s} {n} claims")


if __name__ == "__main__":
    main()
