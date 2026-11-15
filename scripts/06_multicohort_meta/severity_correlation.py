# -*- coding: utf-8 -*-
"""Authoritative gene-vs-severity Spearman correlations, recomputed from source.

This module verifies the severity numbers reported in the manuscript. It reads
the original GEO series matrices and the GSE130970 TPM table, recomputes every
Spearman correlation, and writes a single redistributable table to
``data/processed_small/severity_correlation_verified.csv``.

Inputs (resolved from config / environment, no hard-coded machine paths):
  * GSE130970_series_matrix.txt.gz and GSE130970_TPM.csv.gz (project root)
  * GSE89632_series_matrix.txt.gz (project root)
  * GPL14951 annotation table, located through CHAIHU_GPL_ANNOT or the
    documented default; required only for the GSE89632 probe-to-gene mapping.

Run:
  python scripts/06_multicohort_meta/severity_correlation.py
"""
import gzip
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS  # noqa: E402

PROJECT = PATHS["project_root"]
OUT = PATHS["processed_small"] / "severity_correlation_verified.csv"

# Human Entrez identifiers for the TPM-keyed GSE130970 table.
ENTREZ = {"IL1B": "3553", "PTGS2": "5743", "IL6": "3569", "NGF": "4803"}


def parse_characteristics(sm_path):
    """Return per-sample numeric severity measures, allowing missing values."""
    raw = gzip.open(sm_path, "rt", encoding="utf-8", errors="replace").read().splitlines()
    geo, measures = [], {}
    for line in raw:
        if line.startswith("!Sample_geo_accession"):
            geo = re.findall(r'"(GSM\d+)"', line)
        if line.startswith("!Sample_characteristics_ch1"):
            vals = re.findall(r'"([^"]*)"', line)
            key = vals[0].split(":")[0].strip()
            nums = []
            for v in vals:
                m = re.search(r":\s*([0-9.]+)\s*$", v)
                nums.append(float(m.group(1)) if m else math.nan)
            if len(nums) == len(geo) and not all(math.isnan(a) for a in nums):
                measures[key] = nums
    return geo, measures, raw


def spearman_nan(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = ~(np.isnan(x) | np.isnan(y))
    return spearmanr(x[keep], y[keep])


def load_gpl_probe_map():
    """Illumina probe ID -> gene symbol from the GPL14951 annotation table."""
    gpl = Path(__import__("os").environ.get(
        "CHAIHU_GPL_ANNOT",
        str(PROJECT / "GPL14951_annot.txt")))
    if not gpl.exists():
        raise FileNotFoundError(
            "GPL14951 annotation not found. Set CHAIHU_GPL_ANNOT to the "
            "downloaded platform table (GPL14951 annot).")
    lines = gpl.read_text(encoding="utf-8", errors="replace").splitlines()
    hidx = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    probe_map = {}
    for r in csv.DictReader(lines[hidx:], delimiter="\t"):
        sym = (r.get("Symbol") or "").strip()
        if sym:
            probe_map[r["ID"].strip()] = sym
    return probe_map


def gse89632_gene_vectors(raw, probe_map):
    """Collapse the in-matrix Illumina probes to per-gene sample vectors."""
    start = next(i for i, l in enumerate(raw) if l.startswith("!series_matrix_table_begin"))
    data = []
    for line in raw[start + 1:]:
        if line.startswith("!series_matrix_table_end"):
            break
        data.append(line)
    gene_vec = {}
    for line in data:
        parts = [p.strip().strip('"') for p in line.split("\t")]
        sym = probe_map.get(parts[0])
        if not sym:
            continue
        vals = [float(v) for v in parts[1:]]
        if sym not in gene_vec or sum(vals) > sum(gene_vec[sym]):
            gene_vec[sym] = vals
    return gene_vec


def main():
    rows = []

    def record(cohort, gene, measure, rho, p):
        rows.append({"cohort": cohort, "gene": gene,
                     "severity_measure": measure,
                     "spearman_rho": round(float(rho), 4),
                     "spearman_p": float(p)})

    # GSE130970: TPM expression keyed by Entrez, aligned in sample-table order.
    _, m1, _ = parse_characteristics(PROJECT / "GSE130970_series_matrix.txt.gz")
    with gzip.open(PROJECT / "GSE130970_TPM.csv.gz", "rt",
                   encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.reader(fh)
        next(rdr)
        tpm = {r[0]: r for r in rdr if r and r[0] in ENTREZ.values()}
    for gene, ez in ENTREZ.items():
        x = [float(v) for v in tpm[ez][1:]]
        for key in ("lobular inflammation grade",
                    "nafld activity score", "fibrosis stage"):
            rho, p = spearman_nan(x, m1[key])
            record("GSE130970", gene, key, rho, p)

    # GSE89632: probe-level log2 intensity collapsed to genes.
    _, m2, raw2 = parse_characteristics(PROJECT / "GSE89632_series_matrix.txt.gz")
    gene_vec = gse89632_gene_vectors(raw2, load_gpl_probe_map())
    for gene in ("CALCRL", "TRPV1"):
        for key in ("fibrosis (stage)", "nafld activity score"):
            rho, p = spearman_nan(gene_vec[gene], m2[key])
            record("GSE89632", gene, key, rho, p)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cohort", "gene", "severity_measure",
                                           "spearman_rho", "spearman_p"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
