"""Step 5: export clean CSV tables for the manuscript and repository."""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(__file__))
# Default export target is the repository's processed_small data directory.
_DEFAULT_OUT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "data", "processed_small"))
OUTDIR = os.environ.get("CHAIHU_MR_OUTDIR", _DEFAULT_OUT)

mr = json.load(open(os.path.join(WORKDIR, "mr_results.json")))
coloc = json.load(open(os.path.join(WORKDIR, "coloc_results.json")))
cls = json.load(open(os.path.join(WORKDIR, "instrument_classification.json")))

os.makedirs(OUTDIR, exist_ok=True)

with open(os.path.join(OUTDIR, "mr_results_table.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(mr[0].keys()))
    w.writeheader()
    w.writerows(mr)

with open(os.path.join(OUTDIR, "coloc_table.csv"), "w", newline="") as fh:
    fields = ["gene", "n_overlap", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
              "H4_1e-06", "H4_5e-06", "H4_1e-05", "H4_2e-05"]
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    for gene, info in coloc.items():
        if "pp" not in info:
            continue
        row = {"gene": gene, "n_overlap": info["n_overlap"], **info["pp"]}
        row.update({f"H4_{k}": v for k, v in info["PPH4_prior_sensitivity"].items()})
        w.writerow(row)

with open(os.path.join(OUTDIR, "instrument_table.csv"), "w", newline="") as fh:
    fields = ["gene", "rsid", "chr", "pos", "beta", "se", "p", "F", "cis",
              "distance_to_gene", "ea", "nea", "eaf"]
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    for row in cls:
        w.writerow({k: row.get(k) for k in fields})

print("exported mr_results_table.csv, coloc_table.csv, instrument_table.csv to", OUTDIR)
print("MR rows:", len(mr), "| coloc genes:", sum(1 for v in coloc.values() if "pp" in v))
