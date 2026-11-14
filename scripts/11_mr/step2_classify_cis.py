"""Step 2: classify eQTL instruments as cis (same chr, within +/- 1 Mb) or trans,
using GRCh37 gene coordinates from the Ensembl grch37 REST service."""
import json
import os
import sys
import time

import requests

sys.path.insert(0, os.path.dirname(__file__))
WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(__file__))

GENES = {
    "IL1B": "ENSG00000125538",
    "TNF": "ENSG00000232810",
    "PTGS2": "ENSG00000073756",
    "IL6": "ENSG00000136244",
    "CALCRL": "ENSG00000064989",
    "TACR1": "ENSG00000115350",
    "TRPV1": "ENSG00000196689",
}
WINDOW = 1_000_000
SESSION = requests.Session()
SESSION.trust_env = False


def gene_coords(sym):
    url = (f"https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}"
           "?content-type=application/json")
    for attempt in range(5):
        try:
            r = SESSION.get(url, timeout=30)
            r.raise_for_status()
            j = r.json()
            return int(j["seq_region_name"]), int(j["start"]), int(j["end"])
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2)


def main():
    instruments = json.load(open(os.path.join(WORKDIR, "exposure_instruments.json")))
    coords = {}
    cis_set = {}
    table = []
    for sym in GENES:
        chrom, start, end = gene_coords(sym)
        coords[sym] = {"chr": chrom, "start": start, "end": end}
        cis_rows = []
        for r in instruments[sym]:
            same = int(r["chr"]) == chrom
            dist = 0
            if r["pos"] < start:
                dist = start - r["pos"]
            elif r["pos"] > end:
                dist = r["pos"] - end
            is_cis = same and dist <= WINDOW
            r2 = dict(r)
            r2["cis"] = is_cis
            r2["distance_to_gene"] = dist
            table.append({"gene": sym, **r2})
            if is_cis:
                cis_rows.append(r)
        cis_set[sym] = cis_rows
        print(f"{sym:7s} gene chr{chrom}:{start}-{end} | cis={len(cis_rows)} "
              f"total={len(instruments[sym])}")
        time.sleep(0.4)

    json.dump(coords, open(os.path.join(WORKDIR, "gene_coords_hg19.json"), "w"), indent=1)
    json.dump(cis_set, open(os.path.join(WORKDIR, "cis_instruments.json"), "w"), indent=1)
    json.dump(table, open(os.path.join(WORKDIR, "instrument_classification.json"), "w"), indent=1)
    print("\nsaved cis_instruments.json")
    for sym, rows in cis_set.items():
        mean_f = sum(r["F"] for r in rows) / len(rows)
        print(f"{sym:7s} cis n={len(rows)} SNPs={[r['rsid'] for r in rows]} "
              f"meanF={mean_f:.1f}")


if __name__ == "__main__":
    main()
