"""Step 4: cis-region colocalization between each eQTL and abdominal pain GWAS."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from og_client import OpenGWAS
from coloc import coloc_abf

WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(__file__))
OUTCOME = "ukb-b-11413"
FLANK = 100_000
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
GENES = {
    "IL1B": "ENSG00000125538",
    "TNF": "ENSG00000232810",
    "PTGS2": "ENSG00000073756",
    "IL6": "ENSG00000136244",
    "CALCRL": "ENSG00000064989",
    "TACR1": "ENSG00000115350",
    "TRPV1": "ENSG00000196689",
}


def align(exp_rows, out_rows):
    o_by_id = {r["rsid"]: r for r in out_rows}
    b1, v1, b2, v2, n = [], [], [], [], 0
    for e in exp_rows:
        o = o_by_id.get(e["rsid"])
        if o is None or o.get("beta") is None:
            continue
        ea, nea = e["ea"].upper(), e["nea"].upper()
        oea, onea = o["ea"].upper(), o["nea"].upper()
        if ea == oea and nea == onea:
            flip = False
        elif ea == onea and nea == oea:
            flip = True
        else:
            continue
        if frozenset((ea, nea)) in PALINDROMIC:
            if abs(e.get("eaf", 0.5) - o.get("eaf", 0.5)) > 0.2:
                continue
        b1.append(e["beta"])
        v1.append(e["se"] ** 2)
        b2.append(-o["beta"] if flip else o["beta"])
        v2.append(o["se"] ** 2)
        n += 1
    return np.array(b1), np.array(v1), np.array(b2), np.array(v2), n


def main():
    og = OpenGWAS()
    coords = json.load(open(os.path.join(WORKDIR, "gene_coords_hg19.json")))
    summary = {}
    for gene, ensg in GENES.items():
        c = coords[gene]
        chrom, start, end = c["chr"], c["start"], c["end"]
        lo, hi = max(1, start - FLANK), end + FLANK
        rng = f"{chrom}:{lo}-{hi}"
        exp_rows = og.associations(f"eqtl-a-{ensg}", [rng], proxies=0)
        out_rows = og.associations(OUTCOME, [rng], proxies=0)
        b1, v1, b2, v2, n = align(exp_rows, out_rows)
        if n < 10:
            summary[gene] = {"n_overlap": n, "note": "too few overlapping SNPs"}
            print(gene, "too few SNPs:", n)
            continue
        # standard coloc priors: eQTL W=0.04^2, binary log-odds W=0.15^2
        pp = coloc_abf(b1, v1, b2, v2, w1=0.04 ** 2, w2=0.15 ** 2)
        # prior sensitivity on p12
        sens = {}
        for p12 in (1e-6, 5e-6, 1e-5, 2e-5):
            s = coloc_abf(b1, v1, b2, v2, w1=0.04 ** 2, w2=0.15 ** 2, p12=p12)
            sens[str(p12)] = round(s["PP.H4"], 4)
        summary[gene] = {"n_overlap": n, "pp": {k: round(x, 4) for k, x in pp.items()},
                         "PPH4_prior_sensitivity": sens}
        print(f"{gene:7s} n={n:4d} PP.H4={pp['PP.H4']:.4f} "
              f"PP.H3={pp['PP.H3']:.4f} PP.H4 range={min(sens.values()):.3f}-{max(sens.values()):.3f}")

    json.dump(summary, open(os.path.join(WORKDIR, "coloc_results.json"), "w"), indent=1)
    print("\nsaved coloc_results.json")


if __name__ == "__main__":
    main()
