"""Step 3: harmonize cis instruments with each pain outcome and run MR."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from og_client import OpenGWAS
from mr_methods import run_all

WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(__file__))

OUTCOMES = {
    "stomach_abdominal_pain": "ukb-b-11413",
    "chronic_abdominal_pain_3mo": "ukb-b-19097",
    "upper_abdominal_pain_R101": "ukb-b-6608",
    "abdominal_pain_R104": "ukb-b-6223",
}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}


def harmonize(exposure, outcome):
    by_id = {r["rsid"]: r for r in outcome}
    bx, by, sx, sy, kept, dropped = [], [], [], [], [], []
    for e in exposure:
        o = by_id.get(e["rsid"])
        if o is None or o.get("beta") is None:
            dropped.append((e["rsid"], "missing in outcome"))
            continue
        ea, nea = e["ea"].upper(), e["nea"].upper()
        oea, onea = o["ea"].upper(), o["nea"].upper()
        flip = False
        if ea == oea and nea == onea:
            flip = False
        elif ea == onea and nea == oea:
            flip = True
        else:
            dropped.append((e["rsid"], "allele mismatch"))
            continue
        if frozenset((ea, nea)) in PALINDROMIC:
            # ambiguous strand: rely on EAF agreement, otherwise exclude
            if abs(e["eaf"] - o["eaf"]) > 0.2 and abs(e["eaf"] - (1 - o["eaf"])) > 0.2:
                dropped.append((e["rsid"], "ambiguous palindrome"))
                continue
        beta_y = -o["beta"] if flip else o["beta"]
        bx.append(e["beta"])
        by.append(beta_y)
        sx.append(e["se"])
        sy.append(o["se"])
        kept.append(e["rsid"])
    return (np.array(bx), np.array(by), np.array(sx), np.array(sy),
            kept, dropped)


def main():
    og = OpenGWAS()
    cis = json.load(open(os.path.join(WORKDIR, "cis_instruments.json")))
    all_inst = json.load(open(os.path.join(WORKDIR, "exposure_instruments.json")))

    results, detail = [], {}
    for out_name, out_id in OUTCOMES.items():
        detail[out_name] = {}
        for gene, exp_rows in cis.items():
            snps = [r["rsid"] for r in exp_rows]
            out_rows = og.associations(out_id, snps, proxies=0, align_alleles=1)
            bx, by, sx, sy, kept, dropped = harmonize(exp_rows, out_rows)
            if len(kept) == 0:
                detail[out_name][gene] = {"kept": [], "dropped": dropped}
                continue
            estimates = run_all(bx, by, sx, sy)
            detail[out_name][gene] = {
                "kept": kept, "dropped": dropped, "n_snp": len(kept),
                "estimates": {k: list(v) for k, v in estimates.items()},
            }
            for method, vals in estimates.items():
                if method in ("Cochran Q", "Egger intercept"):
                    continue
                theta, se, p = vals
                results.append({
                    "outcome": out_name, "gene": gene, "instrument": "cis",
                    "method": method, "n_snp": len(kept),
                    "beta": round(theta, 5), "se": round(se, 5), "p": round(p, 4),
                })

    # sensitivity: all cis+trans instruments on the primary outcome
    out_id = OUTCOMES["stomach_abdominal_pain"]
    for gene, exp_rows in all_inst.items():
        snps = [r["rsid"] for r in exp_rows]
        out_rows = og.associations(out_id, snps, proxies=0, align_alleles=1)
        bx, by, sx, sy, kept, dropped = harmonize(exp_rows, out_rows)
        if len(kept) == 0:
            continue
        estimates = run_all(bx, by, sx, sy)
        detail["sensitivity_all_instruments"] = detail.get("sensitivity_all_instruments", {})
        detail["sensitivity_all_instruments"][gene] = {
            "kept": kept, "dropped": dropped, "n_snp": len(kept),
            "estimates": {k: list(v) for k, v in estimates.items()},
        }
        for method, vals in estimates.items():
            if method in ("Cochran Q", "Egger intercept"):
                continue
            theta, se, p = vals
            results.append({
                "outcome": "stomach_abdominal_pain", "gene": gene,
                "instrument": "cis+trans", "method": method, "n_snp": len(kept),
                "beta": round(theta, 5), "se": round(se, 5), "p": round(p, 4),
            })

    json.dump(detail, open(os.path.join(WORKDIR, "mr_detail.json"), "w"), indent=1)
    json.dump(results, open(os.path.join(WORKDIR, "mr_results.json"), "w"), indent=1)
    print("=== cis MR, primary outcome (abdominal pain) ===")
    for r in results:
        if r["outcome"] == "stomach_abdominal_pain" and r["instrument"] == "cis":
            print(f"{r['gene']:7s} {r['method']:14s} n={r['n_snp']} "
                  f"beta={r['beta']:>9} p={r['p']}")
    print("total result rows:", len(results))


if __name__ == "__main__":
    main()
