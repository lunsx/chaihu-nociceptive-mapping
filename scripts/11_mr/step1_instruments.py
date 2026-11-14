"""Step 1: confirm all 7 eQTLGen exposure datasets and pull clumped top hits."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from og_client import OpenGWAS

WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(__file__))

# Core genes with Ensembl IDs (no version suffix, as OpenGWAS indexes them)
GENES = {
    "IL1B": "ENSG00000125538",
    "TNF": "ENSG00000232810",
    "PTGS2": "ENSG00000073756",
    "IL6": "ENSG00000136244",
    "CALCRL": "ENSG00000064989",
    "TACR1": "ENSG00000115350",
    "TRPV1": "ENSG00000196689",
}


def main():
    og = OpenGWAS()
    ids = [f"eqtl-a-{e}" for e in GENES.values()]
    info = og.gwasinfo(ids)
    info_by_id = {rec["id"]: rec for rec in info}
    print("=== gwasinfo availability ===")
    for sym, ensg in GENES.items():
        gid = f"eqtl-a-{ensg}"
        rec = info_by_id.get(gid)
        if rec:
            print(f"{sym:7s} {gid:28s} N={rec['sample_size']:>6} nsnp={rec['nsnp']}")
        else:
            print(f"{sym:7s} {gid:28s} MISSING")

    print("\n=== clumped top hits (p <= 5e-8, r2=0.001, kb=5000) ===")
    instruments = {}
    for sym, ensg in GENES.items():
        gid = f"eqtl-a-{ensg}"
        hits = og.tophits(gid, pval=5e-8, clump=1, preclumped=0)
        rows = []
        for h in hits:
            beta = h.get("beta")
            se = h.get("se")
            p = h.get("p")
            f = (beta / se) ** 2 if beta is not None and se else None
            rows.append({
                "rsid": h.get("rsid") or h.get("target_snp"),
                "chr": h.get("chr"), "pos": h.get("position"),
                "beta": beta, "se": se, "p": p,
                "F": round(f, 2) if f is not None else None,
                "ea": h.get("ea"), "nea": h.get("nea"), "eaf": h.get("eaf"),
            })
        instruments[sym] = rows
        n_pass = sum(1 for r in rows if r["F"] is not None and r["F"] >= 10)
        print(f"{sym:7s} n_snps={len(rows):3d} F>=10:{n_pass:3d}",
              "| lead:", rows[0]["rsid"] if rows else "none",
              f"F={rows[0]['F']}" if rows else "")

    out = os.path.join(WORKDIR, "exposure_instruments.json")
    with open(out, "w") as fh:
        json.dump(instruments, fh, indent=1)
    print("\nsaved exposure_instruments.json")


if __name__ == "__main__":
    main()
