"""SMR (summary-data MR) and HEIDI heterogeneity test, pure Python.

Implements the SMR test at the top cis eQTL SNP and the analytic HEIDI test of
Zhu et al. (2016, Nat Genet). No R is required. All paths are environment
overridable so the same file runs in the curated repository.
"""
import json
import os
import sys

import numpy as np
from scipy import stats

WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR", os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.environ.get("CHAIHU_MR_OUTDIR", WORKDIR)
sys.path.insert(0, WORKDIR)
from og_client import OpenGWAS  # noqa: E402

OUTCOMES = {
    "stomach_abdominal_pain": "ukb-b-11413",
    "chronic_abdominal_pain_3mo": "ukb-b-19097",
    "upper_abdominal_pain_R101": "ukb-b-6608",
    "abdominal_pain_R104": "ukb-b-6223",
}
HEIDI_P = 0.05          # eQTL p threshold for HEIDI SNPs
HEIDI_MAX = 20          # cap on number of SNPs
PROXY_R2 = 0.9          # exclude top-SNP proxies above this
FLANK = 100_000         # region half-width around gene body


def smr_test(bx0, sx0, by0, sy0):
    """Wald-style SMR at a single instrument with delta-method SE."""
    b = by0 / bx0
    var = (sy0 ** 2) / (bx0 ** 2) + (by0 ** 2 * sx0 ** 2) / (bx0 ** 4)
    se = np.sqrt(var)
    z = b / se
    p = 2 * stats.norm.sf(abs(z))
    return b, se, z, p


def heidi(og, exposure_id, outcome_id, top_rsid, chrom, start, end):
    """Analytic HEIDI. Returns dict (statistic, p, n_snp, reason)."""
    region = f"{chrom}:{start}-{end}"
    e_rows = og.associations(exposure_id, [region], proxies=0)
    o_rows = og.associations(outcome_id, [region], proxies=0)
    e = {r["rsid"]: r for r in e_rows if r.get("beta") is not None and r.get("p") is not None}
    o = {r["rsid"]: r for r in o_rows if r.get("beta") is not None}
    if top_rsid not in e or top_rsid not in o:
        return {"n_snp": 0, "reason": "top SNP missing in region data"}
    # candidate HEIDI SNPs: eQTL p<0.05, present in outcome, not the top SNP
    cand = [rs for rs, r in e.items() if rs != top_rsid and rs in o and r["p"] < HEIDI_P]
    if not cand:
        return {"n_snp": 0, "reason": "no additional eQTL-significant SNPs"}
    # order by eQTL p, cap
    cand.sort(key=lambda rs: e[rs]["p"])
    cand = cand[:HEIDI_MAX]
    # LD among top + candidates
    ld = og._request("POST", "ld/matrix",
                     [("rsid", [top_rsid] + cand), ("pop", "EUR")])
    idx = {s.split("_")[0]: i for i, s in enumerate(ld["snplist"])}
    R = np.array(ld["matrix"], dtype=float)
    # align candidate order to rows that exist in LD
    cand = [rs for rs in cand if rs in idx and top_rsid in idx]
    if not cand:
        return {"n_snp": 0, "reason": "LD reference missing"}
    # drop close proxies of the top SNP
    top_i = idx[top_rsid]
    keep = [rs for rs in cand if R[top_i, idx[rs]] ** 2 < PROXY_R2]
    if len(keep) < 2:
        return {"n_snp": len(keep), "reason": "fewer than 2 non-proxy SNPs"}
    # SMR slope at top
    et, ot = e[top_rsid], o[top_rsid]
    b_smr = ot["beta"] / et["beta"]
    # per-SNP d and covariance
    rows = [idx[rs] for rs in keep]
    m = len(keep)
    d = np.empty(m)
    sx = np.empty(m)
    sy = np.empty(m)
    for j, rs in enumerate(keep):
        d[j] = o[rs]["beta"] - b_smr * e[rs]["beta"]
        sx[j] = e[rs]["se"]
        sy[j] = o[rs]["se"]
    Rsub = R[np.ix_(rows, rows)]
    Dx = np.diag(sx ** 2)
    V = np.diag(sy ** 2) + (b_smr ** 2) * Rsub @ Dx @ Rsub
    # guard against singular V
    V = V + np.eye(m) * 1e-12
    T = float(d @ np.linalg.solve(V, d))
    dof = m - 1
    p = float(stats.chi2.sf(T, dof))
    return {"n_snp": m, "statistic": round(T, 4), "df": dof,
            "p_heidi": round(p, 4), "b_smr": round(b_smr, 5)}


def main():
    og = OpenGWAS()
    cis = json.load(open(os.path.join(WORKDIR, "cis_instruments.json")))
    coords = json.load(open(os.path.join(WORKDIR, "gene_coords_hg19.json")))
    exposure = json.load(open(os.path.join(WORKDIR, "exposure_ids.json"))) \
        if os.path.exists(os.path.join(WORKDIR, "exposure_ids.json")) else None
    # map gene -> eQTL dataset id
    if exposure is None:
        ensg = {
            "IL1B": "ENSG00000125538", "TNF": "ENSG00000232810",
            "PTGS2": "ENSG00000073756", "IL6": "ENSG00000136244",
            "CALCRL": "ENSG00000064989", "TACR1": "ENSG00000115350",
            "TRPV1": "ENSG00000196689",
        }
        exposure = {g: f"eqtl-a-{v}" for g, v in ensg.items()}

    rows, detail = [], {}
    for gene, inst in cis.items():
        # cis instruments ordered by strength (F) descending
        cand = sorted(inst, key=lambda r: r.get("F", (r["beta"] / r["se"]) ** 2),
                      reverse=True)
        top = cand[0]
        c = coords[gene]
        chrom, start, end = c["chr"], c["start"] - FLANK, c["end"] + FLANK
        detail[gene] = {"top_snp": top["rsid"], "outcomes": {}}
        for out_name, out_id in OUTCOMES.items():
            # choose the strongest cis SNP that is present in this outcome
            chosen, ot = None, None
            for snp in cand:
                o_rows = og.associations(out_id, [snp["rsid"]], proxies=0)
                if o_rows and o_rows[0].get("beta") is not None:
                    chosen, ot = snp, o_rows[0]
                    break
            if chosen is None:
                continue
            is_fallback = chosen["rsid"] != top["rsid"]
            b, se, z, p = smr_test(chosen["beta"], chosen["se"], ot["beta"], ot["se"])
            h = heidi(og, exposure[gene], out_id, chosen["rsid"], chrom, start, end)
            if is_fallback:
                note = (h.get("reason", "") +
                        f" fallback instrument (top {top['rsid']} absent)").strip()
                h = dict(h)
                h["reason"] = note
            detail[gene]["outcomes"][out_name] = {
                "smr_b": b, "smr_se": se, "smr_p": p,
                "instrument_snp": chosen["rsid"], "is_fallback": is_fallback,
                "heidi": h}
            rows.append({
                "gene": gene, "outcome": out_name, "instrument_snp": chosen["rsid"],
                "is_fallback": int(is_fallback),
                "smr_b": round(b, 5), "smr_se": round(se, 5), "smr_p": round(p, 4),
                "heidi_n_snp": h.get("n_snp", 0),
                "heidi_stat": h.get("statistic"),
                "heidi_p": h.get("p_heidi"),
                "heidi_note": h.get("reason", ""),
            })
    json.dump(detail, open(os.path.join(OUTDIR, "smr_heidi_results.json"), "w"), indent=1)
    import csv
    out_csv = os.path.join(OUTDIR, "smr_heidi_table.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} SMR/HEIDI rows")
    for r in rows:
        if r["outcome"] == "stomach_abdominal_pain":
            print(f"{r['gene']:7s} SMR b={r['smr_b']:>9} p={r['smr_p']:<7} "
                  f"HEIDI n={r['heidi_n_snp']} p={r['heidi_p']} {r['heidi_note']}")


if __name__ == "__main__":
    main()
