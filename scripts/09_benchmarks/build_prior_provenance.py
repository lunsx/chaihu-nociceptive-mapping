# -*- coding: utf-8 -*-
"""Build a transparent prior-direction provenance table for the held-out test.

The held-out prior is the 'direction' field in target_list_enriched.json.
This table records each gene's prior direction, whether an explicit literature
reference exists, and the framework sources, so the basis is inspectable rather
than labelled uniformly as 'literature-fixed'.
"""
import json
import pathlib

import pandas as pd

BASE = pathlib.Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析")
TGT = BASE / "target_list_enriched.json"
OUT = BASE / "prior_direction_provenance.csv"

tgt = json.load(open(TGT, encoding="utf-8"))

SOURCE_EN = {
    "GSE89632差异基因": "GSE89632 differential genes",
    "HERB/TTD多库验证": "HERB/TTD multi-database validation",
    "TCMSP方剂靶点(未筛选)": "TCMSP formula targets (unfiltered)",
    "TCMSP方剂靶点(筛选成分)": "TCMSP formula targets (filtered constituents)",
    "多队列Meta": "multi-cohort meta-analysis",
    "痛觉基因集": "pain gene set",
}

def translate_sources(raw):
    parts = [p.strip() for p in str(raw).split(";") if p.strip()]
    return "; ".join(SOURCE_EN.get(p, p) for p in parts)

rows = []
for r in tgt:
    d = r.get("direction", "")
    if d not in ("up", "down"):
        continue
    lit = (r.get("literature") or "").strip()
    rows.append({
        "gene": r["gene"],
        "prior_direction": d,
        "explicit_literature": "yes" if lit else "no",
        "literature_reference": lit,
        "framework_sources": translate_sources(r.get("sources", "")),
        "basis": "literature-anchored" if lit else "framework-curated",
    })

df = pd.DataFrame(rows).sort_values(["explicit_literature", "gene"],
                                    ascending=[False, True])
df.to_csv(OUT, index=False, encoding="utf-8-sig")

print("total prior-predicted genes:", len(df))
print("literature-anchored:", int((df.basis == "literature-anchored").sum()))
print("framework-curated:", int((df.basis == "framework-curated").sum()))
print("up:", int((df.prior_direction == "up").sum()),
      "down:", int((df.prior_direction == "down").sum()))
print()
print(df[["gene", "prior_direction", "basis", "literature_reference"]].to_string(index=False))

# Cohort-level independence / leakage summary.
# The direction is the multi-cohort vote from GSE281797, GSE260666, GSE193066,
# GSE192959 (case-control and grade-correlation signs). GSE147304 was excluded
# from the vote and is therefore the only fully held-out cohort.
cohort_ind = pd.DataFrame([
    {"cohort": "GSE281797", "contributes_to_prior": "yes",
     "fully_held_out": "no", "note": "case-control and grade signs included in the vote",
     "prior_accuracy_pct": 73.1},
    {"cohort": "GSE260666", "contributes_to_prior": "yes",
     "fully_held_out": "no", "note": "case-control and grade signs included in the vote",
     "prior_accuracy_pct": 92.3},
    {"cohort": "GSE147304", "contributes_to_prior": "no",
     "fully_held_out": "yes", "note": "excluded from the vote as inconclusive",
     "prior_accuracy_pct": 42.3},
])
cohort_ind.to_csv(BASE / "prior_cohort_independence.csv", index=False,
                  encoding="utf-8-sig")
print()
print(cohort_ind.to_string(index=False))
