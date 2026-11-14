# -*- coding: utf-8 -*-
"""Write back the round-6 TrEMBL adjudication into the master mapping.

Actions:
  1. back up target_gene_mapping.csv
  2. map the two REPLACE targets to reviewed human genes
  3. emit the final 18-row accession-level verdict table
  4. verify mapping counts (693 mapped, 253 unmapped, 664 unique genes)
"""
import csv
import shutil
from datetime import datetime
from pathlib import Path

ROOT = Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析")
MAP = ROOT / "target_gene_mapping.csv"
ADJ = ROOT / "补强分析_20261114" / "uniprot_trembl_adjudicated.csv"
FINAL = ROOT / "uniprot_TrEMBL_final_verdict.csv"

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
backup = MAP.with_name(f"target_gene_mapping_备份_{stamp}.csv")
shutil.copy2(MAP, backup)
print("backup:", backup.name)

REPLACE = {
    "6-phosphofructokinase": ("PFKP", "Q01813"),
    "Orotidine 5'-phosphate decarboxylase": ("UMPS", "P11172"),
}

rows = list(csv.reader(MAP.open(encoding="utf-8")))
header, data = rows[0], rows[1:]
changed = []
for r in data:
    if r[0] in REPLACE and r[2] == "unmapped":
        gene, acc = REPLACE[r[0]]
        r[1], r[2] = gene, "replaced_reviewed_r2"
        changed.append((r[0], gene, acc))
with MAP.open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    w.writerows(data)
print("changed rows:", changed)

# final verdict table: merge adjudication with the writeback actions
verdict_rows = list(csv.DictReader(ADJ.open(encoding="utf-8-sig")))
out_header = ["target_name", "original_trembl_acc", "verdict", "final_gene",
              "final_uniprot", "action", "note"]
action_map = {
    "6-phosphofructokinase": ("PFKP", "Q01813", "mapped to reviewed human entry"),
    "Orotidine 5'-phosphate decarboxylase": ("UMPS", "P11172", "mapped to reviewed human entry"),
}
with FINAL.open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(out_header)
    for r in verdict_rows:
        gene, acc, action = "", "", "stays unmapped; see verdict"
        if r["target_name"] in action_map:
            gene, acc, action = action_map[r["target_name"]]
        elif r["verdict"] == "REJECT_WRONG":
            action = "rejected; not a valid human target"
        elif r["verdict"] == "AMBIGUOUS":
            action = "ambiguous; cannot map uniquely"
        elif r["verdict"] == "KEEP_TR":
            action = "no reviewed human entry; kept as TrEMBL reference"
        note = r["note"]
        if r["target_name"] == "Adenylsuccinate synthetase":
            note = ("ASS1 P00966 is argininosuccinate synthase; the correct reviewed "
                    "entries are ADSS1 Q8N142 / ADSS2 P30520")
        w.writerow([r["target_name"], r["trembl_acc"], r["verdict"],
                    gene, acc, action, note])
print("final verdict table:", FINAL.name, "rows:", len(verdict_rows))

# verify counts
mapped = [r for r in data if r[2] != "unmapped" and r[1]]
unmapped = [r for r in data if r[2] == "unmapped"]
genes = {r[1] for r in data if r[1]}
print(f"total rows: {len(data)} | mapped: {len(mapped)} | unmapped: {len(unmapped)} | unique genes: {len(genes)}")
assert len(data) == 946 and len(mapped) == 693 and len(unmapped) == 253 and len(genes) == 664
