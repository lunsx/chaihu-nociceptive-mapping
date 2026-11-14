# -*- coding: utf-8 -*-
"""Replace TrEMBL hits with their reviewed Swiss-Prot counterpart where one exists.

Reads the 18 unreviewed (TrEMBL) rows from unmapped_复核回捞结果.csv, asks UniProt
for the reviewed human entry of the same gene symbol, and decides:
  REPLACE          reviewed entry exists and its protein name matches the target
  KEEP_TR          no reviewed human entry exists (record stays flagged)
  REJECT_WRONG     the gene symbol does not correspond to the TCMSP target name

Output: uniprot_trembl_canonical.csv / .json
"""
import json
import pathlib
import time

import pandas as pd
import requests

ROOT = pathlib.Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析")
OUT = ROOT / "补强分析_20261114"
API = "https://rest.uniprot.org/uniprotkb/search"
PROXIES = {"http": "http://127.0.0.1:8902", "https": "http://127.0.0.1:8902"}
FIELDS = "accession,id,protein_name,gene_names,reviewed,ec"


def query(q, size=3):
    for attempt in range(3):
        try:
            r = requests.get(API, params=dict(query=q, fields=FIELDS,
                                              format="json", size=size),
                             proxies=PROXIES, timeout=60)
            r.raise_for_status()
            return r.json().get("results", [])
        except Exception as exc:  # noqa: BLE001
            print("retry", attempt, exc)
            time.sleep(3)
    return []


def brief(rec):
    return dict(
        accession=rec.get("primaryAccession"),
        id=rec.get("uniProtkbId"),
        name=(rec.get("proteinDescription", {}).get("recommendedName", {})
              .get("fullName", {}).get("value")),
        genes=[g.get("geneName", {}).get("value")
               for g in rec.get("genes", []) if g.get("geneName")],
        reviewed=rec.get("entryType", ""),
    )


def main():
    res = pd.read_csv(ROOT / "unmapped_复核回捞结果.csv")
    tr = res[res["entry_type"].astype(str).str.contains("TrEMBL", na=False)].copy()

    rows = []
    for r in tr.itertuples(index=False):
        gene = r.human_gene if isinstance(r.human_gene, str) else None
        row = dict(target_name=r.target_name, trembl_acc=r.human_uniprot,
                   trembl_name=r.human_protein_name, gene=gene)
        if not gene:
            row["decision"] = "KEEP_TR"
            row["reason"] = "no human gene symbol on the unreviewed entry"
            rows.append(row)
            continue
        hits = query(f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true")
        time.sleep(0.4)
        if not hits:
            row["decision"] = "KEEP_TR"
            row["reason"] = f"no reviewed human entry for {gene}"
            row["canonical_acc"] = None
            row["canonical_name"] = None
            rows.append(row)
            continue
        b = brief(hits[0])
        row["canonical_acc"] = b["accession"]
        row["canonical_id"] = b["id"]
        row["canonical_name"] = b["name"]
        # crude semantic check: target words vs reviewed protein name
        tgt_words = {w.lower() for w in str(r.target_name).replace("-", " ").split()
                     if len(w) > 3}
        name_words = {w.lower() for w in str(b["name"] or "").replace("-", " ").split()}
        overlap = tgt_words & name_words
        if overlap:
            row["decision"] = "REPLACE"
            row["reason"] = "reviewed entry shares target wording: " + \
                            ", ".join(sorted(overlap))[:80]
        else:
            row["decision"] = "REJECT_WRONG"
            row["reason"] = ("gene symbol leads to an entry whose name does not match "
                             "the TCMSP target name: " + str(b["name"]))
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "uniprot_trembl_canonical.csv", index=False, encoding="utf-8-sig")
    (OUT / "uniprot_trembl_canonical.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")

    print(out.decision.value_counts().to_string())
    print()
    print(out[["target_name", "gene", "trembl_acc", "decision", "canonical_acc",
               "canonical_name"]].to_string(index=False))


if __name__ == "__main__":
    main()
