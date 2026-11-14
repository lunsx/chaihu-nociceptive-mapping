# -*- coding: utf-8 -*-
"""Close the remaining KEEP_TR rows that have a resolvable protein name.

Rows without a gene symbol but with an unambiguous human protein name:
  Adenylsuccinate synthetase        -> ASS1
  Methylmalonyl-CoA mutase          -> MMUT
  Branched chain aminotransferase   -> BCAT1 / BCAT2 (two isoforms)
  Alkaline phosphatase              -> ALPL / ALPP / ALPPL2 / ALPI (four isozymes)
Output: uniprot_name_resolved.csv
"""
import json
import pathlib
import time

import pandas as pd
import requests

OUT = pathlib.Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\补强分析_20261114")
API = "https://rest.uniprot.org/uniprotkb/search"
PROXIES = {"http": "http://127.0.0.1:8902", "https": "http://127.0.0.1:8902"}

QUERIES = {
    "Adenylsuccinate synthetase": "gene_exact:ASS1 AND organism_id:9606 AND reviewed:true",
    "Methylmalonyl-CoA mutase": "gene_exact:MMUT AND organism_id:9606 AND reviewed:true",
    "Branched chain amino acid-dependent aminotransferase": "gene_exact:BCAT1 OR gene_exact:BCAT2 AND organism_id:9606 AND reviewed:true",
    "Alkaline phosphatase": "gene_exact:ALPL OR gene_exact:ALPP OR gene_exact:ALPPL2 OR gene_exact:ALPI AND organism_id:9606 AND reviewed:true",
}


def q(query):
    r = requests.get(API, params=dict(query=query, format="json", size=5,
                                      fields="accession,id,protein_name,gene_names"),
                     proxies=PROXIES, timeout=60)
    r.raise_for_status()
    return r.json().get("results", [])


rows = []
for target, query in QUERIES.items():
    hits = q(query)
    time.sleep(0.5)
    for h in hits:
        rows.append(dict(
            target_name=target,
            accession=h.get("primaryAccession"),
            entry_id=h.get("uniProtkbId"),
            name=(h.get("proteinDescription", {}).get("recommendedName", {})
                  .get("fullName", {}).get("value")),
            gene=(h.get("genes", [{}])[0].get("geneName", {}).get("value")
                  if h.get("genes") else None),
        ))

df = pd.DataFrame(rows)
df.to_csv(OUT / "uniprot_name_resolved.csv", index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
