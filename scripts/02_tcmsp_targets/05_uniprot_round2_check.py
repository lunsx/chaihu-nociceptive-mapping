# -*- coding: utf-8 -*-
"""Check the two UNRESOLVED target mappings and the TrEMBL entries against UniProt.

- NO66 / URLC2 : is there a reviewed human entry at all?
- DHRS4        : which human entry is canonical?
- TrEMBL rows  : for any target mapped to an unreviewed (TrEMBL) entry that has a
                 gene symbol, ask UniProt for the reviewed entry of the same gene.

Output: uniprot_round2_check.json + uniprot_round2_check.csv
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


def query(q, fields="accession,id,protein_name,gene_names,reviewed,ec,organism_name",
          size=5):
    for attempt in range(3):
        try:
            r = requests.get(API, params={"query": q, "fields": fields,
                                          "format": "json", "size": size},
                             proxies=PROXIES, timeout=60)
            r.raise_for_status()
            return r.json().get("results", [])
        except Exception as exc:  # noqa: BLE001
            print("retry", attempt, exc)
            time.sleep(3)
    return []


def brief(rec):
    return dict(accession=rec.get("primaryAccession"),
                id=rec.get("uniProtkbId"),
                name=rec.get("proteinDescription", {}).get("recommendedName", {})
                .get("fullName", {}).get("value"),
                genes=[g.get("geneName", {}).get("value")
                       for g in rec.get("genes", []) if g.get("geneName")],
                reviewed=rec.get("entryType", ""),
                organism=rec.get("organism", {}).get("scientificName"))


def main():
    res = {}
    res["NO66_URLC2"] = [brief(x) for x in
                         query('gene:URLC2 AND organism_id:9606')]
    res["NO66_byname"] = [brief(x) for x in
                          query('protein_name:"uricalcalin-2" AND organism_id:9606')]
    res["DHRS4"] = [brief(x) for x in query('gene:DHRS4 AND organism_id:9606')]

    # TrEMBL rows with a gene symbol
    mp = pd.read_csv(ROOT / "target_gene_mapping.csv")
    col = "entry_type" if "entry_type" in mp.columns else None
    trembl = mp[mp[col].astype(str).str.contains("TrEMBL", case=False, na=False)] \
        if col else mp.iloc[0:0]
    genes = sorted({str(g) for g in trembl.get("gene", pd.Series(dtype=str)).dropna()})
    print("TrEMBL rows:", len(trembl), "genes:", genes)

    canon = {}
    for g in genes:
        hits = query(f'gene_exact:{g} AND organism_id:9606 AND reviewed:true')
        canon[g] = [brief(x) for x in hits]
        time.sleep(0.4)
    res["trembl_canonical"] = canon

    (OUT / "uniprot_round2_check.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")

    rows = []
    for key, hits in res.items():
        if key == "trembl_canonical":
            for g, hh in hits.items():
                for h in hh[:2]:
                    rows.append(dict(question=f"TrEMBL->canonical:{g}", **h))
        else:
            for h in hits[:5]:
                rows.append(dict(question=key, **h))
    pd.DataFrame(rows).to_csv(OUT / "uniprot_round2_check.csv", index=False,
                              encoding="utf-8-sig")
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: len(vv) for kk, vv in v.items()})
                      for k, v in res.items()}, ensure_ascii=False, indent=1)[:1500])


if __name__ == "__main__":
    main()
