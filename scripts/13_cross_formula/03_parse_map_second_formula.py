# -*- coding: utf-8 -*-
"""Parse the second-formula TCMSP pages and map targets to human genes.

Reads the two raw HTML pages produced by 02_fetch_second_formula.py, extracts
the embedded kendoGrid JSON arrays, applies the same OB >= 30 and DL >= 0.18
filter used for the first formula, and resolves target names to genes. Target
names already resolved for the first formula are reused, and new names are
matched against the reviewed UniProt table; unmapped rows are reported.

Outputs (project work root, small gene mapping also copied to processed_small):
  jinling_compound_master.csv
  jinling_compound_targets.csv
  jinling_gene_mapping.csv
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

PROJECT = PATHS["project_root"]
RAW_DIR = PROJECT / "xiaoyao_case" / "jinling_raw"
WORK = PROJECT / "xiaoyao_case"
PAGES = [("Toosendan Fructus", "Toosendan_Fructus.html"),
         ("Corydalis Rhizoma", "Corydalis_Rhizoma.html")]

mol_rows, targ_rows = [], []
for herb, file_name in PAGES:
    html = (RAW_DIR / file_name).read_text(encoding="utf-8")
    arrays = [json.loads(m.group(1)) for m in
              re.finditer(r"data:\s*(\[.*?\])\s*,\s*pageSize", html, re.S)]
    molecules, targets = arrays[0], arrays[1]
    for row in molecules:
        row = dict(row); row["herb"] = herb; mol_rows.append(row)
    for row in targets:
        row = dict(row); row["herb"] = herb; targ_rows.append(row)

M = pd.DataFrame(mol_rows)
T = pd.DataFrame(targ_rows)
for col in ["ob", "dl", "mw", "bbb", "caco2"]:
    M[col] = pd.to_numeric(M[col], errors="coerce")
M["pass_filter"] = (M["ob"] >= 30) & (M["dl"] >= 0.18)


def join_herbs(values):
    return ",".join(sorted(set(values)))


master = (M.groupby(["MOL_ID", "molecule_name"], as_index=False)
          .agg(OB=("ob", "first"), DL=("dl", "first"),
               herbs=("herb", join_herbs), n_herbs=("herb", "nunique"),
               pass_filter=("pass_filter", "max")))
tmaster = (T.groupby(["MOL_ID", "molecule_name", "target_name", "target_ID"], as_index=False)
           .agg(herbs=("herb", join_herbs), validated=("validated", "max")))
master.to_csv(WORK / "jinling_compound_master.csv", index=False, encoding="utf-8-sig")
tmaster.to_csv(WORK / "jinling_compound_targets.csv", index=False, encoding="utf-8-sig")

# gene mapping: reuse the first-formula map, then resolve new names from UniProt
auth = pd.read_csv(PROJECT / "target_gene_mapping.csv")
auth_map = {r.target_name: (r.gene if isinstance(r.gene, str) else "") for _, r in auth.iterrows()}
uniprot = pd.read_csv(PROJECT / "uniprot_human_reviewed.tsv", sep="\t", dtype=str)


def norm(text):
    text = text.lower()
    text = re.sub(r"\(ec [\d\.\- ]+\)", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text).strip()
    return re.sub(r"\s+", " ", text)


def recommended(text):
    return re.sub(r"\(ec [\d\.\- ]+\)", " ", text).split("(")[0].strip()


exact, rec = {}, {}
for _, r in uniprot.iterrows():
    gene = r["Gene Names (primary)"]
    if not isinstance(gene, str) or not gene.strip():
        continue
    names = r["Protein names"]
    exact.setdefault(norm(names), gene)
    rec.setdefault(norm(recommended(names)), gene)
    for alt in re.findall(r"\(([^()]+)\)", names):
        if not re.search(r"ec ", alt.lower()):
            rec.setdefault(norm(alt), gene)


def resolve_new(name):
    n = norm(name)
    return exact.get(n) or rec.get(n) or rec.get(norm(recommended(name)))


pass_mol = set(master[master["pass_filter"]]["MOL_ID"])
filtered_targets = sorted(tmaster[tmaster["MOL_ID"].isin(pass_mol)]["target_name"].dropna().unique())
mapping_rows, unmapped = [], []
for name in filtered_targets:
    gene = auth_map.get(name) or resolve_new(name)
    if gene:
        mapping_rows.append((name, gene))
    else:
        unmapped.append(name)

mapping = pd.DataFrame(mapping_rows, columns=["target_name", "gene"])
mapping.to_csv(WORK / "jinling_gene_mapping.csv", index=False, encoding="utf-8-sig")
mapping.to_csv(PATHS["processed_small"] / "jinling_gene_mapping.csv", index=False,
               encoding="utf-8-sig")

print("unique molecules:", len(master), "| pass filter:", int(master["pass_filter"].sum()))
print("filtered genes:", mapping["gene"].nunique())
print("unmapped filtered targets:", len(unmapped))
