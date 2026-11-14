# -*- coding: utf-8 -*-
"""Manual adjudication of the 18 TrEMBL (unreviewed) target mappings.

The automatic name-similarity check in uniprot_trembl_canonical.py is too crude:
it misses synonyms (PFKP, UMPS) and it accepts false friends (IRP1, MAOA, CYTB).
Every row below was read against its UniProt entry and given a verdict:

  REPLACE        a reviewed human entry is the same protein -> use that accession
  KEEP_TR        no reviewed human entry exists; keep the unreviewed record, flagged
  REJECT_WRONG   the gene symbol resolves to a different protein than the target
  AMBIGUOUS      two reviewed paralogues fit the generic target name; needs a call

Output: uniprot_trembl_adjudicated.csv / .json  (supersedes the automatic table)
"""
import json
import pathlib

import pandas as pd

OUT = pathlib.Path(r"C:\Users\Administrator\AppData\Local\Temp") / "uniprot_adj"
TARGET = pathlib.Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\补强分析_20261114")

ROWS = [
    # target_name, trembl_acc, gene, verdict, canonical_acc, canonical_name, note
    ("2-isopropylmalate synthase", "Q9HBB2", "IRP1", "REJECT_WRONG", None, None,
     "IRP1/ACO1 是顺乌头酸酶与铁反应元件结合蛋白，不是异丙基苹果酸合酶；基因符号撞名"),
    ("6-phosphofructokinase", "A0A7I2V3Z0", "PFKP", "REPLACE", "Q01813",
     "ATP-dependent 6-phosphofructokinase, platelet type",
     "PFKP 即血小板型磷酸果糖激酶，同酶不同同工酶，名称同义"),
    ("Adenylsuccinate synthetase", "B3KTV4", None, "KEEP_TR", None, None,
     "无基因符号；片段型 TrEMBL 记录，人类 reviewed 对应为 ASS1（P00966）需按名称另核"),
    ("Alkaline phosphatase", "B7Z1D1", None, "KEEP_TR", None, None,
     "泛名，人类有 ALPL/ALPP/ALPPL2/ALPI 多个 reviewed 同工酶，需指定组织来源才能定"),
    ("Branched chain amino acid-dependent aminotransferase", "B3KSI3", None,
     "KEEP_TR", None, None, "片段型记录；人类 BCAT1/BCAT2 为 reviewed，需按名称另核"),
    ("Cytochrome c4", "Q2HJZ5", "CYTB", "REJECT_WRONG", None, None,
     "细胞色素 c4 与线粒体细胞色素 b（MT-CYB）不是同一蛋白，仅同属细胞色素家族"),
    ("Endoglucanase F", "Q8N909", None, "REJECT_WRONG", None, None,
     "内切葡聚糖酶为微生物酶，人类无对应；Q8N909 注释可疑，判为非人源靶点"),
    ("Gag-Pol polyprotein", "V9H0H3", None, "REJECT_WRONG", None, None,
     "病毒多聚蛋白，非人类靶点"),
    ("L-fucose isomerase", "A0ACI8ST25", "FUOM", "AMBIGUOUS", "A2VDF0",
     "Fucose mutarotase",
     "FUOM 为岩藻糖变旋酶，与 L-岩藻糖异构酶（EC 5.3.1.25）活性不同，历史注释有混用"),
    ("Lysine-specific demethylase NO66", "A0A0A8K8C0", "URLC2", "KEEP_TR", None, None,
     "人类无 reviewed 记录，仅此条 TrEMBL 且名为双功能去甲基化酶/组氨酸羟化酶"),
    ("Membrane primary amine oxidase", "A0A8I5KP93", "MAOA", "REJECT_WRONG", "P21397",
     "Amine oxidase [flavin-containing] A",
     "膜结合原胺氧化酶是 AOC3/VAP-1（铜依赖），MAOA 是 FAD 依赖单胺氧化酶，机制不同"),
    ("Methylmalonyl-CoA mutase", "Q59HB4", None, "KEEP_TR", None, None,
     "片段型记录；人类 reviewed 为 MMUT（P22033），需按名称另核"),
    ("Orotidine 5'-phosphate decarboxylase", "B5LY70", "UMPS", "REPLACE", "P11172",
     "Uridine 5'-monophosphate synthase",
     "UMPS 为双功能酶，含乳清酸核苷酸脱羧酶活性，即该靶点"),
    ("Serine hydroxymethyltransferase", "Q53ET7", None, "KEEP_TR", None, None,
     "片段型记录；人类 SHMT1（P34896）/SHMT2（P34897）为 reviewed"),
    ("Serine hydroxymethyltransferase 1", "Q53ET7", None, "KEEP_TR", "P34896",
     "Serine hydroxymethyltransferase, cytosolic", "SHMT1 即胞质型，按名称可直接定"),
    ("Serine hydroxymethyltransferase 2", "Q53ET7", None, "KEEP_TR", "P34897",
     "Serine hydroxymethyltransferase, mitochondrial", "SHMT2 即线粒体型，按名称可直接定"),
    ("Tryptophanyl-tRNA synthetase 2", "G3V2Y7", "WARS1", "REJECT_WRONG", "Q9UGM6",
     "Tryptophan--tRNA ligase, mitochondrial",
     "靶点写的是 synthetase 2（线粒体 WARS2），TrEMBL 却给了胞质 WARS1，符号错配"),
    ("Uridine phosphorylase", "A0A0S2Z634", "UPP2", "AMBIGUOUS", "O95045",
     "Uridine phosphorylase 2",
     "人类有 UPP1 与 UPP2 两个同工酶，泛名无法唯一定位"),
]

COLS = ["target_name", "trembl_acc", "gene", "verdict", "canonical_acc",
        "canonical_name", "note"]


def main():
    df = pd.DataFrame(ROWS, columns=COLS)
    df.to_csv(TARGET / "uniprot_trembl_adjudicated.csv", index=False,
              encoding="utf-8-sig")
    counts = df.verdict.value_counts().to_dict()
    summary = dict(n_rows=len(df), verdicts=counts,
                   replace_accs=sorted(df.loc[df.verdict == "REPLACE",
                                              "canonical_acc"].dropna().tolist()),
                   resolved_by_name=sorted(df.loc[df.canonical_acc.notna() &
                                                   (df.verdict != "REPLACE"),
                                                   "canonical_acc"].dropna().tolist()))
    (TARGET / "uniprot_trembl_adjudicated.json").write_text(
        json.dumps(dict(summary=summary, rows=ROWS), ensure_ascii=False, indent=1),
        encoding="utf-8")
    print(df[["target_name", "gene", "verdict", "canonical_acc"]].to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
