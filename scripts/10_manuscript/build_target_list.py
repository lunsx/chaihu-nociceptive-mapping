# -*- coding: utf-8 -*-
import os, csv, re, json
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

root = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
OUT = os.path.join(root, "论文写作", "CBM_初审上传包_单独存放", "靶点清单_完整版.docx")
recs = json.load(open(os.path.join(root, "target_list_enriched.json"), encoding="utf-8"))

# ---- clean reference lines from manuscript ----
md = open(os.path.join(root, "论文写作", "paper_CBM_submission.md"), encoding="utf-8").read().split("\n")
ri = next(i for i, l in enumerate(md) if l.strip() == "## References")
raw, cur = {}, ""
for l in md[ri+1:]:
    s = l.strip()
    if not s:
        if cur:
            n = int(re.match(r"\[(\d+)\]", cur).group(1)); raw[n] = cur; cur = ""
        continue
    if re.match(r"^\[\d+\]", s):
        if cur:
            n = int(re.match(r"^\[(\d+)\]", cur).group(1)); raw[n] = cur
        cur = s
    else:
        cur += " " + s
if cur:
    raw[int(re.match(r"^\[(\d+)\]", cur).group(1))] = cur
def clean(t):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", t)).strip()
refs = {n: clean(t) for n, t in raw.items()}

def set_font(run, size, bold=False, italic=False, cn="宋体", en="Times New Roman"):
    run.font.name = en
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = OxmlElement('w:rFonts'); rPr.append(rFonts)
    rFonts.set(qn('w:ascii'), en); rFonts.set(qn('w:hAnsi'), en); rFonts.set(qn('w:eastAsia'), cn)
    run.font.size = Pt(size); run.font.bold = bold; run.font.italic = italic
    run.font.color.rgb = RGBColor(0, 0, 0)

def shade(cell, hexc):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd'); shd.set(qn('w:val'),'clear'); shd.set(qn('w:color'),'auto'); shd.set(qn('w:fill'),hexc)
    tcPr.append(shd)

def borders(table):
    pr = table._tbl.tblPr; b = OxmlElement('w:tblBorders')
    for e in ('top','left','bottom','right','insideH','insideV'):
        x = OxmlElement('w:'+e); x.set(qn('w:val'),'single'); x.set(qn('w:sz'),'4')
        x.set(qn('w:space'),'0'); x.set(qn('w:color'),'808080'); b.append(x)
    pr.append(b)

def fixed_layout(table):
    tblPr = table._tbl.tblPr
    lay = OxmlElement('w:tblLayout'); lay.set(qn('w:type'),'fixed'); tblPr.append(lay)

def hdr_row(row):
    trPr = row._tr.get_or_add_trPr()
    th = OxmlElement('w:tblHeader'); th.set(qn('w:val'),'true'); trPr.append(th)
    cs = OxmlElement('w:cantSplit'); cs.set(qn('w:val'),'true'); trPr.append(cs)

def split_row(row):
    trPr = row._tr.get_or_add_trPr()
    cs = OxmlElement('w:cantSplit'); cs.set(qn('w:val'),'true'); trPr.append(cs)

doc = Document()
sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = Cm(29.7), Cm(21)
sec.top_margin = sec.bottom_margin = Cm(1.6); sec.left_margin = sec.right_margin = Cm(1.8)
doc.styles['Normal'].font.size = Pt(10)

def para(text, size=10, bold=False, cn="宋体", align=None, before=2, after=2):
    p = doc.add_paragraph()
    if align is not None: p.alignment = align
    p.paragraph_format.space_before = Pt(before); p.paragraph_format.space_after = Pt(after)
    r = p.add_run(text); set_font(r, size, bold=bold, cn=cn)
    return p

def heading(text, size=13):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text); set_font(r, size, bold=True, cn="黑体")
    return p

# ---------- title block ----------
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("柴胡疏肝散（醋香附入药）作用靶点清单"); set_font(r, 18, bold=True, cn="黑体")
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Computational target list of Chaihu Shugan San (with vinegar-processed Cyperi Rhizoma)")
set_font(r, 11, bold=True, cn="黑体")
para("学生：刘丰恺    导师：高天慧    单位：齐鲁医药学院药学院    基金：国家自然科学基金 82505032",
     10, align=WD_ALIGN_PARAGRAPH.CENTER, before=4, after=2)

# ---------- core conclusion ----------
heading("一、核心结论")
for t in [
    "1. 严格三者交集（疾病差异基因 × 方剂成分靶点 × 人类蛋白）收敛到 6 个基因：PTGS2、IL6、IL1B、IL1A、PTGES、PTGER3，指向 IL1B–NF-kB–PTGS2–PGE2 外周炎症痛通路。",
    "2. 全部 60 个候选靶点按证据分为 T0–T4 五级，T1 为核心、T2 为强支撑；分级、方向、一致率和细胞定位见下表。",
    "3. TRPV1 在筛选成分中零命中、多队列方向不稳健，证据不支持其作为核心，已在表中如实标注为负向/弱证据。",
]:
    para(t, 10)

# ---------- core gene quantitative table ----------
heading("二、核心基因量化指标")
core = ["PTGS2","IL6","IL1B","IL1A","PTGES","PTGER3","TRPV1"]
cmap = {x["gene"]: x for x in recs}
cols = ["基因","GSE89632 logFC","GSE89632 adj_p","GSE130970 r","Meta 合并 g","Meta BH p","筛选成分数"]
t = doc.add_table(rows=1+len(core), cols=len(cols)); t.alignment = WD_TABLE_ALIGNMENT.CENTER
borders(t); fixed_layout(t); hdr_row(t.rows[0])
for ci, c in enumerate(cols):
    cell = t.rows[0].cells[ci]; p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(c); set_font(r, 9, bold=True, cn="黑体"); shade(cell, "D9D9D9")
for rii, g in enumerate(core, start=1):
    x = cmap[g]; split_row(t.rows[rii])
    vals = [g, x["logFC"], ("%.2e" % x["adj_p"] if x["adj_p"] != "" else ""),
            x["gse130970_r"], x["meta_g"], ("%.3f" % x["meta_bh"] if x["meta_bh"] != "" else ""),
            x["n_filtered_compounds"]]
    for ci, v in enumerate(vals):
        cell = t.rows[rii].cells[ci]; p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(str(v)); set_font(r, 9, bold=(ci == 0))

# ---------- legends ----------
heading("三、分级与来源说明")
para("分级：T0 上游转录因子；T1 核心（多库验证）；T2 强支撑；T3 证据弱或负向；T4 其他候选。", 9.5)
para("数据来源：TCMSP＝方剂成分靶点；GSE89632＝疾病差异表达；痛觉基因集＝外周痛相关基因；HERB/TTD＝多库交叉验证；多队列Meta＝五队列随机效应整合；单细胞定位＝人肝单细胞转录组。", 9.5)

# ---------- full 60 table ----------
heading("四、完整 60 靶点明细")
headers = ["序号","基因","蛋白全称","UniProt","Entrez","TCMSP","分级","方向","一致率","单细胞定位","筛选成分","涉及药味","数据来源","关键文献"]
widths = [0.8,1.5,3.4,1.4,1.2,1.2,1.7,0.9,1.0,2.9,1.0,2.1,3.3,1.3]
order = {"T0":0,"T1":1,"T2":2,"T3":3,"T4":4}
def trank(x):
    for k, v in order.items():
        if x["tier"].startswith(k): return v
    return 9
data = sorted(recs, key=lambda x: (trank(x), x["gene"]))
t = doc.add_table(rows=1+len(data), cols=len(headers)); t.alignment = WD_TABLE_ALIGNMENT.CENTER
borders(t); fixed_layout(t); hdr_row(t.rows[0])
for ci, (h, w) in enumerate(zip(headers, widths)):
    cell = t.rows[0].cells[ci]; cell.width = Cm(w)
    p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(h); set_font(r, 8.5, bold=True, cn="黑体"); shade(cell, "D9D9D9")
for idx, x in enumerate(data, start=1):
    row = t.rows[idx]; split_row(row)
    vals = [idx, x["gene"], x["protein"], x["uniprot"], x["entrez"], x["tcmsp_target_id"],
            x["tier"], x["direction"], x["concordance"], x["celltype"], x["n_filtered_compounds"],
            x["herbs"], x["sources"], x["literature"]]
    for ci, (v, w) in enumerate(zip(vals, widths)):
        cell = row.cells[ci]; cell.width = Cm(w)
        p = cell.paragraphs[0]; p.paragraph_format.line_spacing = 1.0
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if ci in (0,1,3,4,5,7,8,10,13) else WD_ALIGN_PARAGRAPH.LEFT
        r = p.add_run(str(v)); set_font(r, 8, bold=(ci == 1))

# ---------- cited literature ----------
heading("五、关键文献（编号对应投稿论文参考文献）")
for n in [3,5,8,15,16,19,34,35]:
    if n in refs:
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(3); p.paragraph_format.line_spacing = 1.1
        r = p.add_run(refs[n]); set_font(r, 8.5)

doc.save(OUT)
print("saved", OUT)
