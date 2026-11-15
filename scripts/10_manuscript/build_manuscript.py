# -*- coding: utf-8 -*-
import os, re, csv
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
MD = os.path.join(BASE, "论文写作", "paper_CBM_submission.md")
FIGDIR = os.path.join(BASE, "投稿图集")
SUPFIGDIR = os.path.join(BASE, "投稿图集", "补充图")
TABDIR = os.path.join(BASE, "投稿表格")
OUTDIR = os.path.join(BASE, "论文写作", "CBM_初审上传包_单独存放")

FIG_MAP = {
    "2.1 Overall design": [
        (FIGDIR, "Fig1_技术路线图.png",
         "Figure 1. Layered computational workflow, from compound and target cataloging through independent human validation, network and pathway analysis, cellular and spatial localization, structural simulation, and in silico perturbation."),
    ],
    "3.2 Core intersection": [
        (FIGDIR, "Fig4_PPI与富集.png",
         "Figure 4. Protein-protein interaction network and functional enrichment of the core targets."),
    ],
    "3.4 Cross-cohort direction": [
        (FIGDIR, "Fig2_多队列方向矩阵.png",
         "Figure 2. Multi-cohort direction matrix for the candidate genes across the five independent human cohorts."),
    ],
    "3.5 Meta-analysis": [
        (FIGDIR, "Fig3_Meta森林图拼版.png",
         "Figure 3. Random-effects meta-analysis forest plots for the core pain-related genes."),
    ],
    "3.7 Single-cell localization": [
        (FIGDIR, "Fig5_单细胞定位.png",
         "Figure 5. Single-cell localization of the core targets across human liver cell types."),
        (SUPFIGDIR, "S1_Harmony_UMAP.png",
         "Supplementary Figure S1. Harmony integration of the 20-sample MASLD collection, shown before and after batch correction."),
        (SUPFIGDIR, "S2_marker_heatmap.png",
         "Supplementary Figure S2. Marker-gene heatmap confirming the nine cell-type annotations in the MASLD collection."),
        (SUPFIGDIR, "S3_composition_stacked.png",
         "Supplementary Figure S3. Donor-level cell-type composition across the 20 samples, presented without inferential p-values."),
    ],
    "3.8 Spatial": [
        (FIGDIR, "Fig6_空间转录组.png",
         "Figure 6. Spatial transcriptomic localization of the target signals in human and mouse liver."),
    ],
    "3.9 Molecular docking": [
        (FIGDIR, "Fig7_分子对接.png",
         "Figure 7. Molecular docking poses and binding energies of the leading compounds."),
    ],
    "3.10 Molecular dynamics": [
        (FIGDIR, "Fig8_MD稳定性.png",
         "Figure 8. Molecular dynamics stability of the 5KIR-luteolin complex over the production trajectory."),
        (SUPFIGDIR, "S4_pose_clustering.png",
         "Supplementary Figure S4. Backbone and ligand pose clustering, showing a single protein cluster and six ligand clusters."),
    ],
    "3.11 In silico perturbation": [
        (FIGDIR, "Fig9_因果扰动.png",
         "Figure 9. In silico causal perturbation effects of transcription factor knockout and overexpression."),
    ],
    "3.12 The sender-to-receiver": [
        (SUPFIGDIR, "S5_send_receive.png",
         "Supplementary Figure S5. Sender-to-receiver panels linking liver mediators to their sensory-neuron receptors."),
        (SUPFIGDIR, "S6_drg_dotplot.png",
         "Supplementary Figure S6. Receptor-expression dot plot across the human dorsal root ganglion clusters."),
    ],
}

# ---------- low-level helpers ----------
def set_run_font(run, size=11, bold=False, italic=False, color=(0,0,0)):
    run.font.name = "Times New Roman"
    r = run._element
    rPr = r.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = OxmlElement('w:rFonts'); rPr.append(rFonts)
    rFonts.set(qn('w:ascii'), 'Times New Roman')
    rFonts.set(qn('w:hAnsi'), 'Times New Roman')
    rFonts.set(qn('w:eastAsia'), 'Times New Roman')
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = RGBColor(*color)

def add_tagged_text(p, text, size=11, bold=False):
    # handle <i>..</i>; strip all other tags
    parts = re.split(r"(<i>.*?</i>)", text)
    for part in parts:
        if not part:
            continue
        if part.startswith("<i>") and part.endswith("</i>"):
            run = p.add_run(part[3:-4]); set_run_font(run, size=size, bold=bold, italic=True)
        else:
            cleaned = re.sub(r"<[^>]+>", "", part)
            if cleaned:
                run = p.add_run(cleaned); set_run_font(run, size=size, bold=bold)

def shade_cell(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear'); shd.set(qn('w:color'), 'auto'); shd.set(qn('w:fill'), hexcolor)
    tcPr.append(shd)

def set_table_borders(table):
    tbl = table._tbl
    tblPr = tbl.tblPr
    borders = OxmlElement('w:tblBorders')
    for edge in ('top','left','bottom','right','insideH','insideV'):
        el = OxmlElement('w:'+edge)
        el.set(qn('w:val'),'single'); el.set(qn('w:sz'),'4')
        el.set(qn('w:space'),'0'); el.set(qn('w:color'),'808080')
        borders.append(el)
    tblPr.append(borders)

def repeat_header(row):
    trPr = row._tr.get_or_add_trPr()
    th = OxmlElement('w:tblHeader'); th.set(qn('w:val'),'true'); trPr.append(th)

def cant_split(row):
    trPr = row._tr.get_or_add_trPr()
    cs = OxmlElement('w:cantSplit'); cs.set(qn('w:val'),'true'); trPr.append(cs)

def build_table(doc, rows, widths=None, fontsize=9.5):
    ncol = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    set_table_borders(t)
    for ri, rowdata in enumerate(rows):
        row = t.rows[ri]
        cant_split(row)
        if ri == 0:
            repeat_header(row)
        for ci, val in enumerate(rowdata):
            cell = row.cells[ci]
            cell.vertical_alignment = 1  # center
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1.0
            p.paragraph_format.first_line_indent = Pt(0)
            run = p.add_run(str(val))
            set_run_font(run, size=fontsize, bold=(ri == 0))
            if ri == 0:
                shade_cell(cell, "D9D9D9")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if ri == 0 else WD_ALIGN_PARAGRAPH.LEFT
    return t

def add_figure(doc, directory, file, caption):
    path = os.path.join(directory, file)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(8); p.paragraph_format.space_after = Pt(2)
    run = p.add_run()
    run.add_picture(path, width=Cm(15.5))
    cp = doc.add_paragraph(); cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cp.paragraph_format.space_before = Pt(2); cp.paragraph_format.space_after = Pt(10)
    add_tagged_text(cp, caption, size=10)

# ---------- document setup ----------
doc = Document()
sec = doc.sections[0]
sec.page_width = Cm(21); sec.page_height = Cm(29.7)
sec.top_margin = sec.bottom_margin = Cm(2.5)
sec.left_margin = sec.right_margin = Cm(2.5)

normal = doc.styles['Normal']
normal.font.name = 'Times New Roman'; normal.font.size = Pt(11)
normal.paragraph_format.line_spacing = 1.5
normal.paragraph_format.space_before = Pt(0); normal.paragraph_format.space_after = Pt(0)

lines = open(MD, encoding="utf-8").read().split("\n")
# split at references
ref_idx = next(i for i,l in enumerate(lines) if l.strip() == "## References")
body_lines = lines[:ref_idx]
ref_lines = lines[ref_idx+1:]

pending_fig = []
def flush_figure():
    global pending_fig
    for directory, file, caption in pending_fig:
        add_figure(doc, directory, file, caption)
    pending_fig = []

i = 0
in_author_block = True
while i < len(body_lines):
    line = body_lines[i]
    s = line.strip()
    if not s:
        i += 1; continue
    if s.startswith("### "):
        flush_figure()
        head = s[4:]
        p = doc.add_paragraph()
        run = p.add_run(head); set_run_font(run, size=12, bold=True)
        p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(4)
        match_head = re.sub(r'^(\d+(?:\.\d+)*)\.\s+', r'\1 ', head)
        for key, val in FIG_MAP.items():
            if match_head.startswith(key):
                pending_fig = list(val)
        in_author_block = False
    elif s.startswith("## "):
        flush_figure()
        head = s[3:]
        p = doc.add_paragraph()
        run = p.add_run(head); set_run_font(run, size=14, bold=True)
        p.paragraph_format.space_before = Pt(14); p.paragraph_format.space_after = Pt(6)
        in_author_block = False
    elif s.startswith("# "):
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(s[2:]); set_run_font(run, size=16, bold=True)
        p.paragraph_format.space_after = Pt(14)
    elif s.startswith("|"):
        # collect markdown table
        block = []
        while i < len(body_lines) and body_lines[i].strip().startswith("|"):
            block.append(body_lines[i].strip()); i += 1
        rows = []
        for bi, b in enumerate(block):
            if re.fullmatch(r"[\|\s\-:]+", b):
                continue
            cells = [c.strip() for c in b.strip("|").split("|")]
            rows.append(cells)
        build_table(doc, rows, fontsize=9.5)
        doc.add_paragraph().paragraph_format.space_after = Pt(4)
        continue
    else:
        p = doc.add_paragraph()
        if in_author_block:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_tagged_text(p, s, size=11)
        else:
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            add_tagged_text(p, s, size=11)
    i += 1
flush_figure()

# ---------- references ----------
p = doc.add_paragraph()
run = p.add_run("References"); set_run_font(run, size=14, bold=True)
p.paragraph_format.space_before = Pt(14); p.paragraph_format.space_after = Pt(6)
cur = ""
entries = []
for l in ref_lines:
    if not l.strip():
        if cur: entries.append(cur); cur = ""
        continue
    if re.match(r"^\[\d+\]", l.strip()):
        if cur: entries.append(cur)
        cur = l.strip()
    else:
        cur += " " + l.strip()
if cur: entries.append(cur)
for e in entries:
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.line_spacing = 1.15
    p.paragraph_format.space_after = Pt(4)
    add_tagged_text(p, e, size=10)

# ---------- additional end tables (3,4,5) ----------
def read_csv_rows(name):
    with open(os.path.join(TABDIR, name), encoding="utf-8-sig") as f:
        return list(csv.reader(f))

doc.add_page_break()
extra = [
    ("Table 1. Cohort characteristics and the role each dataset played in the analysis.", "Table1_队列特征.csv", 9),
    ("Table 2. Evidence-tier classification of the candidate targets.", "Table3_靶点证据分级.csv", 8.5),
    ("Table 3. Cross-cohort direction voting for the candidate genes.", "Table2_方向投票.csv", 9),
    ("Table 4. Molecular docking binding energies (kcal/mol).", "Table4_对接结合能.csv", 9),
    ("Table 5. Random-effects meta-analysis effects with Benjamini-Hochberg adjusted p values.", "Table5_Meta效应BH.csv", 8.5),
]
for title, fname, fs in extra:
    p = doc.add_paragraph()
    run = p.add_run(title); set_run_font(run, size=11, bold=True)
    p.paragraph_format.space_before = Pt(6); p.paragraph_format.space_after = Pt(4)
    rows = read_csv_rows(fname)
    build_table(doc, rows, fontsize=fs)
    doc.add_paragraph().paragraph_format.space_after = Pt(6)

out = os.path.join(OUTDIR, "Manuscript_CBM.docx")
doc.save(out)
print("saved", out)
print("paragraphs", len(doc.paragraphs), "tables", len(doc.tables))
