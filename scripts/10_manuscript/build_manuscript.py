# -*- coding: utf-8 -*-
import os, re, csv
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from lxml import etree

BASE = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
MD = os.path.join(BASE, "论文写作", "paper_CMPB_submission.md")
FIGDIR = os.path.join(BASE, "投稿图集")
SUPFIGDIR = os.path.join(BASE, "投稿图集", "补充图")
TABDIR = os.path.join(BASE, "投稿表格")
OUTDIR = os.path.join(BASE, "论文写作", "CMPB_初审上传包")

FIG_MAP = {
    "2.1": [
        (FIGDIR, "Fig1_技术路线图.png",
         "Figure 1. Layered computational workflow, from compound and target cataloging through independent human validation, network and pathway analysis, cellular and spatial localization, structural simulation, and in silico perturbation."),
    ],
    "3.2": [
        (FIGDIR, "Fig4_PPI与富集.png",
         "Figure 2. Protein-protein interaction network and functional enrichment of the core targets."),
    ],
    "3.4": [
        (FIGDIR, "Fig2_多队列方向矩阵.png",
         "Figure 3. Multi-cohort direction matrix for the candidate genes across the five independent human cohorts."),
    ],
    "3.5": [
        (FIGDIR, "Fig3_Meta森林图拼版.png",
         "Figure 4. Random-effects meta-analysis forest plots for the core pain-related genes."),
    ],
    "3.6": [
        (FIGDIR, "Fig5_单细胞定位.png",
         "Figure 5. Single-cell localization of the core targets across human liver cell types."),
        (os.path.join(BASE, "healthy_liver"), "fig_external_healthy_reference.png",
         "Figure 6. External healthy human liver reference: expressing-cell percentages versus the MASLD donor median, with the single GSE344087 control shown for comparison."),
        (SUPFIGDIR, "S7_suspension_sensitivity.png",
         "Supplementary Figure S1. Whole-cell versus nuclear suspension sensitivity: the healthy-vs-disease contrast for IL1B, TNF, and PTGS2 holds in both suspensions, while FOS is strongly suspension-dependent."),
        (SUPFIGDIR, "S1_Harmony_UMAP.png",
         "Supplementary Figure S2. Harmony integration of the 20-sample MASLD collection, shown before and after batch correction."),
        (SUPFIGDIR, "S2_marker_heatmap.png",
         "Supplementary Figure S3. Marker-gene heatmap confirming the nine cell-type annotations in the MASLD collection."),
        (SUPFIGDIR, "S3_composition_stacked.png",
         "Supplementary Figure S4. Donor-level cell-type composition across the 20 samples, presented without inferential p-values."),
    ],
    "3.7": [
        (FIGDIR, "Fig6_空间转录组.png",
         "Figure 7. Spatial transcriptomic localization of the target signals in human and mouse liver."),
    ],
    "3.8": [
        (FIGDIR, "Fig7_分子对接.png",
         "Figure 8. Molecular docking poses and binding energies of the leading compounds."),
        (FIGDIR, "Fig8_MD稳定性.png",
         "Figure 9. Molecular dynamics stability of the 5KIR-luteolin complex over the production trajectory."),
        (SUPFIGDIR, "S4_pose_clustering.png",
         "Supplementary Figure S5. Backbone and ligand pose clustering, showing a single protein cluster and six ligand clusters."),
    ],
    "3.9": [
        (FIGDIR, "Fig9_因果扰动.png",
         "Figure 10. In silico causal perturbation effects of transcription factor knockout and overexpression."),
    ],
    "3.10": [
        (SUPFIGDIR, "S5_send_receive.png",
         "Supplementary Figure S6. Sender-to-receiver panels linking liver mediators to their sensory-neuron receptors, with edge evidence graded."),
        (SUPFIGDIR, "S6_drg_dotplot.png",
         "Supplementary Figure S7. Receptor-expression dot plot across the human dorsal root ganglion clusters."),
    ],
    "3.18": [
        (os.path.join(BASE, "chaihu-nociceptive-mapping", "figures"),
         "wp7_cross_replicate.png",
         "Figure 11. Cross-replicate structural robustness over three independent 20 ns molecular dynamics runs, including backbone and ligand RMSD, pocket retention, residue flexibility, and block-based MM-GBSA convergence."),
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

TOK_RE = re.compile(r"(<i>.*?</i>|\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]*\^)")
TOK_AUTH_RE = re.compile(r"(<i>.*?</i>|\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]*\^|\*)")

def add_tagged_text(p, text, size=11, bold=False, author=False):
    for part in (TOK_AUTH_RE if author else TOK_RE).split(text):
        if not part:
            continue
        if part.startswith("<i>") and part.endswith("</i>"):
            run = p.add_run(part[3:-4]); set_run_font(run, size=size, bold=bold, italic=True)
        elif part.startswith("**") and part.endswith("**"):
            run = p.add_run(part[2:-2]); set_run_font(run, size=size, bold=True)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            run = p.add_run(part[1:-1]); set_run_font(run, size=size, bold=bold, italic=True)
        elif part.startswith("^") and part.endswith("^") and len(part) > 2:
            run = p.add_run(part[1:-1]); set_run_font(run, size=size, bold=bold)
            run.font.superscript = True
        elif part == "*":
            run = p.add_run("*"); set_run_font(run, size=size, bold=bold)
            run.font.superscript = True
        else:
            cleaned = re.sub(r"<[^>]+>", "", part)
            if cleaned:
                run = p.add_run(cleaned); set_run_font(run, size=size, bold=bold)

def set_table_borders(table):
    # Three-line table: thick top and bottom rules only. Vertical and
    # internal horizontal rules are removed per CMPB table requirements.
    tbl = table._tbl
    tblPr = tbl.tblPr
    borders = OxmlElement('w:tblBorders')
    for edge in ('top','left','bottom','right','insideH','insideV'):
        el = OxmlElement('w:'+edge)
        if edge in ('top', 'bottom'):
            el.set(qn('w:val'),'single'); el.set(qn('w:sz'),'12')
            el.set(qn('w:space'),'0'); el.set(qn('w:color'),'000000')
        else:
            el.set(qn('w:val'),'nil')
        borders.append(el)
    tblPr.append(borders)

def set_header_rule(row):
    # Thin rule below the header row. w:sz is in eighths of a point,
    # so sz = 6 gives 0.75 pt.
    for cell in row.cells:
        tcPr = cell._tc.get_or_add_tcPr()
        old = tcPr.find(qn('w:tcBorders'))
        if old is not None:
            tcPr.remove(old)
        borders = OxmlElement('w:tcBorders')
        bottom = OxmlElement('w:bottom')
        bottom.set(qn('w:val'),'single'); bottom.set(qn('w:sz'),'6')
        bottom.set(qn('w:space'),'0'); bottom.set(qn('w:color'),'000000')
        borders.append(bottom)
        tcPr.append(borders)

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
            set_header_rule(row)
        for ci, val in enumerate(rowdata):
            cell = row.cells[ci]
            cell.vertical_alignment = 1
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1.0
            p.paragraph_format.first_line_indent = Pt(0)
            run = p.add_run(str(val))
            set_run_font(run, size=fontsize, bold=(ri == 0))
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
        sec_num = re.match(r'^(\d+(?:\.\d+)*)\.', head)
        if sec_num and sec_num.group(1) in FIG_MAP:
            pending_fig = list(FIG_MAP[sec_num.group(1)])
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
            add_tagged_text(p, s, size=11, author=True)
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

# ---------- end tables ----------
_CJK = re.compile(r"[一-鿿]")
def read_csv_rows(directory, name):
    with open(os.path.join(directory, name), encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    bad = [c for r in rows for c in r if _CJK.search(c)]
    if bad:
        raise ValueError("Non-English (CJK) cells found in " + name + ": " + "; ".join(bad[:5]))
    return rows

doc.add_page_break()
extra = [
    ("Table 1. Cohort characteristics and the role each dataset played in the analysis.", TABDIR, "Table1_队列特征.csv", 9),
    ("Table 2. Evidence-tier classification of the candidate targets.", TABDIR, "Table3_靶点证据分级.csv", 8.5),
    ("Table 3. Cross-cohort direction voting for the candidate genes.", TABDIR, "Table2_方向投票.csv", 9),
    ("Table 4. Molecular docking binding energies (kcal/mol).", TABDIR, "Table4_对接结合能.csv", 9),
    ("Table 5. Random-effects meta-analysis effects with Benjamini-Hochberg adjusted p values.", TABDIR, "Table5_Meta效应BH.csv", 8.5),
    ("Table 6. Inference ledger: evidence type, status, corresponding figure or table, and the experiment that would close each step.", TABDIR, "Table6_推断台账.csv", 8),
    ("Table 7. Cross-formula composition and de-circularized benchmark inputs for the two formulas.", OUTDIR, "Table7_CrossFormula_Composition.csv", 9),
    ("Table 8. Convergence of the strict core genes across the two formulas.", OUTDIR, "Table8_CoreConvergence.csv", 9.5),
    ("Table 9. Component-level validation of the framework, including the failed end-to-end generalization test.", OUTDIR, "Table9_ComponentValidation.csv", 8.5),
]
for title, directory, fname, fs in extra:
    p = doc.add_paragraph()
    run = p.add_run(title); set_run_font(run, size=11, bold=True)
    p.paragraph_format.space_before = Pt(6); p.paragraph_format.space_after = Pt(4)
    rows = read_csv_rows(directory, fname)
    build_table(doc, rows, fontsize=fs)
    doc.add_paragraph().paragraph_format.space_after = Pt(6)

def strip_cjk_metadata(document):
    # Picture internal names inherit the source PNG file names, which are Chinese.
    n = 0
    for elem in document.element.iter():
        if etree.QName(elem).localname == "cNvPr":
            n += 1
            elem.set("name", "Image %d" % n)
            descr = elem.get("descr")
            if descr and _CJK.search(descr):
                elem.set("descr", "Image %d" % n)
    # The default template keeps east-asian font slots in the font table and theme.
    for part in document.part.package.iter_parts():
        s = str(part.partname)
        if s == "/word/fontTable.xml" or s.startswith("/word/theme/"):
            el = getattr(part, "_element", None)
            if el is not None:
                for e in el.iter():
                    for a in list(e.attrib):
                        if _CJK.search(e.attrib[a] or ""):
                            e.attrib[a] = "Calibri"
                    if e.text and _CJK.search(e.text):
                        e.text = "Calibri"
            else:
                # Generic parts keep their bytes in .blob rather than an lxml tree.
                part._blob = _CJK.sub("Calibri", part.blob.decode("utf-8")).encode("utf-8")

strip_cjk_metadata(doc)

out = os.path.join(OUTDIR, "Manuscript_CMPB.docx")
doc.save(out)
print("saved", out)
print("paragraphs", len(doc.paragraphs), "tables", len(doc.tables))
