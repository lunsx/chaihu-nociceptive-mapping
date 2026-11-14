# -*- coding: utf-8 -*-
import os
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\论文写作"
OUTDIR = os.path.join(BASE, "CBM_初审上传包_单独存放")

def style_run(run, size=11, bold=False, italic=False):
    run.font.name = "Times New Roman"
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = OxmlElement('w:rFonts'); rPr.append(rFonts)
    for a in ('w:ascii','w:hAnsi','w:eastAsia'):
        rFonts.set(qn(a), 'Times New Roman')
    run.font.size = Pt(size); run.font.bold = bold; run.font.italic = italic
    run.font.color.rgb = RGBColor(0,0,0)

def new_doc():
    d = Document()
    s = d.sections[0]
    s.page_width = Cm(21); s.page_height = Cm(29.7)
    s.top_margin = s.bottom_margin = Cm(2.5); s.left_margin = s.right_margin = Cm(2.5)
    d.styles['Normal'].font.name = 'Times New Roman'
    d.styles['Normal'].font.size = Pt(11)
    return d

# ---------- cover letter ----------
d = new_doc()
for raw in open(os.path.join(BASE, "cover_letter_CBM.md"), encoding="utf-8"):
    line = raw.rstrip("\n")
    if not line.strip():
        continue
    p = d.add_paragraph()
    p.paragraph_format.line_spacing = 1.15
    p.paragraph_format.space_after = Pt(8)
    r = p.add_run(line); style_run(r)
d.save(os.path.join(OUTDIR, "Cover_Letter_CBM.docx"))
print("cover letter saved")

# ---------- highlights ----------
d = new_doc()
h = d.add_paragraph()
r = h.add_run("Highlights"); style_run(r, size=14, bold=True)
h.paragraph_format.space_after = Pt(8)
for raw in open(os.path.join(BASE, "highlights_CBM.md"), encoding="utf-8"):
    line = raw.strip()
    if line.startswith("- "):
        p = d.add_paragraph(style='List Bullet')
        p.paragraph_format.space_after = Pt(4)
        r = p.add_run(line[2:]); style_run(r)
d.save(os.path.join(OUTDIR, "Highlights_CBM.docx"))
print("highlights saved")
