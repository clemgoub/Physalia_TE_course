#!/usr/bin/env python3
"""Draws the Practical 3 cartoons as SVG files.

Each variant cartoon shows the reference and the two HG002 haplotypes around
one GraffiTE record, then the VCF fields that describe it. The sequences,
coordinates and field values come from the 3 Nov 2026 run on the course
server (~ubuntu/practical3_precomputed). Run: python3 make_figures.py
"""
from html import escape
from pathlib import Path

OUT = Path(__file__).parent
W = 960
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"
SANS = "-apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"

C = {
    "flank": ("#e6e9ee", "#9aa3ad"),
    "alu": ("#9fdccf", "#1f7a6f"),
    "l1": ("#cdbdf2", "#5b3fa8"),
    "sat": ("#d9d9d9", "#6b6b6b"),
    "tsd": ("#ffd27a", "#a87300"),
    "polya": ("#f7b8aa", "#b5462f"),
    "plain": ("#ffffff", "#9aa3ad"),
}
TEXT = "#24292f"
MUTED = "#57606a"


def t(x, y, s, size=13, anchor="start", weight="normal", fam=SANS, fill=TEXT):
    return (f'<text x="{x}" y="{y}" font-family="{fam}" font-size="{size}" '
            f'text-anchor="{anchor}" font-weight="{weight}" fill="{fill}">{escape(s)}</text>')


def spans(x, y, parts, size=12.5):
    """One line of monospace text made of (text, colour, bold) parts."""
    out = [f'<text x="{x}" y="{y}" font-family="{MONO}" font-size="{size}" fill="{TEXT}" xml:space="preserve">']
    for s, col, bold in parts:
        out.append(f'<tspan fill="{col or TEXT}" font-weight="{"bold" if bold else "normal"}">{escape(s)}</tspan>')
    out.append("</text>")
    return "".join(out)


def seg(x1, x2, y, kind, label="", sub="", direction=None, h=26, size=12, fam=SANS):
    """A stretch of DNA drawn as a box, or as an arrow when it has a strand."""
    fill, stroke = C[kind]
    p = 9 if direction else 0
    if direction == "+":
        pts = f"{x1},{y} {x2 - p},{y} {x2},{y + h / 2} {x2 - p},{y + h} {x1},{y + h}"
    elif direction == "-":
        pts = f"{x1 + p},{y} {x2},{y} {x2},{y + h} {x1 + p},{y + h} {x1},{y + h / 2}"
    else:
        pts = f"{x1},{y} {x2},{y} {x2},{y + h} {x1},{y + h}"
    out = [f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>']
    if kind == "sat":
        out.append(f'<polygon points="{pts}" fill="url(#hatch)" stroke="none"/>')
    if label:
        out.append(t((x1 + x2) / 2, y + h / 2 + 4.5, label, size, "middle", "bold" if kind not in ("flank", "plain") else "normal", fam))
    if sub:
        out.append(t((x1 + x2) / 2, y + h + 15, sub, 11, "middle", fill=MUTED))
    return "".join(out)


def gap(x1, x2, y, h=26):
    """A stretch missing from this haplotype: a dashed line where the DNA would be."""
    return (f'<line x1="{x1}" y1="{y + h / 2}" x2="{x2}" y2="{y + h / 2}" stroke="{MUTED}" '
            f'stroke-width="1.5" stroke-dasharray="5 4"/>')


def brace(x1, x2, y, label):
    return (f'<path d="M{x1},{y + 8} L{x1},{y} L{x2},{y} L{x2},{y + 8}" fill="none" stroke="{MUTED}" stroke-width="1.2"/>'
            + t((x1 + x2) / 2, y - 6, label, 12, "middle", fill=MUTED))


def svg(height, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {height}" width="{W}" height="{height}" role="img">'
            f"<title>{escape(title)}</title>"
            '<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
            '<line x1="0" y1="0" x2="0" y2="6" stroke="#9a9a9a" stroke-width="1.2"/></pattern>'
            '<marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M1 1L9 5L1 9" fill="none" stroke="{MUTED}" stroke-width="1.6"/></marker></defs>'
            f'<rect width="{W}" height="{height}" fill="#ffffff"/>' + "".join(body) + "</svg>")


def vcf_panel(y, rows, head):
    """rows: list of (parts, explanation). Returns (svg, bottom y)."""
    line = 22
    h = 44 + line * len(rows)
    out = [f'<rect x="16" y="{y}" width="{W - 32}" height="{h}" rx="8" fill="#f6f8fa" stroke="#d0d7de"/>',
           t(32, y + 24, head, 13, weight="bold")]
    yy = y + 48
    for parts, expl in rows:
        out.append(spans(32, yy, parts))
        if expl:
            out.append(t(500, yy, expl, 12.5, fill=MUTED))
        yy += line
    return "".join(out), y + h


def graph_inset(x, y, ref_has_te, kind, te_label, n_ref, n_alt, gt):
    """Graph bubble with the TE path and the no-TE path, and the reads on each."""
    out = [t(x, y, "Step 4: the reads on each path of the graph", 13, weight="bold")]
    yb = y + 78
    out.append(seg(x, x + 80, yb, "flank"))
    out.append(seg(x + 330, x + 410, yb, "flank"))
    # upper path: with the TE; lower path: without
    fill, stroke = C[kind]
    out.append(f'<path d="M{x + 80},{yb + 13} C{x + 110},{yb - 30} {x + 140},{yb - 30} {x + 160},{yb - 30}" fill="none" stroke="{stroke}" stroke-width="2"/>')
    out.append(seg(x + 160, x + 250, yb - 42, kind, te_label, size=11))
    out.append(f'<path d="M{x + 250},{yb - 30} C{x + 270},{yb - 30} {x + 300},{yb - 30} {x + 330},{yb + 13}" fill="none" stroke="{stroke}" stroke-width="2"/>')
    out.append(f'<path d="M{x + 80},{yb + 13} C{x + 130},{yb + 60} {x + 280},{yb + 60} {x + 330},{yb + 13}" fill="none" stroke="{MUTED}" stroke-width="2"/>')
    n_with = n_ref if ref_has_te else n_alt
    n_without = n_alt if ref_has_te else n_ref
    for i in range(n_with):
        yy = yb - 52 - 7 * i
        out.append(f'<line x1="{x + 170 + 12 * i}" y1="{yy}" x2="{x + 235 + 12 * i}" y2="{yy}" stroke="#24292f" stroke-width="2.5"/>')
    for i in range(n_without):
        yy = yb + 56 + 7 * i
        out.append(f'<line x1="{x + 160 + 12 * i}" y1="{yy}" x2="{x + 225 + 12 * i}" y2="{yy}" stroke="#24292f" stroke-width="2.5"/>')
    out.append(t(x + 425, yb - 26, f"path with the TE: {n_with} read{'s' * (n_with != 1)}", 12))
    out.append(t(x + 425, yb + 64, f"path without the TE: {n_without} read{'s' * (n_without != 1)}", 12))
    which = "the reference (0)" if ref_has_te else "the alternative allele (1)"
    out.append(t(x + 425, yb + 20, f"The TE path is {which}: GT={gt}", 12, weight="bold"))
    return "".join(out), yb + 84


def rows_block(y, rows):
    out = []
    for i, (name, items) in enumerate(rows):
        yy = y + i * 46
        out.append(t(20, yy + 18, name, 13, weight="bold" if i == 0 else "normal"))
        out.extend(f(yy) for f in items)
    return out, y + len(rows) * 46


# --- Figure: DEL.71_275, the variant followed through the steps -----------------
def fig_del71():
    b = [t(20, 30, "One variant through the steps: HG002.mat.svim_asm.DEL.71_275 (chr22:5,912,711)", 16, weight="bold"),
         t(20, 52, "An AluYb8 in the reference genome, with its target site duplication (TSD) and poly(A) tail. Both HG002 haplotypes lack it.", 13, fill=MUTED)]
    ref = [lambda y: seg(170, 260, y, "flank", "…TTTTTCA", fam=MONO),
           lambda y: seg(260, 370, y, "tsd", "AAGAAAATTAGTG", size=10.5, fam=MONO),
           lambda y: seg(370, 690, y, "alu", "AluYb8  (+ strand)", direction="+"),
           lambda y: seg(690, 750, y, "polya", "A×24", size=11),
           lambda y: seg(750, 860, y, "tsd", "AAGAAAATTAGTG", size=10.5, fam=MONO),
           lambda y: seg(860, 945, y, "flank", "TTCCTT…", fam=MONO)]
    hap = [lambda y: seg(170, 260, y, "flank", "…TTTTTCA", fam=MONO), lambda y: gap(260, 750, y),
           lambda y: seg(750, 860, y, "tsd", "AAGAAAATTAGTG", size=10.5, fam=MONO),
           lambda y: seg(860, 945, y, "flank", "TTCCTT…", fam=MONO)]
    b.append(brace(260, 750, 96, "REF allele: 327 bp missing from HG002 (plus the base before it, at POS)"))
    body, y = rows_block(110, [("Reference (CHM13)", ref), ("HG002 maternal", hap), ("HG002 paternal", hap)])
    b += body
    b.append(t(505, y - 8, "one TSD copy left: the empty site", 11.5, "middle", fill=MUTED))
    panel, y2 = vcf_panel(y + 14, [
        ([("REF=", None, False), ("T", MUTED, False), ("AAGAAAATTAGTG", C["tsd"][1], True), ("GGCC…", C["alu"][1], True), ("AAAA", C["polya"][1], True), ("  ALT=T", None, False)],
         "REF is the reference sequence (328 bp), ALT what HG002 has"),
        ([("SVLEN=-327", None, True)], "negative: a deletion relative to the reference"),
        ([("repeat_ids=AluYb8 matching_classes=SINE/Alu", C["alu"][1], True)], "RepeatMasker: one AluYb8 (n_hits=1)"),
        ([("RM_hit_strands=+ match_lengths=322", C["alu"][1], True)], "on the + strand; 322 bp of the variant match the Alu"),
        ([("TSD=AAGAAAATTAGTG,AAGAAAATTAGTG", C["tsd"][1], True)], "the same 13 bp on both sides of the element"),
        ([("polyA=TRUE", C["polya"][1], True)], "a poly(A) tail at the 3' end of the Alu"),
        ([("total_repeat_span=0.987805 ULTRA_TR=23", None, False)], "99% of the variant is repeat; 23 bp of tandem repeat"),
        ([("HG002.mat  1/1   HG002.pat  1/1", None, True)], "the deletion is in both assemblies"),
    ], "In tutorial/3_TSD_search/pangenome.vcf")
    b.append(panel)
    return svg(y2 + 16, b, "DEL.71_275: an AluYb8 present in the reference and absent from both HG002 haplotypes")


# --- Example 1: L1HS deletion ----------------------------------------------------
def fig_l1():
    b = [t(20, 30, "Example 1: L1HS deletion, HG002.pat.svim_asm.DEL.140_387 (chr22:11,244,543)", 16, weight="bold"),
         t(20, 52, "The 3' end of an L1HS, on the minus strand, in the reference. The paternal haplotype lacks it; the maternal one has it.", 13, fill=MUTED)]
    # full-length L1HS for scale, with the part present here
    b.append(t(170, 84, "Full-length L1HS consensus, 6,155 bp:", 12, fill=MUTED))
    b.append(seg(400, 900, 72, "plain", "5' UTR · ORF1 · ORF2 · 3' UTR", h=18, size=11))
    b.append(f'<rect x="892" y="72" width="8" height="18" fill="{C["l1"][0]}" stroke="{C["l1"][1]}" stroke-width="1.2"/>')
    b.append(t(896, 106, "positions 6,053–6,155: the part present here", 11.5, "end", fill=MUTED))
    ref = [lambda y: seg(170, 260, y, "flank", "…ACATGAA", fam=MONO),
           lambda y: seg(260, 360, y, "tsd", "TATTGCCTTTT", size=10.5, fam=MONO),
           lambda y: seg(360, 470, y, "polya", "T×42", size=11),
           lambda y: seg(470, 750, y, "l1", "L1HS 3' end  (− strand)", direction="-"),
           lambda y: seg(750, 850, y, "tsd", "TATTGCCTTTT", size=10.5, fam=MONO),
           lambda y: seg(850, 945, y, "flank", "AACAAAT…", fam=MONO)]
    hap_del = [lambda y: seg(170, 260, y, "flank", "…ACATGAA", fam=MONO), lambda y: gap(260, 750, y),
               lambda y: seg(750, 850, y, "tsd", "TATTGCCTTTT", size=10.5, fam=MONO),
               lambda y: seg(850, 945, y, "flank", "AACAAAT…", fam=MONO)]
    b.append(brace(260, 750, 136, "REF allele: 138 bp missing from the paternal haplotype"))
    body, y = rows_block(150, [("Reference (CHM13)", ref), ("HG002 maternal", ref), ("HG002 paternal", hap_del)])
    b += body
    b.append(t(170, y + 4, "On the minus strand, the poly(A) tail reads as a run of T on the + strand, at the left end.", 11.5, fill=MUTED))
    inset, y = graph_inset(20, y + 34, True, "l1", "L1HS", 2, 2, "1/0")
    b.append(inset)
    panel, y2 = vcf_panel(y + 6, [
        ([("SVLEN=-138", None, True)], "a deletion of 138 bp relative to the reference"),
        ([("repeat_ids=L1HS matching_classes=LINE/L1", C["l1"][1], True)], "RepeatMasker: one L1HS (n_hits=1)"),
        ([("RM_hit_strands=C match_lengths=130", C["l1"][1], True)], "C = complementary (minus) strand; 130 bp match"),
        ([("TSD=TATTGCCTTTT,TATTGCCTTTT", C["tsd"][1], True)], "the same 11 bp on both sides"),
        ([("polyA=TRUE", C["polya"][1], True)], "the tail, here the run of T"),
        ([("ULTRA_TR=47 ULTRA_TR_span=0.338129", None, False)], "47 bp of tandem repeat, mostly the run of T"),
        ([("GT:DP:AD  1/0:4:2,2", None, True)], "4 reads: 2 on each path, so heterozygous"),
    ], "In tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz")
    b.append(panel)
    return svg(y2 + 16, b, "Example 1: an L1HS 3' end deleted from the paternal HG002 haplotype")


# --- Example 2: AluYa5 insertion -------------------------------------------------
def fig_alu():
    b = [t(20, 30, "Example 2: AluYa5 insertion, HG002.mat.svim_asm.INS.182_786 (chr22:25,461,790)", 16, weight="bold"),
         t(20, 52, "An AluYa5 on the minus strand, absent from the reference and present in both HG002 assemblies, inserted in a run of T.", 13, fill=MUTED)]
    ref = [lambda y: seg(170, 300, y, "flank", "…CCTGTTTTTTTTTTT", size=10.5, fam=MONO),
           lambda y: gap(300, 680, y),
           lambda y: seg(680, 800, y, "tsd", "T×18", size=11),
           lambda y: seg(800, 945, y, "flank", "AAATGTAG…", fam=MONO)]
    hap = [lambda y: seg(170, 300, y, "flank", "…CCTGTTTTTTTTTTT", size=10.5, fam=MONO),
           lambda y: seg(300, 380, y, "polya", "T×17", size=11),
           lambda y: seg(380, 680, y, "alu", "AluYa5  (− strand)", direction="-"),
           lambda y: seg(680, 800, y, "tsd", "T×18", size=11),
           lambda y: seg(800, 945, y, "flank", "AAATGTAG…", fam=MONO)]
    b.append(brace(300, 680, 96, "ALT allele: 297 bp present in HG002, absent from the reference"))
    body, y = rows_block(110, [("Reference (CHM13)", ref), ("HG002 maternal", hap), ("HG002 paternal", hap)])
    b += body
    b.append(t(170, y + 4, "On the + strand, the Alu's tail is a run of T, and it sits next to a run of T in the reference.", 11.5, fill=MUTED))
    b.append(t(170, y + 20, "The TSD that GraffiTE found (18 T) overlaps the tail.", 11.5, fill=MUTED))
    inset, y = graph_inset(20, y + 50, False, "alu", "AluYa5", 1, 1, "0/1")
    b.append(inset)
    panel, y2 = vcf_panel(y + 6, [
        ([("SVLEN=297", None, True)], "positive: an insertion relative to the reference"),
        ([("repeat_ids=AluYa5 matching_classes=SINE/Alu", C["alu"][1], True)], "RepeatMasker: one AluYa5 (n_hits=1)"),
        ([("RM_hit_strands=C match_lengths=296", C["alu"][1], True)], "on the minus strand; 296 bp match"),
        ([("TSD=TTTTTTTTTTTTTTTTTT,TTTTTTTTTTTTTTTTTT", C["tsd"][1], True)], "18 T on both sides"),
        ([("polyA=TRUE", C["polya"][1], True)], "the tail, here the run of T"),
        ([("GT:DP:AD  0/1:2:1,1", None, True)], "only 2 reads, one per path: heterozygous call"),
        ([("SVs.vcf: HG002.mat 1/1  HG002.pat 1/1", None, False)], "both assemblies carry the insertion"),
    ], "In tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz (last line: tutorial/1_SV_search/SVs.vcf)")
    b.append(panel)
    return svg(y2 + 16, b, "Example 2: an AluYa5 insertion present in both HG002 assemblies")


# --- Example 3: not an MEI -------------------------------------------------------
def fig_sat():
    b = [t(20, 30, "Example 3: not an MEI, HG002.mat.svim_asm.DEL.83_330 (chr22:8,522,276)", 16, weight="bold"),
         t(20, 52, "341 bp of alpha satellite (ALR/Alpha) missing from the maternal haplotype, inside a satellite array.", 13, fill=MUTED)]
    units_ref = [(170, 300), (300, 430), (430, 560), (560, 690), (690, 820), (820, 945)]
    ref = [lambda y, a=a, z=z, i=i: seg(a, z, y, "sat", "ALR/Alpha" if i in (2, 3) else "") for i, (a, z) in enumerate(units_ref)]
    hap_del = [lambda y: seg(170, 300, y, "sat"), lambda y: seg(300, 430, y, "sat"), lambda y: gap(430, 690, y),
               lambda y: seg(690, 820, y, "sat"), lambda y: seg(820, 945, y, "sat")]
    b.append(brace(430, 690, 96, "REF allele: 341 bp, about two units of the array"))
    body, y = rows_block(110, [("Reference (CHM13)", ref), ("HG002 maternal", hap_del), ("HG002 paternal", ref)])
    b += body
    b.append(t(170, y + 4, "The same 341 bp occur 3.4 kb upstream and 3.4 kb downstream at 99–100% identity (blastn on the reference).", 11.5, fill=MUTED))
    panel, y2 = vcf_panel(y + 22, [
        ([("SVLEN=-341.0", None, True)], "a deletion of 341 bp relative to the reference"),
        ([("repeat_ids=ALR/Alpha matching_classes=Satellite/centr", C["sat"][1], True)], "RepeatMasker: centromeric alpha satellite, not a TE"),
        ([("RM_hit_strands=C match_lengths=339", C["sat"][1], True)], "339 of the 341 bp match"),
        ([("TSD=CAATCAA,CAATCAA", C["tsd"][1], True)], "the same 7 bp at both ends"),
        ([("polyA=FALSE", C["polya"][1], True)], "no poly(A) tail"),
        ([("HG002.mat  1/1   HG002.pat  0/0", None, True)], "missing from the maternal assembly only"),
    ], "In tutorial/3_TSD_search/pangenome.vcf")
    b.append(panel)
    return svg(y2 + 16, b, "Example 3: a deletion in an alpha-satellite array, not a mobile element insertion")


# --- Background figure 2: how short reads find a TE insertion ---------------------
def mate(x, y, d, col, L=38):
    """One read of a pair, as a short arrow pointing in its direction."""
    if d == "+":
        return (f'<line x1="{x}" y1="{y}" x2="{x + L - 6}" y2="{y}" stroke="{col}" stroke-width="3"/>'
                f'<path d="M{x + L - 7},{y - 5} L{x + L},{y} L{x + L - 7},{y + 5} z" fill="{col}"/>')
    return (f'<line x1="{x + 6}" y1="{y}" x2="{x + L}" y2="{y}" stroke="{col}" stroke-width="3"/>'
            f'<path d="M{x + 7},{y - 5} L{x},{y} L{x + 7},{y + 5} z" fill="{col}"/>')


def link(x1, y1, x2, y2):
    """Dotted line joining the two reads of a pair."""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{MUTED}" stroke-width="1.2" '
            f'stroke-dasharray="2 3"/>')


def fig_shortreads():
    te = C["alu"][1]
    b = [t(20, 30, "What short reads show at a TE insertion", 16, weight="bold")]
    for x, ti in [(20, "TE in the sample, not in the reference"), (490, "TE in the reference, not in the sample")]:
        b.append(f'<rect x="{x}" y="50" width="450" height="300" rx="10" fill="#ffffff" stroke="#d0d7de"/>')
        b.append(t(x + 16, 78, ti, 14, weight="bold"))
    # left: reads that point to a TE absent from the reference
    x, y = 20, 230
    xi = x + 225
    b.append(t(x + 16, 124, "TE library", 12, fill=MUTED))
    b.append(seg(x + 110, x + 340, 110, "alu", "TE consensus", h=18, size=11))
    b.append(t(x + 16, y + 13, "reference", 12, fill=MUTED))
    b.append(seg(x + 90, x + 430, y, "flank", h=18))
    b.append(f'<line x1="{xi}" y1="{y - 50}" x2="{xi}" y2="{y + 30}" stroke="#b5462f" stroke-width="1.5" stroke-dasharray="4 3"/>')
    b.append(t(xi, y + 44, "insertion point", 11.5, "middle", fill="#b5462f"))
    b.append(mate(x + 150, y - 10, "+", TEXT) + mate(x + 160, 102, "-", te) + link(x + 169, y - 14, x + 179, 106))
    b.append(mate(x + 262, y - 10, "-", TEXT) + mate(x + 250, 102, "+", te) + link(x + 281, y - 14, x + 269, 106))
    b.append(f'<line x1="{xi - 40}" y1="{y - 34}" x2="{xi}" y2="{y - 34}" stroke="{TEXT}" stroke-width="3"/>'
             f'<line x1="{xi}" y1="{y - 34}" x2="{xi + 34}" y2="{y - 58}" stroke="{te}" stroke-width="3"/>')
    b.append(t(xi + 40, y - 56, "split read", 11.5, fill=MUTED))
    b.append(t(x + 16, 310, "One read maps next to the insertion point and its mate matches", 12, fill=MUTED))
    b.append(t(x + 16, 328, "the TE library, or one read is split between flank and TE.", 12, fill=MUTED))
    # right: reads that skip a TE the sample lacks
    x = 490
    b.append(t(x + 16, y + 13, "reference", 12, fill=MUTED))
    b.append(seg(x + 90, x + 190, y, "flank", h=18))
    b.append(seg(x + 190, x + 320, y, "alu", "TE", h=18, size=11))
    b.append(seg(x + 320, x + 430, y, "flank", h=18))
    for i, (a, z) in enumerate([(x + 140, x + 330), (x + 110, x + 360)]):
        yy = y - 10 - 16 * i
        b.append(mate(a, yy, "+", TEXT) + mate(z, yy, "-", TEXT))
        b.append(f'<path d="M{a + 38},{yy} C{a + 70},{yy - 45 - 15 * i} {z - 32},{yy - 45 - 15 * i} {z},{yy}" '
                 f'fill="none" stroke="{MUTED}" stroke-width="1.2" stroke-dasharray="2 3"/>')
    b.append(brace(x + 148, x + 368, 130, "mates farther apart than the insert size"))
    b.append(f'<line x1="{x + 150}" y1="{y + 34}" x2="{x + 190}" y2="{y + 34}" stroke="{TEXT}" stroke-width="3"/>'
             f'<line x1="{x + 190}" y1="{y + 34}" x2="{x + 320}" y2="{y + 34}" stroke="{MUTED}" stroke-width="1.2" stroke-dasharray="5 4"/>'
             f'<line x1="{x + 320}" y1="{y + 34}" x2="{x + 360}" y2="{y + 34}" stroke="{TEXT}" stroke-width="3"/>')
    b.append(t(x + 255, y + 52, "split read: jumps over the TE", 11.5, "middle", fill=MUTED))
    b.append(t(x + 16, 310, "The two mates of a pair land on either side of the TE, too far", 12, fill=MUTED))
    b.append(t(x + 16, 328, "apart, or one read is split across it.", 12, fill=MUTED))
    return svg(370, b, "What short reads show at a TE insertion: discordant read pairs and split reads")


# --- Overview: the four steps ----------------------------------------------------
def fig_overview():
    b = [t(20, 30, "What GraffiTE does, in four steps", 16, weight="bold")]
    px = [(20, 60), (490, 60), (20, 330), (490, 330)]
    titles = ["1. Compare each assembly to the reference", "2. Keep the differences made of repeats",
              "3. Look for the signs of a transposition", "4. Genotype with short reads on a graph"]
    for (x, y), ti in zip(px, titles):
        b.append(f'<rect x="{x}" y="{y}" width="450" height="250" rx="10" fill="#ffffff" stroke="#d0d7de"/>')
        b.append(t(x + 16, y + 28, ti, 14, weight="bold"))
    # 1: a deletion and an insertion
    x, y = 20, 60
    b.append(t(x + 16, y + 74, "reference", 12, fill=MUTED))
    b.append(seg(x + 100, x + 190, y + 60, "flank", h=20))
    b.append(seg(x + 190, x + 260, y + 60, "plain", "DEL", h=20, size=11))
    b.append(seg(x + 260, x + 430, y + 60, "flank", h=20))
    b.append(f'<path d="M{x + 340},{y + 84} l-6,9 l12,0 z" fill="{MUTED}"/>')
    b.append(t(x + 16, y + 144, "assembly", 12, fill=MUTED))
    b.append(seg(x + 100, x + 190, y + 130, "flank", h=20))
    b.append(gap(x + 190, x + 260, y + 130, h=20))
    b.append(seg(x + 260, x + 340, y + 130, "flank", h=20))
    b.append(seg(x + 340, x + 400, y + 130, "plain", "INS", h=20, size=11))
    b.append(seg(x + 400, x + 430, y + 130, "flank", h=20))
    b.append(t(x + 225, y + 172, "DEL: only in the reference", 11.5, "middle", fill=MUTED))
    b.append(t(x + 370, y + 172, "INS: only in the assembly", 11.5, "middle", fill=MUTED))
    b.append(t(x + 16, y + 225, "Each insertion (INS) and deletion (DEL) is a candidate.", 12, fill=MUTED))
    # 2: which differences are repeats
    x, y = 490, 60
    for i, (kind, lab, verdict) in enumerate([("alu", "TE", "kept"), ("sat", "satellite", "kept (filtered later)"), ("plain", "unique sequence", "dropped")]):
        yy = y + 56 + i * 50
        b.append(seg(x + 20, x + 200, yy, kind, lab, h=26, size=12))
        b.append(t(x + 220, yy + 18, verdict, 13, weight="bold" if verdict == "kept" else "normal",
                   fill=TEXT if verdict != "dropped" else MUTED))
    b.append(t(x + 16, y + 225, "RepeatMasker labels each candidate with the TE library.", 12, fill=MUTED))
    # 3: hallmarks
    x, y = 20, 330
    b.append(seg(x + 30, x + 80, y + 90, "flank", h=26))
    b.append(seg(x + 80, x + 140, y + 90, "tsd", "TSD", h=26, size=11))
    b.append(seg(x + 140, x + 300, y + 90, "alu", "TE", direction="+", h=26))
    b.append(seg(x + 300, x + 350, y + 90, "polya", "AAAA", h=26, size=11))
    b.append(seg(x + 350, x + 410, y + 90, "tsd", "TSD", h=26, size=11))
    b.append(seg(x + 410, x + 430, y + 90, "flank", h=26))
    b.append(t(x + 110, y + 140, "same short sequence", 11.5, "middle", fill=MUTED))
    b.append(t(x + 380, y + 140, "on both sides", 11.5, "middle", fill=MUTED))
    b.append(t(x + 325, y + 80, "poly(A) tail", 11.5, "middle", fill=MUTED))
    b.append(t(x + 16, y + 200, "A target site duplication (TSD) and, for Alu, L1 and SVA,", 12, fill=MUTED))
    b.append(t(x + 16, y + 218, "a poly(A) tail are left behind by the insertion.", 12, fill=MUTED))
    # 4: graph and reads
    x, y = 490, 330
    yb = y + 110
    fill, stroke = C["alu"]
    b.append(seg(x + 20, x + 100, yb, "flank", h=22))
    b.append(seg(x + 350, x + 430, yb, "flank", h=22))
    b.append(f'<path d="M{x + 100},{yb + 11} C{x + 130},{yb - 30} {x + 160},{yb - 30} {x + 185},{yb - 30}" fill="none" stroke="{stroke}" stroke-width="2"/>')
    b.append(seg(x + 185, x + 265, yb - 41, "alu", "TE", h=22, size=11))
    b.append(f'<path d="M{x + 265},{yb - 30} C{x + 290},{yb - 30} {x + 320},{yb - 30} {x + 350},{yb + 11}" fill="none" stroke="{stroke}" stroke-width="2"/>')
    b.append(f'<path d="M{x + 100},{yb + 11} C{x + 150},{yb + 55} {x + 300},{yb + 55} {x + 350},{yb + 11}" fill="none" stroke="{MUTED}" stroke-width="2"/>')
    for i in range(2):
        b.append(f'<line x1="{x + 190 + 10 * i}" y1="{yb - 50 - 7 * i}" x2="{x + 250 + 10 * i}" y2="{yb - 50 - 7 * i}" stroke="#24292f" stroke-width="2.5"/>')
        b.append(f'<line x1="{x + 190 + 10 * i}" y1="{yb + 50 + 7 * i}" x2="{x + 250 + 10 * i}" y2="{yb + 50 + 7 * i}" stroke="#24292f" stroke-width="2.5"/>')
    b.append(t(x + 275, yb - 50, "reads with the TE", 11.5, fill=MUTED))
    b.append(t(x + 275, yb + 58, "reads without", 11.5, fill=MUTED))
    b.append(t(x + 16, y + 225, "2 reads on each path: one copy with the TE, one without (0/1).", 12, fill=MUTED))
    return svg(600, b, "GraffiTE in four steps: compare, keep the repeats, check the hallmarks, genotype")


if __name__ == "__main__":
    for name, f in [("overview_steps.svg", fig_overview), ("variant_DEL71_AluYb8.svg", fig_del71),
                    ("example1_L1HS.svg", fig_l1), ("example2_AluYa5.svg", fig_alu),
                    ("example3_satellite.svg", fig_sat),
                    ("background_shortreads.svg", fig_shortreads)]:
        (OUT / name).write_text(f())
        print("wrote", OUT / name)
