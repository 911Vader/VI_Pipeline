"""
pdf_report.py — Generate a PDF report using reportlab Platypus.

Report structure:
  1. Title + sample/case identification
  2. Mandatory disclaimer (bold, prominent)
  3. Summary: variant counts by ACMG category
  4. Top-N candidate variants table
  5. VUS in listed genes section (only when such variants exist)
  6. Footer with generation timestamp

The mandatory disclaimer must appear verbatim:
  "Top genomic candidate(s) requiring expert review.
   Variant pathogenicity alone does not prove cause of death."

max_candidates_in_report is read from settings.reporting.max_candidates_in_report
(default 20).
"""

import datetime
import logging
from pathlib import Path
from typing import List, Optional

from analysis.parser import Variant

logger = logging.getLogger(__name__)

#: Mandatory disclaimer — must not be altered.
MANDATORY_DISCLAIMER = (
    "Top genomic candidate(s) requiring expert review. "
    "Variant pathogenicity alone does not prove cause of death."
)

PIPELINE_VERSION = "1.0.0"


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_pdf(
    variants: List[Variant],
    output_path: str,
    metadata: dict,
    settings: Optional[dict] = None,
) -> str:
    """
    Generate a PDF candidate report.

    Parameters
    ----------
    variants:    Ranked list from ranker.rank_variants().
    output_path: Destination .pdf path.
    metadata:    Dict with sample_id, case_id, etc.
    settings:    Parsed settings.yaml dict. Controls max_candidates_in_report.

    Returns
    -------
    Path to the written PDF file (str).
    """
    try:
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.units import cm
        from reportlab.lib import colors
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            HRFlowable,
        )
        from reportlab.lib.enums import TA_CENTER, TA_LEFT
    except ImportError:
        raise RuntimeError(
            "reportlab is required for PDF output.\n"
            "  Install it:  pip install reportlab"
        )

    settings    = settings or {}
    max_cands   = int(
        settings.get("reporting", {}).get("max_candidates_in_report", 20)
    )
    sample_id   = metadata.get("sample_id", "UNKNOWN")
    case_id     = metadata.get("case_id", "")
    generated   = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    retained = [v for v in variants if not v.excluded]
    candidates = retained[:max_cands]

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(out),
        pagesize=landscape(A4),
        leftMargin=1.5 * cm, rightMargin=1.5 * cm,
        topMargin=2 * cm,    bottomMargin=2 * cm,
        title=f"VI_Pipeline Report — {sample_id}",
        author="VI_Pipeline",
    )

    styles = getSampleStyleSheet()
    story  = []

    # --- Title ---
    title_style = ParagraphStyle(
        "Title2",
        parent=styles["Title"],
        fontSize=16,
        spaceAfter=6,
    )
    story.append(Paragraph("VI_Pipeline — Variant Prioritisation Report", title_style))

    # --- Sample / case info ---
    sub_style = ParagraphStyle(
        "Sub", parent=styles["Normal"], fontSize=10, spaceAfter=2,
    )
    story.append(Paragraph(f"<b>Sample ID:</b> {sample_id}", sub_style))
    if case_id:
        story.append(Paragraph(f"<b>Case ID:</b> {case_id}", sub_style))
    input_vcf = metadata.get("input_vcf", "")
    if input_vcf:
        story.append(Paragraph(f"<b>Input:</b> {input_vcf}", sub_style))
    story.append(Paragraph(f"<b>Generated:</b> {generated}", sub_style))
    story.append(Spacer(1, 0.3 * cm))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.grey))
    story.append(Spacer(1, 0.3 * cm))

    # --- Mandatory disclaimer ---
    disclaimer_style = ParagraphStyle(
        "Disclaimer",
        parent=styles["Normal"],
        fontSize=11,
        textColor=colors.red,
        borderColor=colors.red,
        borderWidth=1,
        borderPadding=6,
        backColor=colors.Color(1, 0.95, 0.95),
        spaceAfter=8,
        spaceBefore=4,
        leading=16,
    )
    story.append(Paragraph(f"⚠ {MANDATORY_DISCLAIMER}", disclaimer_style))
    story.append(Spacer(1, 0.4 * cm))

    # --- Summary table ---
    story.append(Paragraph("<b>Summary</b>", styles["Heading2"]))
    story.append(_build_summary_table(variants, colors, Table, TableStyle))
    story.append(Spacer(1, 0.5 * cm))

    # --- Top candidates table ---
    n_label = f"Top {len(candidates)}" if len(candidates) < len(retained) else "All"
    story.append(
        Paragraph(
            f"<b>{n_label} Retained Candidates</b> "
            f"({len(retained)} total retained, {len(variants) - len(retained)} excluded)",
            styles["Heading2"],
        )
    )

    if candidates:
        story.append(
            _build_candidates_table(candidates, styles, colors, Table, TableStyle, Paragraph)
        )
    else:
        story.append(Paragraph("No retained variants to display.", styles["Normal"]))

    # --- VUS in listed genes ---
    vus_gene_list = [
        v for v in variants
        if not v.excluded
        and v.acmg_classification == "VUS"
        and getattr(v, "in_gene_list", False)
        and v not in candidates
    ]
    if vus_gene_list:
        story.append(Spacer(1, 0.5 * cm))
        story.append(
            Paragraph(
                "<b>Variants of Uncertain Significance in Listed Genes</b>",
                styles["Heading2"],
            )
        )
        vus_note_style = ParagraphStyle(
            "VusNote", parent=styles["Normal"], fontSize=9,
            textColor=colors.HexColor("#7B4A00"), spaceAfter=4,
        )
        story.append(
            Paragraph(
                "These variants are classified as VUS under ACMG/AMP 2015 criteria — "
                "neither pathogenic nor benign evidence is sufficient for reclassification. "
                "They are listed here because the gene is present in the supplied gene list "
                "and may be clinically relevant. Expert review is required.",
                vus_note_style,
            )
        )
        story.append(
            _build_candidates_table(vus_gene_list, styles, colors, Table, TableStyle, Paragraph)
        )

    story.append(Spacer(1, 0.5 * cm))
    story.append(HRFlowable(width="100%", thickness=0.5, color=colors.lightgrey))

    # --- Footer ---
    footer_style = ParagraphStyle(
        "Footer", parent=styles["Normal"], fontSize=8,
        textColor=colors.grey, spaceAfter=0,
    )
    story.append(
        Paragraph(
            f"VI_Pipeline {PIPELINE_VERSION} | GRCh38 | Single-sample SNV/indel | "
            f"ACMG/AMP classification | This report requires expert clinical review.",
            footer_style,
        )
    )

    doc.build(story)
    logger.info("PDF report written: %s (%d candidates)", out, len(candidates))
    return str(out)


# ---------------------------------------------------------------------------
# Table builders
# ---------------------------------------------------------------------------

def _build_summary_table(variants, colors, Table, TableStyle):
    """Build a compact summary table of variant counts by ACMG class."""
    from collections import Counter
    from reportlab.lib.units import cm

    retained = [v for v in variants if not v.excluded]
    excluded = [v for v in variants if v.excluded]
    counts   = Counter(v.acmg_classification for v in retained if v.acmg_classification)

    rows = [
        ["Category", "Count"],
        ["Total variants",          str(len(variants))],
        ["Retained (passed filter)", str(len(retained))],
        ["Excluded (filtered out)", str(len(excluded))],
        ["Pathogenic",              str(counts.get("Pathogenic", 0))],
        ["Likely Pathogenic",       str(counts.get("Likely Pathogenic", 0))],
        ["VUS",                     str(counts.get("VUS", 0))],
        ["Likely Benign",           str(counts.get("Likely Benign", 0))],
        ["Benign",                  str(counts.get("Benign", 0))],
    ]

    tbl = Table(rows, colWidths=[8 * cm, 4 * cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND",  (0, 0), (-1, 0), colors.HexColor("#4472C4")),
        ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
        ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",    (0, 0), (-1, -1), 9),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.Color(0.96, 0.96, 0.96)]),
        ("GRID",        (0, 0), (-1, -1), 0.5, colors.lightgrey),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
    ]))
    return tbl


# Column spec for the candidates table: (header, accessor, width_cm)
# Total target: <= 26.7 cm (A4 landscape minus 1.5 cm margins each side).
# Scores are shown as "-" when not available (plugin not configured or no match).
_CAND_COLS = [
    ("Gene",        lambda v, i: v.gene_symbol or "",                                               2.0),
    ("Chr:Pos",     lambda v, i: f"{v.chrom}:{v.pos}",                                             2.8),
    ("HGVS c.",     lambda v, i: _trunc(v.hgvs_c, 18),                                             3.0),
    ("HGVS p.",     lambda v, i: _trunc(v.hgvs_p, 14),                                             2.5),
    ("Consequence", lambda v, i: _trunc(v.consequence, 16),                                         2.5),
    ("ACMG",        lambda v, i: v.acmg_classification or "",                                       2.5),
    ("gnomAD AF",   lambda v, i: f"{v.gnomad_af:.4g}" if v.gnomad_af is not None else "-",         1.5),
    ("ClinVar",     lambda v, i: _trunc(v.clinvar_significance, 14),                                2.2),
    ("CADD",        lambda v, i: f"{v.cadd_phred:.1f}" if v.cadd_phred is not None else "-",       1.5),
    ("REVEL",       lambda v, i: f"{v.revel_score:.3f}" if v.revel_score is not None else "-",     1.3),
    ("SpliceAI",    lambda v, i: f"{v.spliceai_max_score:.3f}" if v.spliceai_max_score is not None else "-", 1.5),
    ("LOFTEE",      lambda v, i: v.loftee_lof or "-",                                              1.3),
    ("ClinGen",     lambda v, i: _trunc(v.clingen_validity, 12) if v.clingen_validity else "-",   1.8),
]
# Total width: 2.0+2.8+3.0+2.5+2.5+2.5+1.5+2.2+1.5+1.3+1.5+1.3+1.8 = 26.4 cm

_ACMG_PDF_COLOURS = {
    "Pathogenic":        (1.0, 0.8, 0.8),
    "Likely Pathogenic": (1.0, 0.88, 0.7),
    "VUS":               (1.0, 0.99, 0.87),
    "Likely Benign":     (0.89, 0.95, 1.0),
    "Benign":            (0.91, 0.97, 0.91),
}


def _build_candidates_table(candidates, styles, colors, Table, TableStyle, Paragraph):
    """Build the top-candidates table for the PDF."""
    from reportlab.lib.units import cm
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_LEFT

    # Cell styles — created fresh (never mutate the shared styles dict)
    cell_style = ParagraphStyle(
        "CandCell",
        parent=styles["Normal"],
        fontSize=7.5,
        leading=9,
        wordWrap="CJK",   # wraps at any char when no spaces exist (HGVS etc.)
        alignment=TA_LEFT,
    )
    hdr_style = ParagraphStyle(
        "CandHdr",
        parent=styles["Normal"],
        fontSize=7.5,
        leading=9,
        fontName="Helvetica-Bold",
        textColor=colors.white,
        wordWrap="CJK",
        alignment=TA_LEFT,
    )

    widths = [col[2] * cm for col in _CAND_COLS]

    # Wrap every cell in Paragraph so text reflows within column width
    header_row = [Paragraph(col[0], hdr_style) for col in _CAND_COLS]
    data = [header_row]
    for i, v in enumerate(candidates):
        row = [Paragraph(str(col[1](v, i)), cell_style) for col in _CAND_COLS]
        data.append(row)

    tbl = Table(data, colWidths=widths, repeatRows=1)

    style_cmds = [
        ("BACKGROUND",    (0, 0), (-1, 0),  colors.HexColor("#4472C4")),
        ("VALIGN",        (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING",   (0, 0), (-1, -1), 4),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 4),
        ("TOPPADDING",    (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("GRID",          (0, 0), (-1, -1), 0.4, colors.lightgrey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.Color(0.97, 0.97, 0.97)]),
    ]

    # Per-row ACMG colour on the ACMG column (index 5)
    acmg_col_idx = 5
    for row_idx, v in enumerate(candidates, start=1):
        rgb = _ACMG_PDF_COLOURS.get(v.acmg_classification or "")
        if rgb:
            style_cmds.append(
                ("BACKGROUND", (acmg_col_idx, row_idx), (acmg_col_idx, row_idx),
                 colors.Color(*rgb))
            )

    tbl.setStyle(TableStyle(style_cmds))
    return tbl


def _trunc(value: Optional[str], max_len: int) -> str:
    """Truncate a string to max_len chars, adding '…' if cut."""
    if not value:
        return ""
    if len(value) <= max_len:
        return value
    return value[: max_len - 1] + "…"
