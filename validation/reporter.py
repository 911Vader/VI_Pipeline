"""
reporter.py — Generate all validation output files.

Produces four artefacts per validation run:
  1. validation_summary.pdf     — metrics, pass/fail, run metadata
  2. discordance_log.xlsx        — one row per discordant variant
  3. tool_concordance.xlsx       — pairwise concordance matrix
  4. raw_results.tsv             — full results table
  5. raw_results.json            — same data, JSON format
"""

import datetime
import json
import logging
from pathlib import Path
from typing import List, Optional

from validation.concordance import (
    VariantResult, ValidationMetrics,
    CONCORDANT, BORDERLINE, DIRECTION_ERROR, CRITICAL_FAILURE, NOT_FOUND,
)
from validation.discordance_log import DiscordanceEntry
from validation.tool_concordance import ToolConcordanceReport, TOOLS

logger = logging.getLogger(__name__)

PIPELINE_VERSION = "1.0.0"

# ---------------------------------------------------------------------------
# PDF summary
# ---------------------------------------------------------------------------

def write_pdf_summary(
    output_path: str,
    metrics: ValidationMetrics,
    tool_concordance: ToolConcordanceReport,
    metadata: dict,
) -> str:
    """Generate the PDF validation summary report."""
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import cm
        from reportlab.lib import colors
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            HRFlowable,
        )
        from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    except ImportError:
        raise RuntimeError("reportlab is required: pip install reportlab")

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(out), pagesize=A4,
        leftMargin=2*cm, rightMargin=2*cm,
        topMargin=2*cm, bottomMargin=2*cm,
        title="VI_Pipeline Validation Summary",
        author="VI_Pipeline",
    )

    styles  = getSampleStyleSheet()
    story   = []
    ts_now  = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    title_style = ParagraphStyle(
        "VTitle", parent=styles["Title"], fontSize=18, spaceAfter=8,
    )
    h2 = ParagraphStyle(
        "VH2", parent=styles["Heading2"], fontSize=13, spaceAfter=4, spaceBefore=12,
    )
    body = ParagraphStyle(
        "VBody", parent=styles["Normal"], fontSize=10, spaceAfter=3,
    )

    verdict, verdict_reasons = metrics.pass_fail_verdict()
    verdict_colour = {
        "PASS": colors.HexColor("#2e7d32"),
        "CONDITIONAL PASS": colors.HexColor("#e65100"),
        "FAIL": colors.HexColor("#c62828"),
    }.get(verdict, colors.black)

    verdict_style = ParagraphStyle(
        "Verdict", parent=styles["Title"],
        fontSize=22, textColor=verdict_colour, spaceAfter=6,
    )

    # --- Header ---
    story.append(Paragraph("VI_Pipeline — Internal Validation Report", title_style))
    story.append(Paragraph(f"<b>Pipeline version:</b> {PIPELINE_VERSION}", body))
    story.append(Paragraph(f"<b>Run date:</b> {ts_now}", body))
    story.append(Paragraph(
        f"<b>Sample/run:</b> {metadata.get('sample_id', 'N/A')}  |  "
        f"<b>Truth set:</b> {metadata.get('truth_set_path', 'N/A')}",
        body,
    ))
    story.append(HRFlowable(width="100%", thickness=1, spaceAfter=12))

    # --- Verdict banner ---
    story.append(Paragraph(f"Verdict: {verdict}", verdict_style))
    for r in verdict_reasons:
        story.append(Paragraph(f"• {r}", body))
    story.append(Spacer(1, 0.4*cm))

    # --- Classification performance ---
    story.append(Paragraph("1. Classification Performance", h2))

    def pct(v: Optional[float]) -> str:
        return f"{v*100:.1f}%" if v is not None else "N/A"

    perf_data = [
        ["Metric", "Value"],
        ["Truth set variants",           str(metrics.n_truth_total)],
        ["Matched in VCF",               str(metrics.n_matched)],
        ["Not found in VCF",             str(metrics.n_not_found)],
        ["Concordant",                   f"{metrics.n_concordant}  ({pct(metrics.concordance_rate)})"],
        ["Borderline (VUS vs LP/LB)",    str(metrics.n_borderline)],
        ["Direction errors",             str(metrics.n_direction_error)],
        ["Critical failures (P↔B)",      str(metrics.n_critical_failure)],
        ["Sensitivity (P/LP detection)", pct(metrics.sensitivity)],
        ["Specificity (B/LB detection)", pct(metrics.specificity)],
        ["Pipeline VUS rate",            pct(metrics.vus_rate)],
        ["Critical failure rate",        pct(metrics.critical_failure_rate)],
    ]
    story.append(_make_table(perf_data, col_widths=[10*cm, 6*cm]))
    story.append(Spacer(1, 0.4*cm))

    # --- Tool concordance ---
    story.append(Paragraph("2. Tool Concordance Summary", h2))
    tc_data = [["Tool Pair", "Variants Both Called", "Concordant", "Rate"]]
    for stat in tool_concordance.pairwise:
        tc_data.append([
            stat.pair_label,
            str(stat.n_both_called),
            str(stat.n_concordant),
            pct(stat.concordance_rate),
        ])
    story.append(_make_table(tc_data, col_widths=[7*cm, 4*cm, 3*cm, 3*cm]))
    story.append(Spacer(1, 0.3*cm))
    story.append(Paragraph(
        f"Variants with majority in-silico discordance: "
        f"<b>{tool_concordance.n_majority_discordant}</b>",
        body,
    ))
    story.append(Paragraph(
        f"ClinVar vs pipeline direction conflicts flagged: "
        f"<b>{tool_concordance.n_clinvar_conflicts}</b>",
        body,
    ))

    # --- Pass/fail per dimension ---
    story.append(Paragraph("3. Pass / Fail Per Validation Dimension", h2))
    dim_data = [["Dimension", "Status", "Notes"]]

    # Classification accuracy
    concordance_rate = metrics.concordance_rate or 0.0
    cls_status = "PASS" if concordance_rate >= 0.95 else (
        "CONDITIONAL PASS" if concordance_rate >= 0.90 else "FAIL"
    )
    dim_data.append(["Classification Accuracy", cls_status,
                     f"Concordance: {pct(metrics.concordance_rate)}"])

    # Critical failures
    cf_status = "PASS" if metrics.n_critical_failure == 0 else "FAIL"
    dim_data.append(["Critical Failures", cf_status,
                     f"{metrics.n_critical_failure} P↔B direction errors"])

    # Not-found rate
    nf_rate = (metrics.n_not_found / max(metrics.n_truth_total, 1))
    nf_status = "PASS" if nf_rate <= 0.01 else "FAIL"
    dim_data.append(["Workflow Integrity (variant loss)",
                     nf_status,
                     f"{metrics.n_not_found} variants not found in VCF ({pct(nf_rate)} loss)"])

    story.append(_make_table(dim_data, col_widths=[6.5*cm, 3.5*cm, 7*cm]))

    # --- Footer note ---
    story.append(Spacer(1, 0.6*cm))
    story.append(HRFlowable(width="100%", thickness=0.5))
    story.append(Paragraph(
        "This report is generated for internal quality control. "
        "No result may be manually overridden without documentation and second-reviewer sign-off. "
        "All discordances are detailed in the accompanying discordance_log.xlsx.",
        ParagraphStyle("Footer", parent=styles["Normal"], fontSize=8, textColor=colors.grey),
    ))

    doc.build(story)
    logger.info("Validation PDF written to %s", output_path)
    return output_path


def _make_table(data, col_widths):
    from reportlab.platypus import Table, TableStyle
    from reportlab.lib import colors

    t = Table(data, colWidths=col_widths, repeatRows=1)
    style = [
        ("BACKGROUND",  (0, 0), (-1, 0), colors.HexColor("#37474f")),
        ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
        ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",    (0, 0), (-1, -1), 9),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f5f5")]),
        ("GRID",        (0, 0), (-1, -1), 0.4, colors.HexColor("#bdbdbd")),
        ("VALIGN",      (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING",  (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 4),
    ]
    t.setStyle(TableStyle(style))
    return t


# ---------------------------------------------------------------------------
# Discordance log XLSX
# ---------------------------------------------------------------------------

_DISCORD_COLUMNS = [
    ("VARIANT",             lambda e: e.variant_hgvs),
    ("GENE",                lambda e: e.gene),
    ("VARIANT ID",          lambda e: e.variant_id),
    ("EXPECTED",            lambda e: e.expected),
    ("PIPELINE OUTPUT",     lambda e: e.pipeline_output),
    ("CRITERIA APPLIED",    lambda e: e.criteria_applied),
    ("EXPECTED CRITERIA",   lambda e: e.expected_acmg_codes),
    ("gnomAD AF (global)",  lambda e: e.gnomad_af_global),
    ("ClinVar Sig.",        lambda e: e.clinvar_significance),
    ("ClinVar Review",      lambda e: e.clinvar_review_status),
    ("CADD phred",          lambda e: e.cadd),
    ("REVEL",               lambda e: e.revel),
    ("SpliceAI",            lambda e: e.spliceai),
    ("AlphaMissense",       lambda e: e.alphamissense),
    ("LOFTEE",              lambda e: e.loftee),
    ("CONFLICT TYPE",       lambda e: e.conflict_type),
    ("ROOT CAUSE",          lambda e: e.root_cause),
    ("ACTION REQUIRED",     lambda e: e.action_required),
    ("CURATOR NOTES",       lambda e: e.curator_notes),
]

_CONFLICT_FILL = {
    CRITICAL_FAILURE: "FFCCCC",
    DIRECTION_ERROR:  "FFE0B2",
    BORDERLINE:       "FFF9C4",
    NOT_FOUND:        "E3F2FD",
}


def write_discordance_xlsx(
    output_path: str,
    entries: List[DiscordanceEntry],
) -> str:
    """Write the discordance log XLSX."""
    import openpyxl
    from openpyxl.styles import PatternFill, Font, Alignment

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Discordance Log"

    # Header row
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="37474F")
    headers = [col[0] for col in _DISCORD_COLUMNS]
    for col_idx, hdr in enumerate(headers, start=1):
        cell = ws.cell(row=1, column=col_idx, value=hdr)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(wrap_text=True, horizontal="center")

    # Data rows
    for row_idx, entry in enumerate(entries, start=2):
        fill_hex = _CONFLICT_FILL.get(entry.conflict_type, "FFFFFF")
        row_fill = PatternFill("solid", fgColor=fill_hex)
        for col_idx, (_, accessor) in enumerate(_DISCORD_COLUMNS, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=accessor(entry))
            cell.fill = row_fill
            cell.alignment = Alignment(wrap_text=True)

    # Column widths
    col_widths = [
        20, 10, 28, 16, 16, 22, 22,
        14, 10, 10, 10, 10, 10, 10, 10, 10,
        30, 30, 50,
    ]
    assert len(col_widths) == len(_DISCORD_COLUMNS)
    for col_idx, width in enumerate(col_widths, start=1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(col_idx)].width = width

    ws.freeze_panes = "A2"

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_path)
    logger.info("Discordance log written to %s (%d entries)", output_path, len(entries))
    return output_path


# ---------------------------------------------------------------------------
# Tool concordance matrix XLSX
# ---------------------------------------------------------------------------

def write_tool_concordance_xlsx(
    output_path: str,
    report: ToolConcordanceReport,
) -> str:
    """Write the pairwise tool concordance matrix XLSX."""
    import openpyxl
    from openpyxl.styles import PatternFill, Font, Alignment

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Tool Concordance"

    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="37474F")

    # Build matrix dict
    matrix = report.matrix()

    # Write header row: blank + tool names
    ws.cell(row=1, column=1, value="").fill = header_fill
    for col_idx, tool in enumerate(TOOLS, start=2):
        c = ws.cell(row=1, column=col_idx, value=tool)
        c.font = header_font
        c.fill = header_fill
        c.alignment = Alignment(horizontal="center")

    # Write data rows
    for row_idx, tool_a in enumerate(TOOLS, start=2):
        c = ws.cell(row=row_idx, column=1, value=tool_a)
        c.font = Font(bold=True)
        c.fill = PatternFill("solid", fgColor="ECEFF1")
        for col_idx, tool_b in enumerate(TOOLS, start=2):
            if tool_a == tool_b:
                cell = ws.cell(row=row_idx, column=col_idx, value="—")
                cell.fill = PatternFill("solid", fgColor="ECEFF1")
                continue
            # Concordance is symmetric — look up in either order
            rate = matrix.get((tool_a, tool_b)) or matrix.get((tool_b, tool_a))
            if rate is None:
                val = "N/A"
                fill_hex = "F5F5F5"
            else:
                val = f"{rate*100:.1f}%"
                # Green ≥90%, yellow 70-90%, red <70%
                fill_hex = (
                    "C8E6C9" if rate >= 0.90 else
                    "FFF9C4" if rate >= 0.70 else
                    "FFCDD2"
                )
            cell = ws.cell(row=row_idx, column=col_idx, value=val)
            cell.fill = PatternFill("solid", fgColor=fill_hex)
            cell.alignment = Alignment(horizontal="center")

    # Summary rows
    n_row = len(TOOLS) + 3
    ws.cell(row=n_row, column=1, value="Majority in-silico discordant variants:")
    ws.cell(row=n_row, column=2, value=report.n_majority_discordant)
    ws.cell(row=n_row + 1, column=1, value="ClinVar vs pipeline direction conflicts:")
    ws.cell(row=n_row + 1, column=2, value=report.n_clinvar_conflicts)

    # Pairwise detail sheet
    ws2 = wb.create_sheet("Pairwise Detail")
    detail_headers = ["Tool A", "Tool B", "Both Called", "Concordant", "Discordant", "Rate"]
    for col_idx, hdr in enumerate(detail_headers, start=1):
        c = ws2.cell(row=1, column=col_idx, value=hdr)
        c.font = header_font
        c.fill = header_fill
    for row_idx, stat in enumerate(report.pairwise, start=2):
        ws2.cell(row=row_idx, column=1, value=stat.tool_a)
        ws2.cell(row=row_idx, column=2, value=stat.tool_b)
        ws2.cell(row=row_idx, column=3, value=stat.n_both_called)
        ws2.cell(row=row_idx, column=4, value=stat.n_concordant)
        ws2.cell(row=row_idx, column=5, value=stat.n_discordant)
        ws2.cell(row=row_idx, column=6,
                 value=f"{stat.concordance_rate*100:.1f}%" if stat.concordance_rate is not None else "N/A")

    for col_idx in range(1, 7):
        ws2.column_dimensions[openpyxl.utils.get_column_letter(col_idx)].width = 18

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_path)
    logger.info("Tool concordance matrix written to %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Raw results TSV + JSON
# ---------------------------------------------------------------------------

_RAW_FIELDS = [
    "variant_id", "gene", "hgvs_c", "hgvs_p",
    "expected_classification", "expected_acmg_codes",
    "pipeline_classification", "pipeline_acmg_codes",
    "concordance_status",
    "gnomad_af",
    "cadd_phred", "revel_score", "spliceai_max_score", "alphamissense_score",
    "clinvar_significance", "clinvar_review_status",
    "loftee_lof",
    "classification_confidence",
    "reviewer_flags",
    "excluded", "exclusion_reason",
]


def _result_to_dict(r: VariantResult) -> dict:
    v = r.variant
    te = r.truth_entry
    return {
        "variant_id":               r.variant_id,
        "gene":                     r.gene,
        "hgvs_c":                   (v.hgvs_c or te.hgvs_c) if v else te.hgvs_c,
        "hgvs_p":                   (v.hgvs_p or te.hgvs_p) if v else te.hgvs_p,
        "expected_classification":  te.expected_classification,
        "expected_acmg_codes":      " ".join(te.expected_acmg_codes),
        "pipeline_classification":  r.pipeline_classification,
        "pipeline_acmg_codes":      " ".join(r.pipeline_acmg_codes),
        "concordance_status":       r.concordance_status,
        "gnomad_af":                v.gnomad_af          if v else None,
        "cadd_phred":               v.cadd_phred          if v else None,
        "revel_score":              v.revel_score          if v else None,
        "spliceai_max_score":       v.spliceai_max_score   if v else None,
        "alphamissense_score":      v.alphamissense_score  if v else None,
        "clinvar_significance":     v.clinvar_significance if v else None,
        "clinvar_review_status":    v.clinvar_review_status if v else None,
        "loftee_lof":               v.loftee_lof           if v else None,
        "classification_confidence":v.classification_confidence if v else None,
        "reviewer_flags":           " | ".join(v.reviewer_flags) if (v and v.reviewer_flags) else "",
        "excluded":                 v.excluded             if v else None,
        "exclusion_reason":         v.exclusion_reason     if v else None,
    }


def write_raw_results(
    output_dir: str,
    results: List[VariantResult],
) -> tuple:
    """Write raw_results.tsv and raw_results.json. Returns (tsv_path, json_path)."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = [_result_to_dict(r) for r in results]

    tsv_path = out / "raw_results.tsv"
    with tsv_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(_RAW_FIELDS) + "\n")
        for row in rows:
            fh.write("\t".join(
                "" if row[f] is None else str(row[f]) for f in _RAW_FIELDS
            ) + "\n")

    json_path = out / "raw_results.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2, default=str)

    logger.info(
        "Raw results written: %s (%d rows), %s",
        tsv_path, len(rows), json_path,
    )
    return str(tsv_path), str(json_path)
