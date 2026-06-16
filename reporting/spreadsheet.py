"""
spreadsheet.py — Export ranked candidates to XLSX using openpyxl.

Column order is designed for non-bioinformatics reviewers:
  - Clinical context first (rank, gene, ACMG, ClinVar)
  - Genomic coordinates next
  - Transcript annotation
  - Population frequency and in-silico scores last

The ACMG Classification column is colour-coded for at-a-glance review:
  Pathogenic       → red
  Likely Pathogenic → orange
  VUS              → yellow
  Likely Benign    → light blue
  Benign           → green
"""

import logging
from pathlib import Path
from typing import List, Optional

from analysis.parser import Variant
from analysis.acmg.evidence_codes import PATHOGENIC, LIKELY_PATHOGENIC, VUS
from reporting import PIPELINE_VERSION

logger = logging.getLogger(__name__)

# Mandatory disclaimer — must not be altered.
MANDATORY_DISCLAIMER = (
    "Top genomic candidate(s) requiring expert review. "
    "Variant pathogenicity alone does not prove cause of death."
)

# Column definitions: (header_label, field_accessor_function)
# Ordered for clinical readability.
# Public name so tsv_export.py can share the exact same column set.
VARIANT_COLUMNS = [
    ("Sample ID",             lambda v, rank, sid, cid: sid),
    ("Case ID",               lambda v, rank, sid, cid: cid),
    ("Gene",                  lambda v, rank, sid, cid: v.gene_symbol or ""),
    ("ACMG Classification",   lambda v, rank, sid, cid: v.acmg_classification or ""),
    ("ACMG Evidence Codes",   lambda v, rank, sid, cid: " ".join(v.acmg_evidence_codes) if v.acmg_evidence_codes else ""),
    ("Classification Confidence", lambda v, rank, sid, cid: v.classification_confidence or ""),
    ("Gene Relevance",        lambda v, rank, sid, cid: v.gene_relevance or ""),
    ("Reviewer Flags",        lambda v, rank, sid, cid: " | ".join(v.reviewer_flags) if v.reviewer_flags else ""),
    ("ClinVar Significance",  lambda v, rank, sid, cid: v.clinvar_significance or ""),
    ("ClinVar Review Status", lambda v, rank, sid, cid: v.clinvar_review_status or ""),
    ("Chromosome",            lambda v, rank, sid, cid: v.chrom),
    ("Position",              lambda v, rank, sid, cid: v.pos),
    ("Ref",                   lambda v, rank, sid, cid: v.ref),
    ("Alt",                   lambda v, rank, sid, cid: v.alt),
    ("Transcript",            lambda v, rank, sid, cid: v.transcript_id or ""),
    ("HGVS c.",               lambda v, rank, sid, cid: v.hgvs_c or ""),
    ("HGVS p.",               lambda v, rank, sid, cid: v.hgvs_p or ""),
    ("Consequence",           lambda v, rank, sid, cid: v.consequence or ""),
    ("Zygosity",              lambda v, rank, sid, cid: v.zygosity or ""),
    ("DP",                    lambda v, rank, sid, cid: v.depth),
    ("GQ",                    lambda v, rank, sid, cid: v.genotype_quality),
    ("AD",                    lambda v, rank, sid, cid: v.allele_depth or ""),
    ("VAF",                   lambda v, rank, sid, cid: v.vaf),
    ("gnomAD AF",             lambda v, rank, sid, cid: v.gnomad_af),
    ("CADD phred",            lambda v, rank, sid, cid: v.cadd_phred),
    ("REVEL",                 lambda v, rank, sid, cid: v.revel_score),
    ("SpliceAI",              lambda v, rank, sid, cid: v.spliceai_max_score),
    ("AlphaMissense",         lambda v, rank, sid, cid: v.alphamissense_score),
    ("SIFT",                  lambda v, rank, sid, cid: v.sift_score),
    ("SIFT Pred",             lambda v, rank, sid, cid: v.sift_pred or ""),
    ("PolyPhen",              lambda v, rank, sid, cid: v.polyphen_score),
    ("PolyPhen Pred",         lambda v, rank, sid, cid: v.polyphen_pred or ""),
    ("LOFTEE",                lambda v, rank, sid, cid: v.loftee_lof or ""),
    ("LOFTEE Filter",         lambda v, rank, sid, cid: v.loftee_filter or ""),
    ("LOFTEE Flags",          lambda v, rank, sid, cid: v.loftee_flags or ""),
    ("ClinGen Validity",      lambda v, rank, sid, cid: v.clingen_validity or ""),
    ("In Gene List",          lambda v, rank, sid, cid: "Yes" if v.in_gene_list else "No"),
    ("Excluded",               lambda v, rank, sid, cid: "Yes" if v.excluded else "No"),
    ("Exclusion Reason",       lambda v, rank, sid, cid: v.exclusion_reason or ""),
    ("Deprioritized",          lambda v, rank, sid, cid: "Yes" if v.deprioritized else "No"),
    ("Deprioritization Reason",lambda v, rank, sid, cid: v.deprioritization_reason or ""),
    ("ClinVar Disease",       lambda v, rank, sid, cid: v.clinvar_disease or ""),
    ("Variant ID",            lambda v, rank, sid, cid: v.variant_id or ""),
]

_HIGH_IMPACT_CSQ = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "start_lost",
    "missense_variant", "splice_region_variant",
    "inframe_insertion", "inframe_deletion",
})

def _is_high_impact(variant: Variant) -> bool:
    if not variant.consequence:
        return False
    parts = {c.strip() for c in variant.consequence.split("&")}
    return bool(parts & _HIGH_IMPACT_CSQ)


_ACMG_COLOURS = {
    "Pathogenic":        "FFCCCC",   # red-ish
    "Likely Pathogenic": "FFE0B2",   # orange-ish
    "VUS":               "FFFDE7",   # yellow-ish
    "Likely Benign":     "E3F2FD",   # light blue
    "Benign":            "E8F5E9",   # light green
}

_HEADER_FILL   = "4472C4"   # blue
_EXCLUDED_FILL = "F5F5F5"   # light grey for excluded rows

# Excel hard limit: 1,048,576 rows per sheet (row 1 = header → 1,048,575 data rows max).
_XL_MAX_DATA_ROWS = 1_048_575

# Separate cap for the Excluded and Deprioritized sheets.
# WGS datasets routinely produce millions of excluded variants (intergenic, high-AF, etc.),
# so writing them all to openpyxl is extremely slow and memory-intensive.
# The full excluded variant list is available in the JSON export.
# Override via settings.yaml: reporting.max_xlsx_excluded_rows
_EXCLUDED_SHEET_DEFAULT_MAX = 100_000


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_spreadsheet(
    variants: List[Variant],
    output_path: str,
    metadata: dict,
    settings: Optional[dict] = None,
) -> str:
    """
    Write ranked variants to an XLSX file.

    Retained variants appear on the first sheet in ranked order.
    Excluded variants appear on a second sheet.
    A 'Info' sheet contains sample metadata and the mandatory disclaimer.

    Parameters
    ----------
    variants:    Ranked list from ranker.rank_variants().
    output_path: Destination .xlsx path.
    metadata:    Dict with sample_id, case_id, etc.
    settings:    Parsed settings.yaml (optional; used for metadata display).

    Returns
    -------
    Path to the written XLSX file (str).
    """
    try:
        import openpyxl
        from openpyxl.styles import PatternFill, Font, Alignment
        from openpyxl.utils import get_column_letter
    except ImportError:
        raise RuntimeError(
            "openpyxl is required for XLSX output.\n"
            "  Install it:  pip install openpyxl"
        )

    sample_id = metadata.get("sample_id", "UNKNOWN")
    case_id   = metadata.get("case_id", "")

    retained = [v for v in variants if not v.excluded and not v.deprioritized]
    excluded = [v for v in variants if v.excluded]

    wb = openpyxl.Workbook()

    # Sheet 1: Info  (written after sheets so truncation_notes is populated)
    ws_info = wb.active
    ws_info.title = "Info"

    # Classify retained variants into report sections
    priority   = [v for v in retained if v.acmg_classification in (PATHOGENIC, LIKELY_PATHOGENIC)]
    vus_list   = [v for v in retained if v.acmg_classification == VUS]
    in_gl      = [v for v in retained if v.in_gene_list and v.acmg_classification not in (PATHOGENIC, LIKELY_PATHOGENIC, VUS)]
    high_out   = [v for v in retained if not v.in_gene_list
                  and v.acmg_classification not in (PATHOGENIC, LIKELY_PATHOGENIC, VUS)
                  and _is_high_impact(v)]
    deprio     = [v for v in variants if not v.excluded and v.deprioritized]

    truncation_notes: list = []   # collected for the Info sheet

    # Excluded/Deprioritized sheets use a lower row cap so XLSX stays fast for WGS datasets.
    # Priority/VUS/gene-list sheets are unlikely to have millions of rows, so they use
    # Excel's hard limit.  Full excluded data is always in the JSON export.
    excluded_cap = int(
        (settings or {}).get("reporting", {}).get("max_xlsx_excluded_rows", _EXCLUDED_SHEET_DEFAULT_MAX)
    )

    for base_title, rows, start_rank in [
        ("Priority Review",  priority,  1),
        ("VUS",              vus_list,  1),
        ("In Gene List",     in_gl,     1),
        ("High Impact",      high_out,  1),
        ("Deprioritized",    deprio,    None),
        ("Excluded",         excluded,  None),
    ]:
        is_large_sheet = base_title in ("Excluded", "Deprioritized")
        if is_large_sheet:
            # 0 means unlimited (up to Excel's hard limit)
            cap = _XL_MAX_DATA_ROWS if excluded_cap == 0 else min(excluded_cap, _XL_MAX_DATA_ROWS)
        else:
            cap = _XL_MAX_DATA_ROWS
        total_rows = len(rows)
        if total_rows > cap:
            logger.info(
                "Sheet '%s': %d rows — showing first %d. Full data available in JSON export.",
                base_title, total_rows, cap,
            )
            rows = rows[:cap]
            title = f"{base_title} (truncated)"
            truncation_notes.append(
                f"{base_title}: showing {cap:,} of {total_rows:,} rows "
                f"— see JSON export for complete data."
            )
        else:
            title = base_title

        ws = wb.create_sheet(title)
        _write_variant_sheet(
            ws, rows, sample_id, case_id, start_rank,
            Font, PatternFill, Alignment, get_column_letter,
            shade_excluded=(base_title in ("Deprioritized", "Excluded")),
        )

    _write_info_sheet(ws_info, metadata, len(retained), len(excluded),
                      Font, PatternFill, Alignment,
                      truncation_notes=truncation_notes)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out))
    logger.info(
        "Spreadsheet written: %s (%d retained, %d excluded)",
        out, len(retained), len(excluded),
    )
    return str(out)


# ---------------------------------------------------------------------------
# Sheet builders
# ---------------------------------------------------------------------------

def _write_info_sheet(ws, metadata, n_retained, n_excluded,
                      Font, PatternFill, Alignment,
                      truncation_notes=None):
    """Write the metadata / disclaimer Info sheet."""
    ws.column_dimensions["A"].width = 22
    ws.column_dimensions["B"].width = 60

    rows = [
        ("VI_Pipeline Report", ""),
        ("Pipeline Version",   PIPELINE_VERSION),
        ("", ""),
        ("MANDATORY DISCLAIMER", MANDATORY_DISCLAIMER),
        ("", ""),
        ("Sample ID",    metadata.get("sample_id", "")),
        ("Case ID",      metadata.get("case_id", "")),
        ("Input VCF",    metadata.get("input_vcf", "")),
        ("Gene List",    metadata.get("gene_list_file", "")),
        ("", ""),
        ("Retained variants", n_retained),
        ("Excluded variants", n_excluded),
    ]
    pheno = metadata.get("phenotype_terms", [])
    if pheno:
        rows.append(("Phenotype terms", "; ".join(str(p) for p in pheno)))

    if truncation_notes:
        rows.append(("", ""))
        rows.append(("NOTE — Excel row limit", ""))
        for note in truncation_notes:
            rows.append(("", note))

    for r, (label, value) in enumerate(rows, start=1):
        cell_a = ws.cell(row=r, column=1, value=label)
        cell_b = ws.cell(row=r, column=2, value=value)
        if r == 1:
            cell_a.font = Font(bold=True, size=14)
        elif label == "MANDATORY DISCLAIMER":
            cell_a.font = Font(bold=True, color="FF0000")
            cell_b.font = Font(bold=True, color="FF0000")
            cell_b.alignment = Alignment(wrap_text=True)
            ws.row_dimensions[r].height = 40
        elif label == "NOTE — Excel row limit":
            cell_a.font = Font(bold=True, color="FF6600")
            cell_b.alignment = Alignment(wrap_text=True)

    ws.sheet_view.showGridLines = False


def _write_variant_sheet(
    ws, variants, sample_id, case_id, start_rank,
    Font, PatternFill, Alignment, get_column_letter,
    shade_excluded,
):
    """Write a variant data sheet with header row and one row per variant."""
    headers = [col[0] for col in VARIANT_COLUMNS]

    # --- Header row ---
    for col_idx, header in enumerate(headers, start=1):
        cell = ws.cell(row=1, column=col_idx, value=header)
        cell.font      = Font(bold=True, color="FFFFFF")
        cell.fill      = PatternFill("solid", fgColor=_HEADER_FILL)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)

    ws.freeze_panes = "A2"

    # --- Data rows ---
    for row_idx, variant in enumerate(variants, start=2):
        rank = (start_rank + row_idx - 2) if start_rank is not None else None

        acmg_class  = variant.acmg_classification or ""
        row_fill    = _EXCLUDED_FILL if shade_excluded else None
        acmg_fill   = _ACMG_COLOURS.get(acmg_class)

        for col_idx, (header, accessor) in enumerate(VARIANT_COLUMNS, start=1):
            value = accessor(variant, rank, sample_id, case_id)
            cell  = ws.cell(row=row_idx, column=col_idx, value=value)
            cell.alignment = Alignment(wrap_text=False)

            if header == "ACMG Classification" and acmg_fill:
                cell.fill = PatternFill("solid", fgColor=acmg_fill)
            elif row_fill:
                cell.fill = PatternFill("solid", fgColor=row_fill)

    # --- Auto-size columns (approximate) ---
    col_widths = {
        "Sample ID": 14, "Case ID": 14,
        "Gene": 12, "ACMG Classification": 18, "ACMG Evidence Codes": 20,
        "Gene Relevance": 22, "Reviewer Flags": 50,
        "ClinVar Significance": 22, "ClinVar Review Status": 28,
        "Chromosome": 10, "Position": 12, "Ref": 6, "Alt": 6,
        "Transcript": 16, "HGVS c.": 28, "HGVS p.": 20,
        "Consequence": 28, "Zygosity": 12,
        "DP": 6, "GQ": 6, "AD": 12, "VAF": 8,
        "gnomAD AF": 12, "CADD phred": 10, "REVEL": 8,
        "SpliceAI": 10, "AlphaMissense": 14,
        "SIFT": 8, "SIFT Pred": 12, "PolyPhen": 10, "PolyPhen Pred": 14,
        "LOFTEE": 8, "LOFTEE Filter": 18, "LOFTEE Flags": 18,
        "ClinGen Validity": 18,
        "In Gene List": 12, "Excluded": 10, "Exclusion Reason": 30,
        "Deprioritized": 14, "Deprioritization Reason": 40,
        "ClinVar Disease": 28, "Variant ID": 16,
    }
    for col_idx, (header, _) in enumerate(VARIANT_COLUMNS, start=1):
        letter = get_column_letter(col_idx)
        ws.column_dimensions[letter].width = col_widths.get(header, 14)
