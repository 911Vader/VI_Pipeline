"""
tsv_export.py — Export retained variants to a tab-delimited text file.

One row per retained (non-excluded) variant, in ranked order.
Designed for direct loading into R, pandas, Excel, or any tabular tool.

Format:
  - UTF-8, tab-delimited, Unix line endings.
  - Line 1: mandatory disclaimer as a # comment.
  - Line 2: column headers (tab-separated).
  - Lines 3+: one variant per row, retained variants only, in rank order.

Column set is identical to the XLSX Candidates sheet (VARIANT_COLUMNS).

Output consumed by: downstream analysts, R/Python pipelines, data review.
"""

import csv
import logging
from pathlib import Path
from typing import List, Optional

from analysis.parser import Variant
from reporting import PIPELINE_VERSION
from reporting.spreadsheet import VARIANT_COLUMNS

logger = logging.getLogger(__name__)

#: Mandatory disclaimer — must appear verbatim. Do not alter.
MANDATORY_DISCLAIMER = (
    "Top genomic candidate(s) requiring expert review. "
    "Variant pathogenicity alone does not prove cause of death."
)


def write_tsv(
    variants: List[Variant],
    output_path: str,
    metadata: dict,
    settings: Optional[dict] = None,
) -> str:
    """
    Write retained variants to a TSV file.

    Parameters
    ----------
    variants:    Ranked list from ranker.rank_variants().  Excluded variants
                 are silently omitted from the output.
    output_path: Destination .tsv path.  Parent directory is created if needed.
    metadata:    Dict with sample_id, case_id, etc.
    settings:    Parsed settings.yaml (accepted for API compatibility; unused).

    Returns
    -------
    Path to the written TSV file (str).
    """
    sample_id = metadata.get("sample_id", "UNKNOWN")
    case_id   = metadata.get("case_id", "")

    retained = [v for v in variants if not v.excluded]

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# VI_Pipeline {PIPELINE_VERSION}\n")
        # Mandatory disclaimer as a comment line
        fh.write(f"# {MANDATORY_DISCLAIMER}\n")

        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL,
                            lineterminator="\n")

        # Header row
        writer.writerow([col[0] for col in VARIANT_COLUMNS])

        # One data row per retained variant, in ranked order
        for row_idx, variant in enumerate(retained, start=1):
            rank = row_idx
            row = [
                accessor(variant, rank, sample_id, case_id)
                for _, accessor in VARIANT_COLUMNS
            ]
            writer.writerow(row)

    logger.info(
        "TSV written: %s  (%d retained variants)", out, len(retained)
    )
    return str(out)
