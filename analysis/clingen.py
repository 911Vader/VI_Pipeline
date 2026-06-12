"""
clingen.py — Load ClinGen Gene Disease Validity Curations and annotate variants.

Source: https://search.clinicalgenome.org/kb/gene-validity/download
File:   ClinGen_gene_validity.csv  (refreshed each run from paths.yaml)

The CSV has 4 metadata header lines, then the column header row, then data.
Each row is one gene-disease pair. A gene can appear multiple times (different
diseases). This module retains only the strongest classification per gene.

Validity classification hierarchy (strongest → weakest):
  Definitive > Strong > Moderate > Limited > Disputed > Refuted >
  No Known Disease Relationship

Only runs when clingen_gene_validity_csv is configured in paths.yaml.
Output: sets Variant.clingen_validity on each variant in-place.
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional

from analysis.parser import Variant

logger = logging.getLogger(__name__)

# Ordered strongest → weakest. Lower rank = stronger evidence.
_CLASSIFICATION_RANK: Dict[str, int] = {
    "Definitive":                    0,
    "Strong":                        1,
    "Moderate":                      2,
    "Limited":                       3,
    "Disputed":                      4,
    "Refuted":                       5,
    "No Known Disease Relationship": 6,
}

# Column indices in the data rows (0-based, after stripping CSV quoting)
_COL_GENE   = 0
_COL_CLASS  = 6


def load_clingen(csv_path: str) -> Dict[str, str]:
    """
    Parse the ClinGen gene validity CSV into {gene_symbol: best_classification}.

    For genes with multiple disease curations, only the strongest classification
    is kept. Returns an empty dict if the file cannot be read.
    """
    path = Path(csv_path)
    if not path.exists():
        logger.warning("ClinGen file not found: %s — ClinGen lookup skipped", csv_path)
        return {}

    gene_best: Dict[str, str] = {}
    n_rows = 0

    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.reader(fh)
            header_found = False
            for row in reader:
                if not row:
                    continue
                # Skip metadata lines until the actual column header row
                if not header_found:
                    if row[0].strip().strip('"') == "GENE SYMBOL":
                        header_found = True
                    continue

                if len(row) <= _COL_CLASS:
                    continue

                gene           = row[_COL_GENE].strip().strip('"')
                classification = row[_COL_CLASS].strip().strip('"')

                if not gene or classification not in _CLASSIFICATION_RANK:
                    continue

                n_rows += 1
                existing = gene_best.get(gene)
                if existing is None or (
                    _CLASSIFICATION_RANK[classification] < _CLASSIFICATION_RANK[existing]
                ):
                    gene_best[gene] = classification

    except Exception as exc:
        logger.error("Failed to load ClinGen file %s: %s", csv_path, exc)
        return {}

    logger.info(
        "ClinGen loaded: %d curations → %d unique genes (strongest classification kept per gene)",
        n_rows, len(gene_best),
    )
    return gene_best


def annotate_clingen(variants: List[Variant], gene_dict: Dict[str, str]) -> int:
    """
    Set Variant.clingen_validity by looking up each variant's gene symbol.

    Modifies variants in-place. Returns the number of variants that matched.
    Variants with no gene symbol or no ClinGen entry keep clingen_validity=None.
    """
    if not gene_dict:
        return 0

    n_matched = 0
    for v in variants:
        cls = gene_dict.get(v.gene_symbol or "")
        if cls:
            v.clingen_validity = cls
            n_matched += 1

    logger.info(
        "ClinGen annotation: %d / %d variants matched a curated gene",
        n_matched, len(variants),
    )
    return n_matched
