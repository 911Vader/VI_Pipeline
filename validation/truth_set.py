"""
truth_set.py — Load and validate the curated variant truth set.

Truth set file format (tab-separated, UTF-8, first row = header):

  Required columns:
    variant_id              chr:pos:ref:alt  (GRCh38 coordinates, exactly as in VCF)
    gene                    Gene symbol (e.g. BRCA1)
    expected_classification Pathogenic | Likely Pathogenic | VUS |
                            Likely Benign | Benign

  Optional columns:
    hgvs_c                  HGVS coding notation   (e.g. c.5266dupC)
    hgvs_p                  HGVS protein notation  (e.g. p.Gln1756ProfsTer74)
    expected_acmg_codes     Space-separated codes  (e.g. "PVS1 PM2 PS3")
    curator_notes           Free text

Lines beginning with '#' are treated as comments and skipped.
"""

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_VALID_CLASSIFICATIONS = {
    "Pathogenic",
    "Likely Pathogenic",
    "VUS",
    "Likely Benign",
    "Benign",
}

_REQUIRED_COLUMNS = {"variant_id", "gene", "expected_classification"}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class TruthEntry:
    """One curated variant from the truth set."""
    variant_id: str                    # chr:pos:ref:alt
    gene: str
    expected_classification: str
    hgvs_c: str = ""
    hgvs_p: str = ""
    expected_acmg_codes: List[str] = field(default_factory=list)
    curator_notes: str = ""

    # Parsed from variant_id
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""

    def __post_init__(self) -> None:
        parts = self.variant_id.split(":")
        if len(parts) == 4:
            self.chrom, pos_s, self.ref, self.alt = parts
            try:
                self.pos = int(pos_s)
            except ValueError:
                logger.warning("Cannot parse pos from variant_id %s", self.variant_id)
        else:
            logger.warning(
                "variant_id %r does not match chr:pos:ref:alt format",
                self.variant_id,
            )

    @property
    def match_key(self) -> Tuple[str, int, str, str]:
        """(chrom, pos, ref, alt) tuple for matching against Variant objects."""
        return (self.chrom, self.pos, self.ref, self.alt)


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

def load_truth_set(path: str) -> List[TruthEntry]:
    """
    Parse a truth set TSV and return a list of TruthEntry objects.

    Raises ValueError if required columns are missing or any row has
    an unrecognised expected_classification value.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Truth set file not found: {path}")

    entries: List[TruthEntry] = []
    errors: List[str] = []

    with p.open(newline="", encoding="utf-8") as fh:
        # Skip comment lines before header
        lines = [ln for ln in fh if not ln.startswith("#")]

    reader = csv.DictReader(lines, delimiter="\t")

    # Validate header
    if reader.fieldnames is None:
        raise ValueError(f"Truth set file {path} appears to be empty")
    missing = _REQUIRED_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(
            f"Truth set file {path} is missing required columns: {sorted(missing)}"
        )

    for row_num, row in enumerate(reader, start=2):
        vid  = row.get("variant_id", "").strip()
        gene = row.get("gene", "").strip()
        cls  = row.get("expected_classification", "").strip()

        if not vid or not gene or not cls:
            errors.append(f"Row {row_num}: required field(s) empty — skipped")
            continue

        if cls not in _VALID_CLASSIFICATIONS:
            errors.append(
                f"Row {row_num}: unrecognised expected_classification {cls!r} — skipped"
            )
            continue

        codes_raw = row.get("expected_acmg_codes", "").strip()
        codes = [c.strip() for c in codes_raw.split() if c.strip()] if codes_raw else []

        entries.append(TruthEntry(
            variant_id             = vid,
            gene                   = gene,
            expected_classification= cls,
            hgvs_c                 = row.get("hgvs_c", "").strip(),
            hgvs_p                 = row.get("hgvs_p", "").strip(),
            expected_acmg_codes    = codes,
            curator_notes          = row.get("curator_notes", "").strip(),
        ))

    if errors:
        for msg in errors:
            logger.warning("truth_set: %s", msg)

    logger.info(
        "Loaded %d truth entries from %s (%d rows skipped with warnings)",
        len(entries), path, len(errors),
    )
    return entries


def build_truth_index(entries: List[TruthEntry]) -> Dict[Tuple, TruthEntry]:
    """
    Return a dict keyed by (chrom, pos, ref, alt) for O(1) lookup.

    Duplicate variant_ids in the truth set produce a warning; the last
    entry wins.
    """
    index: Dict[Tuple, TruthEntry] = {}
    for e in entries:
        key = e.match_key
        if key in index:
            logger.warning(
                "Duplicate truth set entry for %s — later row wins", e.variant_id
            )
        index[key] = e
    return index
