"""
ranker.py — Order variants by established clinical criteria for VI_Pipeline.

Purpose
-------
Variants are sorted into a deterministic, clinically meaningful order so that
the most significant findings appear at the top of every report.  No composite
score or invented weighting is applied — the sort order is derived entirely
from the data that is actually present for each variant.

Sort order (applied in sequence; ties cascade to the next level):
  1. Retained variants before excluded ones.
  2. ACMG classification:
       Pathogenic > Likely Pathogenic > VUS > Likely Benign > Benign > unclassified
  3. ClinVar significance:
       Pathogenic/Likely Pathogenic first, then all others.
  4. Consequence severity tier:
       tier1 (null variants) > tier2 (missense/splice) > tier3 (UTR/intronic) > unknown
  5. gnomAD allele frequency (ascending — rarer variants rank higher):
       absent from databases (None) → lowest AF → highest AF
  6. Genomic coordinates (chrom → pos → ref → alt) for determinism.

All sort keys are derived directly from the variant's annotation data.
No weights, scores, or thresholds invented by this module are used.

Output consumed by:
  reporting/ — top-N candidate selection and report generation
"""

import logging
from typing import Dict, List, Optional

from analysis.parser import Variant
from analysis.acmg.evidence_codes import (
    PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN,
)

logger = logging.getLogger(__name__)

# ACMG class → sort key (lower = more significant, sorts first).
# Unclassified (None) sorts last.
_ACMG_ORDER: Dict[str, int] = {
    PATHOGENIC:        0,
    LIKELY_PATHOGENIC: 1,
    VUS:               2,
    LIKELY_BENIGN:     3,
    BENIGN:            4,
}
_ACMG_UNCLASSIFIED = 5

# Consequence tier → sort key (lower sorts first).
_TIER_ORDER: Dict[str, int] = {
    "tier1":    0,
    "tier2":    1,
    "tier3":    2,
    "excluded": 3,
}
_TIER_UNKNOWN = 4

# Consequence tier membership — mirrors settings.yaml defaults.
# Used only for ordering; the authoritative tier assignment lives in filter_engine.
_TIER1_CSQ = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "start_lost", "stop_lost",
})
_TIER2_CSQ = frozenset({
    "missense_variant", "splice_region_variant",
    "inframe_insertion", "inframe_deletion",
})
_TIER3_CSQ = frozenset({
    "synonymous_variant", "5_prime_UTR_variant",
    "3_prime_UTR_variant", "intron_variant",
})
_EXCL_CSQ = frozenset({
    "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant",
})


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def rank_variants(variants: List[Variant], settings: dict) -> List[Variant]:
    """
    Sort variants into clinical priority order and return the list.

    Sorting is purely by established clinical criteria — ACMG classification,
    ClinVar significance, consequence tier, gnomAD rarity, and genomic
    coordinates.  No composite score is calculated.

    Parameters
    ----------
    variants: List of Variant objects with acmg_classification already set.
    settings: Parsed config/settings.yaml dict (used to read consequence
              tiers if provided; falls back to module defaults otherwise).

    Returns
    -------
    The same list, sorted in priority order.
    """
    tier_config = settings.get("filtering", {}).get("consequence_tiers", {})

    variants.sort(key=lambda v: _sort_key(v, tier_config))

    n_retained = sum(1 for v in variants if not v.excluded)
    n_excluded = len(variants) - n_retained
    top = next((v for v in variants if not v.excluded), None)
    logger.info(
        "Ordering complete: %d retained / %d excluded — top finding: %s %s",
        n_retained, n_excluded,
        top.gene_symbol if top else "—",
        top.acmg_classification if top else "—",
    )
    return variants


# ---------------------------------------------------------------------------
# Sort key
# ---------------------------------------------------------------------------

def _sort_key(v: Variant, tier_config: dict) -> tuple:
    """
    Build a tuple sort key for one variant.

    All elements are chosen so that Python's default tuple comparison produces
    the correct clinical priority order (most significant variant first).
    """
    # 1. Retained before excluded
    excluded_key = 1 if v.excluded else 0

    # 2. ACMG classification order
    acmg_key = _ACMG_ORDER.get(v.acmg_classification or "", _ACMG_UNCLASSIFIED)

    # 3. ClinVar P/LP first (0), everything else (1)
    clinvar_key = 0 if _is_clinvar_plp(v.clinvar_significance) else 1

    # 4. Consequence tier
    tier_key = _consequence_tier_key(v.consequence, tier_config)

    # 5. gnomAD AF ascending (None = absent from databases = most rare = sorts first)
    #    Use a large sentinel so None sorts before any real AF value.
    af_key = v.gnomad_af if v.gnomad_af is not None else -1.0

    # 6. Genomic coordinates for determinism
    chrom_key = _chrom_sort_key(v.chrom)

    return (excluded_key, acmg_key, clinvar_key, tier_key, af_key, chrom_key, v.pos, v.ref, v.alt)


# ---------------------------------------------------------------------------
# Sort key helpers
# ---------------------------------------------------------------------------

def _is_clinvar_plp(clinvar_significance: Optional[str]) -> bool:
    """Return True if ClinVar significance is Pathogenic or Likely Pathogenic."""
    if not clinvar_significance:
        return False
    sig = clinvar_significance.lower().replace("_", " ")
    if "conflicting" in sig or "uncertain" in sig or "benign" in sig:
        return False
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


def _consequence_tier_key(consequence_str: str, tier_config: dict) -> int:
    """
    Return the tier sort key for a (possibly compound) consequence string.

    Uses tier_config from settings.yaml when available; falls back to
    module-level defaults so ordering works even without settings.
    """
    if not consequence_str:
        return _TIER_UNKNOWN

    parts = {c.strip() for c in consequence_str.split("&") if c.strip()}

    # Use settings-provided tiers if available, otherwise use module defaults
    if tier_config:
        best = _TIER_UNKNOWN
        for tier_name, tier_list in tier_config.items():
            if isinstance(tier_list, list):
                if parts & set(tier_list):
                    key = _TIER_ORDER.get(tier_name, _TIER_UNKNOWN)
                    if key < best:
                        best = key
        return best

    # Module-default fallback
    if parts & _TIER1_CSQ:
        return 0
    if parts & _TIER2_CSQ:
        return 1
    if parts & _TIER3_CSQ:
        return 2
    if parts & _EXCL_CSQ:
        return 3
    return _TIER_UNKNOWN


def _chrom_sort_key(chrom: str) -> tuple:
    """
    Return a sort key for a chromosome string that orders numerically then
    alphabetically: 1, 2, …, 22, X, Y, MT, then anything else.
    """
    c = chrom.lstrip("chr").upper()
    try:
        return (0, int(c), "")
    except ValueError:
        order = {"X": 23, "Y": 24, "MT": 25, "M": 25}
        return (1, order.get(c, 99), c)


# ---------------------------------------------------------------------------
# Reporting helpers (used by run_log and JSON export)
# ---------------------------------------------------------------------------

def get_ranking_summary(variants: List[Variant]) -> dict:
    """
    Return a count summary of ranked variants by ACMG class.

    Useful for run-log reporting.
    """
    retained = [v for v in variants if not v.excluded]
    excluded = [v for v in variants if v.excluded]
    by_class: Dict[str, int] = {}
    for v in retained:
        cat = v.acmg_classification or "unclassified"
        by_class[cat] = by_class.get(cat, 0) + 1
    top_n = [
        {"gene": v.gene_symbol, "acmg": v.acmg_classification}
        for v in retained[:5]
    ]
    return {
        "total":         len(variants),
        "retained":      len(retained),
        "excluded":      len(excluded),
        "by_acmg_class": by_class,
        "top_5":         top_n,
    }
