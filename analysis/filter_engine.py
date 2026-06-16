"""
filter_engine.py — Apply configurable filtering rules to annotated Variant objects.

Design decisions (v1):
  - AF source: variant.gnomad_af is used directly.  This field is set by
    parser.py as MAX_AF (primary — most conservative; max across all gnomAD
    populations) with fallback to gnomAD_AF when MAX_AF is absent.
    No secondary AF lookup is performed here.
  - Transcript scope: filtering uses variant.consequence (primary transcript
    only).  This is explicit v1 behavior; per-transcript consequence logic
    can be added in a future version when needed.
  - Gene-list matching: case-insensitive exact match on variant.gene_symbol.
    Partial and prefix matching are not performed.
  - ClinVar exclusion: variants with a configured benign/likely-benign
    ClinVar significance are excluded.  The reason is recorded in
    variant.exclusion_reason.
  - PASS filter: when require_pass=True, only filter_status=="PASS" is
    retained.  filter_status=="." is treated as excluded in v1 — this is
    the stricter, unambiguous behavior and is noted here explicitly so it
    can be relaxed later.
  - Consequence-tier logic (from settings.yaml):
      tier1 ("Always retain")        → never excluded by consequence filter
      tier2 ("Retain unless benign") → not excluded by consequence filter
                                       (ClinVar benign filter applies separately)
      tier3 ("Gene list / ClinVar")  → excluded UNLESS variant.in_gene_list
                                       OR ClinVar significance is P / LP
      excluded ("Drop")              → always excluded
  - Variants are NEVER removed from the list.  variant.excluded and
    variant.exclusion_reason are set for downstream consumers.

Rules are read from config/settings.yaml.  No thresholds are hardcoded here.
"""

import logging
from typing import Dict, List, Optional, Set

from analysis.parser import Variant

logger = logging.getLogger(__name__)

# Consequence tier priority (lower integer = more favorable / retained first).
# Used to select the best tier when a compound consequence string matches
# multiple tiers (e.g. "missense_variant&intron_variant" → tier2 wins).
_TIER_PRIORITY: Dict[str, int] = {
    "tier1":    1,
    "tier2":    2,
    "tier3":    3,
    "excluded": 4,
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def apply_filters(
    variants: List[Variant],
    settings: dict,
    gene_list: Optional[List[str]] = None,
) -> List[Variant]:
    """
    Apply all configured filtering rules to a list of Variant objects.

    Variants are never removed from the list.  Each variant's .excluded and
    .exclusion_reason fields are set according to the filter outcome.

    Filtering order (first failing check determines the exclusion reason):
      1. PASS-only filter (if require_pass=True in settings)
      2. Allele-frequency filter   (gnomad_af > max_af → exclude)
      3. Consequence-tier filter   (excluded tier or tier3 without justification)
      4. ClinVar benign filter     (significance in exclude_clinvar_significances)

    Gene-list membership (variant.in_gene_list) is set for ALL variants before
    any filter is applied, so downstream code can always read it regardless of
    exclusion status.

    Parameters
    ----------
    variants:
        List of Variant objects from parser.parse_vep_output().
    settings:
        Parsed config/settings.yaml (dict).
    gene_list:
        Optional list of gene symbols for gene-list matching.  Matching is
        case-insensitive and exact.  If None or empty, variant.in_gene_list
        remains False for all variants.

    Returns
    -------
    The same list (mutations applied in place).  Retained variants have
    excluded=False; excluded variants have excluded=True and exclusion_reason
    set to a human-readable string.
    """
    filter_cfg       = settings.get("filtering", {})
    require_pass     = bool(filter_cfg.get("require_pass", True))
    max_af           = float(filter_cfg.get("max_af", 0.01))
    ba1_threshold    = float(
        settings.get("acmg", {}).get("ba1_af_threshold", 0.05)
    )
    tier_config      = filter_cfg.get("consequence_tiers", {})
    exclude_sigs_raw: List[str] = filter_cfg.get("exclude_clinvar_significances", [])
    exclude_sigs: Set[str] = {s.lower().replace("_", " ") for s in exclude_sigs_raw}

    gene_set: Set[str] = set()
    if gene_list:
        gene_set = {g.strip().upper() for g in gene_list if g.strip()}

    n_retained = n_excluded = n_deprioritized = 0

    for variant in variants:
        _mark_gene_list(variant, gene_set)
        _set_gene_relevance(variant, gene_set)

        # Hard exclusions — only these produce excluded=True
        excl_reason = _check_pass_filter(variant, require_pass) or \
                      _check_consequence_excluded(variant, tier_config) or \
                      _check_af_exclusion(variant, max_af)

        if excl_reason:
            variant.excluded = True
            variant.exclusion_reason = excl_reason
            n_excluded += 1
            variant.decision_trace["filter"] = {
                "outcome": "excluded",
                "reason": excl_reason,
                "consequence_tier": _get_consequence_tier(variant.consequence, tier_config),
                "gnomad_af_used": variant.gnomad_af,
                "in_gene_list": variant.in_gene_list,
            }
            continue

        # Soft deprioritizations — retained but ranked lower / shown separately
        deprio_reason = _check_tier3_deprioritize(variant, tier_config) or \
                        _check_clinvar_deprioritize(variant, exclude_sigs)

        if deprio_reason:
            variant.deprioritized = True
            variant.deprioritization_reason = deprio_reason
            n_deprioritized += 1
            outcome = "deprioritized"
        else:
            n_retained += 1
            outcome = "retained"

        variant.decision_trace["filter"] = {
            "outcome": outcome,
            "reason": deprio_reason,
            "consequence_tier": _get_consequence_tier(variant.consequence, tier_config),
            "gnomad_af_used": variant.gnomad_af,
            "in_gene_list": variant.in_gene_list,
        }

    logger.info(
        "Filtering complete: %d retained, %d deprioritized, %d excluded (total %d)",
        n_retained, n_deprioritized, n_excluded, len(variants),
    )
    return variants


# ---------------------------------------------------------------------------
# Gene-list membership
# ---------------------------------------------------------------------------

def _mark_gene_list(variant: Variant, gene_set: Set[str]) -> None:
    """
    Set variant.in_gene_list=True if the variant's gene symbol is in gene_set.

    Matching is case-insensitive and exact (no partial/prefix matching).
    gene_set must be pre-normalised to uppercase before calling this function.
    Variants with an empty gene_symbol are never marked as in the gene list.
    """
    if gene_set and variant.gene_symbol:
        variant.in_gene_list = variant.gene_symbol.strip().upper() in gene_set


# ---------------------------------------------------------------------------
# Filter predicates
# Each returns a non-empty exclusion-reason string if the variant is excluded,
# or None if the variant passes this check.
# ---------------------------------------------------------------------------

def _check_pass_filter(variant: Variant, require_pass: bool) -> Optional[str]:
    """
    Exclude variants that do not have FILTER=PASS when require_pass=True.

    v1 behavior: filter_status=="." is treated as excluded (not PASS).
    This is the stricter interpretation and is documented here explicitly
    so it can be changed to "retain" in a future version if needed.
    """
    if not require_pass:
        return None
    if variant.filter_status == "PASS":
        return None
    if variant.filter_status == ".":
        return (
            "FILTER=. (missing filter status); require_pass=True excludes "
            "variants without an explicit PASS in v1"
        )
    return (
        f"FILTER={variant.filter_status!r}; require_pass=True in settings — "
        f"only PASS variants are retained"
    )



def _set_gene_relevance(variant: Variant, gene_set: Set[str]) -> None:
    """Set variant.gene_relevance based on gene list membership."""
    if not gene_set:
        variant.gene_relevance = ""
        return
    variant.gene_relevance = (
        "In supplied gene list"
        if variant.in_gene_list
        else "Outside supplied gene list"
    )


def _check_consequence_excluded(
    variant: Variant,
    tier_config: dict,
) -> Optional[str]:
    """Return exclusion reason if consequence is in the 'excluded' tier."""
    tier = _get_consequence_tier(variant.consequence, tier_config)
    if tier == "excluded":
        return (
            f"Consequence {variant.consequence!r} is in the excluded tier "
            f"(settings: filtering.consequence_tiers.excluded)"
        )
    return None


def _check_af_exclusion(variant: Variant, max_af: float) -> Optional[str]:
    """Exclude variants with gnomAD AF above the SOP common-variant cutoff (default 1%)."""
    if variant.gnomad_af is None:
        return None
    if variant.gnomad_af > max_af:
        return (
            f"gnomAD AF {variant.gnomad_af:.6g} > {max_af} "
            f"(common variant — excluded per SOP)"
        )
    return None


def _check_tier3_deprioritize(
    variant: Variant,
    tier_config: dict,
) -> Optional[str]:
    """Deprioritize tier3 consequences (synonymous, UTR, intronic) rather than exclude them."""
    tier = _get_consequence_tier(variant.consequence, tier_config)
    if tier == "tier3":
        return (
            f"Consequence {variant.consequence!r} is tier3 "
            f"(synonymous/UTR/intronic — deprioritized; retained for review)"
        )
    return None


def _check_clinvar_deprioritize(
    variant: Variant,
    exclude_sigs: Set[str],
) -> Optional[str]:
    """Deprioritize (not exclude) variants with ClinVar benign/likely-benign significance."""
    if not variant.clinvar_significance:
        return None
    # ClinVar CLNSIG can be slash-separated (e.g. "Benign/Likely_benign").
    # Normalize each component and use exact set membership — substring matching
    # would false-positive on values like "not_provided" containing a listed term.
    components = {
        part.strip().lower().replace("_", " ")
        for part in variant.clinvar_significance.replace(",", "/").split("/")
        if part.strip()
    }
    if components & exclude_sigs:
        return (
            f"ClinVar significance {variant.clinvar_significance!r} indicates benign; "
            f"deprioritized but retained for traceability"
        )
    return None


# Alias preserved for backward compatibility with any existing test imports.
_check_clinvar_filter = _check_clinvar_deprioritize


# ---------------------------------------------------------------------------
# Consequence-tier helpers
# ---------------------------------------------------------------------------

def _get_consequence_tier(
    consequence_str: str,
    tier_config: dict,
) -> Optional[str]:
    """
    Return the most favorable tier for a (possibly compound) consequence string.

    consequence_str may contain multiple consequences joined by '&'
    (e.g. 'missense_variant&splice_region_variant').  Each individual consequence
    is looked up independently; the most favorable (lowest _TIER_PRIORITY) tier
    is returned.

    Example: 'missense_variant&intron_variant' → tier2 (missense wins over tier3)

    Returns None if no individual consequence maps to any configured tier.
    """
    if not consequence_str:
        return None

    individuals = [c.strip() for c in consequence_str.split("&") if c.strip()]

    best_tier: Optional[str] = None
    best_priority = 9999

    for csq in individuals:
        for tier_name, tier_consequences in tier_config.items():
            if not isinstance(tier_consequences, list):
                continue
            if csq in tier_consequences:
                priority = _TIER_PRIORITY.get(tier_name, 9999)
                if priority < best_priority:
                    best_priority = priority
                    best_tier = tier_name

    return best_tier


def _is_clinvar_pathogenic(clinvar_significance: Optional[str]) -> bool:
    """
    Return True if ClinVar significance indicates Pathogenic or Likely_pathogenic.

    Handles common ClinVar formats including underscore-separated values and
    mixed case.

    Explicitly NOT treated as pathogenic:
      - "Conflicting_interpretations_of_pathogenicity" (contains "pathogenic" as
        a substring but is not a pathogenic classification)
      - "Uncertain_significance"
      - Any "Benign" or "Likely_benign" value
    """
    if not clinvar_significance:
        return False
    sig = clinvar_significance.lower().replace("_", " ")
    # Reject ambiguous / conflicting classifications even if they mention pathogenic
    if "conflicting" in sig or "uncertain" in sig or "benign" in sig:
        return False
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


# ---------------------------------------------------------------------------
# Reporting helper
# ---------------------------------------------------------------------------

def get_filter_summary(variants: List[Variant], top_excluded_n: int = 20) -> dict:
    """
    Return a summary dict of filtering outcomes.

    Useful for run-log reporting.  Exclusion reasons are bucketed by their
    leading keyword for compact display.  The top_excluded list makes the
    first N excluded variants visible without requiring log inspection.

    Returns
    -------
    dict with keys: total, retained, deprioritized, excluded,
    exclusion_reasons (bucketed-reason → count),
    top_excluded (list of up to top_excluded_n excluded variant summaries).
    """
    excluded_list      = [v for v in variants if v.excluded]
    deprioritized_list = [v for v in variants if not v.excluded and v.deprioritized]
    retained_count     = len(variants) - len(excluded_list) - len(deprioritized_list)

    reason_counts: Dict[str, int] = {}
    for v in excluded_list:
        raw    = v.exclusion_reason or "unknown"
        bucket = raw.split(";")[0].split(" ")[0].rstrip("=").strip()
        reason_counts[bucket] = reason_counts.get(bucket, 0) + 1

    top_excluded = [
        {
            "variant":     f"{v.chrom}:{v.pos} {v.ref}>{v.alt}",
            "gene":        v.gene_symbol or "",
            "consequence": v.consequence or "",
            "gnomad_af":   v.gnomad_af,
            "reason":      v.exclusion_reason or "unknown",
        }
        for v in excluded_list[:top_excluded_n]
    ]

    return {
        "total":             len(variants),
        "retained":          retained_count,
        "deprioritized":     len(deprioritized_list),
        "excluded":          len(excluded_list),
        "exclusion_reasons": reason_counts,
        "top_excluded":      top_excluded,
    }
