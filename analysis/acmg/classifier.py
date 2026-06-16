"""
classifier.py — Apply ACMG/AMP rules to annotated Variant objects.

Orchestrates:
  1. Evaluating each implemented ACMG/AMP criterion via rules.py
  2. Collecting triggered evidence codes
  3. Combining codes into a final classification via evidence_codes.classify()
  4. Populating variant.acmg_evidence_codes and variant.acmg_classification

Design decisions (v1):
  - Only the rule functions listed in rules.IMPLEMENTED_RULES are evaluated.
    Unimplementable criteria (rules.UNIMPLEMENTED_V1_RULES) are never triggered.
  - A rule that raises an exception logs a warning and is skipped for that
    variant — it does not propagate the error.  Missing data returns no code.
  - The classify_all() function operates on all variants, including those marked
    excluded=True by filter_engine.  ACMG classification is applied regardless
    of exclusion status so that the information is available for reporting.
  - variant.acmg_evidence_codes is always a list (empty if no codes triggered).
  - variant.acmg_classification is always one of the five ACMG categories.

Output consumed by:
  analysis/ranker.py — uses acmg_classification as a ranking signal
  reporting/          — includes codes and classification in all output formats
"""

import logging
from typing import List, Optional, Tuple

from analysis.parser import Variant
from analysis.acmg.evidence_codes import (
    EVIDENCE_CATALOGUE,
    ClassificationCategory,
    PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN,
    classify,
)
from analysis.acmg.rules import IMPLEMENTED_RULES, UNIMPLEMENTED_V1_RULES

logger = logging.getLogger(__name__)

# ClinGen validity categories that lower confidence (gene-disease evidence insufficient)
_CLINGEN_LOW_CONFIDENCE = frozenset({"limited", "disputed", "refuted"})


def _compute_confidence(
    category: ClassificationCategory,
    variant: Variant,
) -> str:
    """
    Derive a classification confidence tier (High / Moderate / Low).

    Base tier from ACMG category:
      Pathogenic / Benign        → High
      Likely Pathogenic / Benign → Moderate
      VUS                        → Low  (inherently uncertain)

    Override to Low when ClinGen gene-disease validity is Limited, Disputed,
    or Refuted — the gene-disease relationship is insufficiently established
    for the ACMG codes to be interpreted reliably.
    """
    if category in (PATHOGENIC, BENIGN):
        base = "High"
    elif category in (LIKELY_PATHOGENIC, LIKELY_BENIGN):
        base = "Moderate"
    else:
        base = "Low"

    if variant.clingen_validity:
        validity_norm = variant.clingen_validity.lower().strip()
        if validity_norm in _CLINGEN_LOW_CONFIDENCE:
            return "Low"

    return base


def classify_variant(
    variant: Variant,
    settings: dict,
) -> Tuple[List[str], ClassificationCategory, str]:
    """
    Evaluate all active ACMG/AMP criteria for a single Variant.

    Iterates over IMPLEMENTED_RULES, evaluates each rule function, collects
    triggered evidence codes, and combines them via evidence_codes.classify().

    Parameters
    ----------
    variant:  An annotated Variant object (from parser.py, optionally filtered).
    settings: Parsed config/settings.yaml dict (supplies ACMG thresholds).

    Returns
    -------
    (evidence_codes, classification_category, confidence)
      evidence_codes:          List of triggered code strings (may be empty).
      classification_category: One of "Pathogenic", "Likely Pathogenic", "VUS",
                               "Likely Benign", "Benign".
      confidence:              "High", "Moderate", or "Low".
    """
    triggered: List[str] = []
    skipped: List[str] = []

    for code, rule_fn in IMPLEMENTED_RULES.items():
        try:
            if rule_fn(variant, settings):
                triggered.append(code)
                logger.debug(
                    "%s:%s — %s triggered",
                    variant.chrom, variant.pos, code,
                )
        except Exception as exc:
            # Rule evaluation errors must never crash the pipeline.
            # Missing or malformed fields produce no code, not an exception.
            # Tracked separately so the output distinguishes "did not fire"
            # from "threw an exception and was skipped".
            logger.warning(
                "Rule %s raised an unexpected error for %s:%s (%s: %s) — skipping",
                code, variant.chrom, variant.pos, type(exc).__name__, exc,
            )
            skipped.append(code)

    # Propagate skipped rules so they are visible in per-variant output.
    variant.acmg_skipped_rules = skipped
    if skipped:
        logger.warning(
            "%s:%s — %d ACMG rule(s) skipped due to errors: %s",
            variant.chrom, variant.pos, len(skipped), skipped,
        )

    category = classify(triggered)
    confidence = _compute_confidence(category, variant)

    # Populate decision trace for this variant.
    variant.decision_trace["acmg"] = {
        "rules_evaluated": list(IMPLEMENTED_RULES.keys()),
        "rules_triggered": triggered,
        "rules_skipped_due_to_error": skipped,
        "classification": str(category),
        "confidence": confidence,
        "scores_used": {
            "gnomad_af":          variant.gnomad_af,
            "cadd_phred":         variant.cadd_phred,
            "revel_score":        variant.revel_score,
            "spliceai_max_score": variant.spliceai_max_score,
            "alphamissense_score":variant.alphamissense_score,
        },
    }

    logger.debug(
        "%s:%s — codes=%s classification=%s confidence=%s skipped=%s",
        variant.chrom, variant.pos, triggered, category, confidence, skipped,
    )
    return triggered, category, confidence


def classify_all(
    variants: List[Variant],
    settings: dict,
) -> List[Variant]:
    """
    Run ACMG/AMP classification over a list of Variant objects.

    Populates variant.acmg_evidence_codes and variant.acmg_classification
    in place for every variant, including those marked excluded=True.

    Parameters
    ----------
    variants: List of Variant objects (from parser and/or filter_engine).
    settings: Parsed config/settings.yaml dict.

    Returns
    -------
    The same list with acmg_evidence_codes and acmg_classification set on
    each element.
    """
    counts = {PATHOGENIC: 0, LIKELY_PATHOGENIC: 0, VUS: 0, LIKELY_BENIGN: 0, BENIGN: 0}

    for variant in variants:
        codes, category, confidence = classify_variant(variant, settings)
        variant.acmg_evidence_codes       = codes
        variant.acmg_classification       = category
        variant.classification_confidence = confidence
        counts[category] = counts.get(category, 0) + 1

    logger.info(
        "ACMG classification complete: %d variants — P=%d LP=%d VUS=%d LB=%d B=%d",
        len(variants),
        counts.get(PATHOGENIC, 0),
        counts.get(LIKELY_PATHOGENIC, 0),
        counts.get(VUS, 0),
        counts.get(LIKELY_BENIGN, 0),
        counts.get(BENIGN, 0),
    )
    return variants


# ---------------------------------------------------------------------------
# Reporting helper
# ---------------------------------------------------------------------------

def get_classification_summary(variants: List[Variant]) -> dict:
    """
    Return a count summary of ACMG classifications across a variant list.

    Useful for run-log reporting.  Variants without acmg_classification set
    (None) are counted separately.

    Returns
    -------
    dict with keys: total, by_category (dict of category → count),
    unclassified (count of variants with acmg_classification=None).
    """
    by_category: dict = {}
    unclassified = 0

    for v in variants:
        cat = v.acmg_classification
        if cat is None:
            unclassified += 1
        else:
            by_category[cat] = by_category.get(cat, 0) + 1

    # Aggregate skipped-rule counts so the run log can flag abnormal runs.
    skipped_counts: dict = {}
    variants_with_skipped = 0
    for v in variants:
        if v.acmg_skipped_rules:
            variants_with_skipped += 1
            for code in v.acmg_skipped_rules:
                skipped_counts[code] = skipped_counts.get(code, 0) + 1

    return {
        "total":                    len(variants),
        "by_category":              by_category,
        "unclassified":             unclassified,
        "variants_with_rule_errors":variants_with_skipped,
        "rule_error_counts":        skipped_counts,
    }
