"""
rules.py — Individual ACMG/AMP criterion evaluation functions (v1 subset).

Each public function evaluates one ACMG/AMP criterion for a single Variant
and returns True (criterion fires) or False (criterion does not fire).
All functions are independently testable.

v1 design principles:
  - Only trigger a code when the underlying variant data supports it.
  - Missing fields (None) must never cause a code to fire or an exception.
  - Every implemented rule has a comment describing exactly what data it reads
    and what limitation applies.

v1 implemented rules
────────────────────
  PVS1  Pathogenic Very Strong  — consequence-based null-variant detection
  PS3   Pathogenic Strong       — ClinVar P/LP ≥2-star review-status proxy
  PM2   Pathogenic Moderate     — absent / extremely rare in gnomAD
  PM4   Pathogenic Moderate     — protein-length change (in-frame indel / stop-loss)
  PP3   Pathogenic Supporting   — ≥2 in-silico tools predict damaging
  BA1   Benign Stand-Alone      — gnomAD AF > 5 %
  BS1   Benign Strong           — gnomAD AF > disorder-specific threshold
  BP4   Benign Supporting       — ≥2 in-silico tools predict benign
  BP7   Benign Supporting       — synonymous variant with no predicted splice impact

v1 NOT implemented — data not available from single-sample VCF + VEP output
────────────────────────────────────────────────────────────────────────────
  PS1   Requires per-amino-acid-position ClinVar database lookup (not available)
  PM1   Requires curated mutational-hotspot / functional-domain database
  PM5   Requires per-codon ClinVar lookup (novel-missense at known-pathogenic codon)
  PP2   Requires per-gene missense-constraint data (gnomAD pLoF/missense metrics)
  BS2   Requires healthy-adult zygosity data (gnomAD homozygous counts by phenotype)

Thresholds used by rules are read from settings.yaml under acmg.insilico_thresholds
and acmg.*_af_threshold keys.  Defaults are shown in each function docstring.

Output consumed by:
  analysis/acmg/classifier.py — collects triggered codes per variant
"""

import logging
from typing import List, Optional

from analysis.parser import Variant

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Consequence sets
# ---------------------------------------------------------------------------

# PVS1: canonical null-variant consequence types.
# Limitation (v1): we do NOT check whether LoF is the established disease
# mechanism for the specific gene — that requires a curated gene LoF database
# (e.g., ClinGen LoF curation).  PVS1 will fire on null-variant types
# regardless of gene.  Downstream users must apply clinical judgement.
_PVS1_CONSEQUENCES = frozenset({
    "stop_gained",           # nonsense
    "frameshift_variant",    # frameshift insertion/deletion
    "splice_acceptor_variant",  # canonical ±1/2 splice site
    "splice_donor_variant",
    "start_lost",            # initiation codon loss
})

# PM4: protein-length-changing variant consequence types.
# stop_lost is listed separately from PVS1 because stop-loss creates a
# longer protein (read-through) rather than a truncated one.
_PM4_CONSEQUENCES = frozenset({
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",             # read-through / protein elongation
})

# BP7: synonymous-variant consequence strings.
_BP7_CONSEQUENCES = frozenset({
    "synonymous_variant",
})


# ---------------------------------------------------------------------------
# Helper: ClinVar significance and review status
# ---------------------------------------------------------------------------

def _is_clinvar_p_or_lp(clinvar_significance: Optional[str]) -> bool:
    """
    Return True if ClinVar significance is Pathogenic or Likely_pathogenic.

    Handles underscore-separated ClinVar formats and mixed case.
    Explicitly excludes "Conflicting_interpretations_of_pathogenicity" —
    this string contains "pathogenic" but is not a pathogenic classification.
    """
    if not clinvar_significance:
        return False
    sig = clinvar_significance.lower().replace("_", " ")
    # Reject ambiguous classifications even if they contain "pathogenic"
    if "conflicting" in sig or "uncertain" in sig or "benign" in sig:
        return False
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


_CLINVAR_TWO_STAR_MARKERS = (
    "multiple_submitters",      # criteria_provided,_multiple_submitters,_no_conflicts
    "expert_panel",             # reviewed_by_expert_panel
    "practice_guideline",       # practice_guideline
)

def _clinvar_review_stars(review_status: Optional[str]) -> int:
    """
    Map a ClinVar review-status string to an approximate star rating (0–4).

    Source: ClinVar review-status documentation
      0 stars — no assertion criteria / no classification
      1 star  — criteria provided, single submitter
      2 stars — criteria provided, multiple submitters, no conflicts
      3 stars — reviewed by expert panel
      4 stars — practice guideline
    """
    if not review_status:
        return 0
    rs = review_status.lower().replace(" ", "_")
    if "practice_guideline" in rs:
        return 4
    if "expert_panel" in rs:
        return 3
    if "multiple_submitters" in rs and "no_conflicts" in rs:
        return 2
    if "criteria_provided" in rs and "single_submitter" in rs:
        return 1
    return 0


# ---------------------------------------------------------------------------
# In-silico scoring helpers (shared by PP3 and BP4)
# ---------------------------------------------------------------------------

def _insilico_thresholds(settings: dict) -> dict:
    """Extract in-silico thresholds from settings with safe defaults."""
    t = settings.get("acmg", {}).get("insilico_thresholds", {})
    return {
        "cadd_pathogenic":   float(t.get("cadd_phred_pathogenic",   20)),
        "revel_pathogenic":  float(t.get("revel_pathogenic",        0.75)),
        "revel_benign":      float(t.get("revel_benign",            0.15)),
        "spliceai_high":     float(t.get("spliceai_high",           0.5)),
        "spliceai_low":      float(t.get("spliceai_low",            0.1)),
        "am_pathogenic":     float(t.get("alphamissense_pathogenic", 0.564)),
        "am_benign":         float(t.get("alphamissense_benign",     0.340)),
    }


def _af_thresholds(settings: dict) -> dict:
    """Extract AF thresholds for ACMG BA1, BS1, PM2 from settings with defaults."""
    acmg_cfg = settings.get("acmg", {})
    return {
        "ba1": float(acmg_cfg.get("ba1_af_threshold", 0.05)),    # 5% — stand-alone benign
        "bs1": float(acmg_cfg.get("bs1_af_threshold", 0.01)),    # 1% — above expected
        "pm2": float(acmg_cfg.get("pm2_af_threshold", 0.0001)),  # 0.01% — extremely rare
    }


# ---------------------------------------------------------------------------
# Pathogenic criteria
# ---------------------------------------------------------------------------

def eval_PVS1(variant: Variant, settings: dict) -> bool:
    """
    PVS1 — Pathogenic Very Strong.

    Fires when BOTH conditions are met:
      1. The primary transcript consequence is a canonical null-variant type
         (stop_gained, frameshift_variant, splice_acceptor_variant,
         splice_donor_variant, start_lost).
      2. LOFTEE classifies the variant as high-confidence (HC) LoF.

    Requiring LOFTEE HC prevents PVS1 from firing when:
      - LoF is not the established disease mechanism for the gene
      - The variant falls in a low-complexity or end-truncation region
      - LOFTEE assigns low confidence (LC) due to NMD escape, short transcript, etc.

    Without LOFTEE HC confirmation, consequence alone is insufficient for PVS1.
    The variant may still accumulate supporting evidence through PP3, PM2, etc.

    Data used: variant.consequence, variant.loftee_lof
    """
    if not variant.consequence:
        return False
    csq_parts = {c.strip() for c in variant.consequence.split("&") if c.strip()}
    if not (csq_parts & _PVS1_CONSEQUENCES):
        return False
    return variant.loftee_lof == "HC"


def eval_PS3(variant: Variant, settings: dict) -> bool:
    """
    PS3 (ClinVar proxy) — Pathogenic Strong.

    Fires when:
      1. ClinVar significance is Pathogenic or Likely_pathogenic (not Conflicting).
      2. ClinVar review status corresponds to ≥ 2 stars (criteria provided,
         multiple submitters, no conflicts — or better).

    Data used: variant.clinvar_significance, variant.clinvar_review_status

    Limitation (v1): this is a proxy for "well-established functional studies"
    (the original PS3 criterion).  ClinVar review status is used as an indicator
    of evidence quality.  Strong pathogenic ClinVar entries with only 1-star
    review are NOT counted (too little submitter consensus to be "well-established").
    """
    if not _is_clinvar_p_or_lp(variant.clinvar_significance):
        return False
    stars = _clinvar_review_stars(variant.clinvar_review_status)
    return stars >= 2


def eval_PM2(variant: Variant, settings: dict) -> bool:
    """
    PM2 — Pathogenic Moderate.

    Fires when:
      - variant.gnomad_af is None (variant not observed in population databases), OR
      - variant.gnomad_af < pm2_af_threshold (extremely rare; default 0.0001 / 0.01%)

    Data used: variant.gnomad_af
      (set by parser as MAX_AF primary / gnomAD_AF fallback)

    Threshold: read from settings as acmg.pm2_af_threshold (default 0.0001).
    This is deliberately more stringent than the variant filter threshold (max_af)
    to distinguish "rare enough to pass filtering" from "extremely rare / novel"
    as required by PM2.
    """
    af_t = _af_thresholds(settings)
    if variant.gnomad_af is None:
        return True   # absent from population databases
    return variant.gnomad_af < af_t["pm2"]


def eval_PM4(variant: Variant, settings: dict) -> bool:
    """
    PM4 — Pathogenic Moderate.

    Fires when the primary transcript consequence is an in-frame insertion,
    in-frame deletion, or stop-lost variant — all of which change protein length
    without disrupting the reading frame (i.e., not frameshift).

    Data used: variant.consequence (primary transcript, split on '&')

    Note: frameshift_variant triggers PVS1 (more severe), not PM4.
    stop_lost produces protein elongation (read-through) and is a protein-length
    change, hence PM4 not PVS1.
    """
    if not variant.consequence:
        return False
    csq_parts = {c.strip() for c in variant.consequence.split("&") if c.strip()}
    return bool(csq_parts & _PM4_CONSEQUENCES)


def eval_PP3(variant: Variant, settings: dict) -> bool:
    """
    PP3 — Pathogenic Supporting.

    Fires when ≥ 2 INDEPENDENT in-silico signals predict a damaging effect.

    Anti-correlation grouping (SOP requirement):
      CADD and REVEL share training data and are treated as ONE correlated
      group (Group 1). They contribute a single damaging vote only when at
      least one predicts damaging AND none in the group predicts benign
      (i.e., no intra-group conflict). CADD alone or REVEL alone still
      contribute one vote if the other tool has no data.

      SpliceAI (Group 2) and AlphaMissense (Group 3) are independent and
      each contribute their own vote.

    PP3 fires when ≥ 2 of the 3 independent groups vote damaging.
    This prevents CADD+REVEL alone (correlated) from triggering PP3.

    Data used: variant.cadd_phred, variant.revel_score,
               variant.spliceai_max_score, variant.alphamissense_score
    Thresholds: read from settings.acmg.insilico_thresholds.
    """
    t = _insilico_thresholds(settings)

    # Group 1 — correlated: CADD + REVEL (count as ONE vote)
    g1_dam, g1_ben = False, False
    if variant.cadd_phred is not None:
        if variant.cadd_phred >= t["cadd_pathogenic"]:
            g1_dam = True
        else:
            g1_ben = True
    if variant.revel_score is not None:
        if variant.revel_score >= t["revel_pathogenic"]:
            g1_dam = True
        else:
            g1_ben = True
    # Intra-group conflict cancels the vote
    group1_damaging = g1_dam and not g1_ben

    # Group 2 — independent: SpliceAI
    group2_damaging = (
        variant.spliceai_max_score is not None
        and variant.spliceai_max_score >= t["spliceai_high"]
    )

    # Group 3 — independent: AlphaMissense
    group3_damaging = (
        variant.alphamissense_score is not None
        and variant.alphamissense_score >= t["am_pathogenic"]
    )

    return sum([group1_damaging, group2_damaging, group3_damaging]) >= 2


# ---------------------------------------------------------------------------
# Benign criteria
# ---------------------------------------------------------------------------

def eval_BA1(variant: Variant, settings: dict) -> bool:
    """
    BA1 — Benign Stand-Alone.

    Fires when gnomAD AF > 5 % (allele frequency above the stand-alone benign
    threshold).  A single common variant cannot be disease-causing for a rare
    Mendelian disorder.

    Data used: variant.gnomad_af
    Threshold: read from settings as acmg.ba1_af_threshold (default 0.05).
    """
    if variant.gnomad_af is None:
        return False
    af_t = _af_thresholds(settings)
    return variant.gnomad_af > af_t["ba1"]


def eval_BS1(variant: Variant, settings: dict) -> bool:
    """
    BS1 — Benign Strong.

    Fires when gnomAD AF exceeds the disorder-specific expected frequency
    (typically much lower than 5 %).

    Data used: variant.gnomad_af
    Threshold: read from settings as acmg.bs1_af_threshold (default 0.01 / 1 %).

    Limitation (v1): without a gene- or disorder-specific expected frequency,
    a fixed threshold is used.  This is a simplification — proper BS1 evaluation
    requires knowing the maximum credible allele frequency for the disorder.
    """
    if variant.gnomad_af is None:
        return False
    af_t = _af_thresholds(settings)
    return variant.gnomad_af > af_t["bs1"]


def eval_BP4(variant: Variant, settings: dict) -> bool:
    """
    BP4 — Benign Supporting.

    Fires when ≥ 2 INDEPENDENT in-silico signals predict a benign effect.

    Anti-correlation grouping (SOP requirement, mirrors PP3):
      CADD and REVEL share training data and are treated as ONE correlated
      group (Group 1). They contribute a single benign vote only when at
      least one predicts benign AND none in the group predicts damaging.
      Intra-group conflict (e.g., CADD benign + REVEL damaging) cancels the vote.

      REVEL uses a three-zone model: ≥ revel_pathogenic → damaging,
      ≤ revel_benign → benign, in between → no signal (gray zone).
      CADD uses a binary cutoff: ≥ cadd_pathogenic → damaging, else benign.

      SpliceAI (Group 2) and AlphaMissense (Group 3) are independent and
      each contribute their own vote.

    BP4 fires when ≥ 2 of the 3 independent groups vote benign.
    This prevents CADD+REVEL alone (correlated) from triggering BP4.

    Data used: variant.cadd_phred, variant.revel_score,
               variant.spliceai_max_score, variant.alphamissense_score
    Thresholds: read from settings.acmg.insilico_thresholds.
    """
    t = _insilico_thresholds(settings)

    # Group 1 — correlated: CADD + REVEL (count as ONE vote)
    g1_dam, g1_ben = False, False
    if variant.cadd_phred is not None:
        if variant.cadd_phred >= t["cadd_pathogenic"]:
            g1_dam = True
        else:
            g1_ben = True
    if variant.revel_score is not None:
        if variant.revel_score >= t["revel_pathogenic"]:
            g1_dam = True
        elif variant.revel_score <= t["revel_benign"]:
            g1_ben = True
        # Gray zone (revel_benign < score < revel_pathogenic): no signal
    # Intra-group conflict cancels the vote
    group1_benign = g1_ben and not g1_dam

    # Group 2 — independent: SpliceAI
    group2_benign = (
        variant.spliceai_max_score is not None
        and variant.spliceai_max_score <= t["spliceai_low"]
    )

    # Group 3 — independent: AlphaMissense
    group3_benign = (
        variant.alphamissense_score is not None
        and variant.alphamissense_score <= t["am_benign"]
    )

    return sum([group1_benign, group2_benign, group3_benign]) >= 2


def eval_BP7(variant: Variant, settings: dict) -> bool:
    """
    BP7 — Benign Supporting.

    Fires when:
      1. The primary transcript consequence is synonymous_variant (or includes
         synonymous_variant in a compound consequence string), AND
      2. SpliceAI max delta score is either absent (not run) or below the
         spliceai_low threshold (default 0.1), indicating no predicted splice impact.

    Data used: variant.consequence, variant.spliceai_max_score

    Rationale for allowing absent SpliceAI: if no splice prediction tool was run,
    there is no positive evidence of splice disruption.  Combined with a synonymous
    consequence annotation (which implies the variant does not alter the reading
    frame or amino acid sequence), absence of splice evidence supports BP7.
    Users running the pipeline with SpliceAI data will get more precise BP7 calls.
    """
    if not variant.consequence:
        return False
    csq_parts = {c.strip() for c in variant.consequence.split("&") if c.strip()}
    if not csq_parts & _BP7_CONSEQUENCES:
        return False
    t = _insilico_thresholds(settings)
    # If SpliceAI was run, require score to be low; if absent, allow BP7
    if variant.spliceai_max_score is not None:
        return variant.spliceai_max_score <= t["spliceai_low"]
    return True   # no SpliceAI data: no positive evidence of splice impact


# ---------------------------------------------------------------------------
# Catalogue of all implemented v1 rule functions
# Used by classifier.py to iterate over active rules.
# ---------------------------------------------------------------------------

#: Maps ACMG code → evaluation function for v1-implemented rules only.
#: Add new rule functions here when implementing additional criteria.
IMPLEMENTED_RULES = {
    "PVS1": eval_PVS1,
    "PS3":  eval_PS3,
    "PM2":  eval_PM2,
    "PM4":  eval_PM4,
    "PP3":  eval_PP3,
    "BA1":  eval_BA1,
    "BS1":  eval_BS1,
    "BP4":  eval_BP4,
    "BP7":  eval_BP7,
}

#: ACMG codes present in the evidence catalogue with implemented_v1=True but
#: NOT yet evaluable in v1 due to missing data sources.  Listed here for
#: transparency — callers should not assume these will ever be triggered.
UNIMPLEMENTED_V1_RULES = {
    "PS1": "Requires per-amino-acid-position ClinVar database lookup",
    "PM1": "Requires curated mutational-hotspot / functional-domain database",
    "PM5": "Requires per-codon ClinVar lookup (novel missense at pathogenic codon)",
    "PP2": "Requires per-gene missense-constraint data (gnomAD pLoF metrics)",
    "BS2": "Requires healthy-adult zygosity counts (gnomAD homozygous by phenotype)",
}
