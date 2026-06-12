"""
discordance_log.py — Build structured discordance log entries.

One entry is created for every non-concordant VariantResult
(Borderline, Direction Error, Critical Failure, Not Found).

Each entry captures the full evidence trail required by the validation SOP:
  - Variant identity (HGVS, gene)
  - Expected vs pipeline classification
  - All ACMG codes applied
  - All tool outputs (gnomAD global + subpopulation AFs, ClinVar, in-silico)
  - Conflict type and root cause hypothesis
  - Required action
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

from analysis.parser import Variant
from validation.concordance import (
    VariantResult,
    CONCORDANT, BORDERLINE, DIRECTION_ERROR, CRITICAL_FAILURE, NOT_FOUND,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Conflict type → action mapping
# ---------------------------------------------------------------------------

_CONFLICT_ACTIONS = {
    CRITICAL_FAILURE: "Pipeline fix — investigate immediately; block go-live until resolved",
    DIRECTION_ERROR:  "Pipeline fix or curator review — root cause must be documented",
    BORDERLINE:       "Curator review — VUS boundary case; flag for expert review",
    NOT_FOUND:        "Pipeline fix — variant absent from VCF; investigate annotation step",
}

_ROOT_CAUSE_HINTS = {
    CRITICAL_FAILURE: (
        "Possible causes: incorrect AF filter boundary; wrong ClinVar source; "
        "BA1 firing incorrectly; PVS1 not triggering when required"
    ),
    DIRECTION_ERROR: (
        "Possible causes: missing evidence code (annotation gap); "
        "threshold misconfiguration; conflicting tool predictions suppressing PP3/BP4; "
        "ClinVar significance not meeting star threshold"
    ),
    BORDERLINE: (
        "Borderline VUS/LP or VUS/LB: one additional piece of evidence may resolve; "
        "check ClinVar star rating, gnomAD subpopulation AFs, and LOFTEE status"
    ),
    NOT_FOUND: (
        "Variant absent from VCF: check AF filter (may have been excluded); "
        "check consequence tier; check if variant is in excluded tier; "
        "verify VEP annotation completed successfully for this locus"
    ),
}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class DiscordanceEntry:
    """Structured log entry for one discordant variant."""
    # Identity
    variant_hgvs: str          # Variant.hgvs_c or variant_id fallback
    gene: str
    variant_id: str            # chr:pos:ref:alt

    # Classification
    expected: str
    pipeline_output: str
    criteria_applied: str      # space-separated ACMG codes

    # Tool outputs
    gnomad_af_global: str
    clinvar_significance: str
    clinvar_review_status: str
    cadd: str
    revel: str
    spliceai: str
    alphamissense: str
    loftee: str

    # Concordance
    conflict_type: str
    root_cause: str
    action_required: str

    # Optional extras
    curator_notes: str = ""
    expected_acmg_codes: str = ""


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

def build_discordance_entries(
    results: List[VariantResult],
) -> List[DiscordanceEntry]:
    """
    Filter results to non-concordant outcomes and build DiscordanceEntry objects.
    """
    entries: List[DiscordanceEntry] = []

    for r in results:
        if r.concordance_status == CONCORDANT:
            continue

        v = r.variant
        entry = DiscordanceEntry(
            variant_hgvs        = _hgvs_or_id(v, r.truth_entry.variant_id),
            gene                = r.gene,
            variant_id          = r.variant_id,
            expected            = r.expected,
            pipeline_output     = r.pipeline_classification,
            criteria_applied    = " ".join(r.pipeline_acmg_codes) if r.pipeline_acmg_codes else "—",
            gnomad_af_global    = _fmt_af(v.gnomad_af if v else None),
            clinvar_significance= (v.clinvar_significance or "—") if v else "—",
            clinvar_review_status=(v.clinvar_review_status or "—") if v else "—",
            cadd                = _fmt_score(v.cadd_phred          if v else None, ".1f"),
            revel               = _fmt_score(v.revel_score          if v else None, ".3f"),
            spliceai            = _fmt_score(v.spliceai_max_score   if v else None, ".3f"),
            alphamissense       = _fmt_score(v.alphamissense_score  if v else None, ".3f"),
            loftee              = _loftee(v),
            conflict_type       = r.concordance_status,
            root_cause          = _ROOT_CAUSE_HINTS.get(r.concordance_status, "Unknown"),
            action_required     = _CONFLICT_ACTIONS.get(r.concordance_status, "Review required"),
            curator_notes       = r.truth_entry.curator_notes,
            expected_acmg_codes = " ".join(r.truth_entry.expected_acmg_codes),
        )
        entries.append(entry)

        logger.info(
            "Discordance [%s] %s: expected=%s pipeline=%s codes=%s",
            r.concordance_status, r.variant_id,
            r.expected, r.pipeline_classification,
            " ".join(r.pipeline_acmg_codes) or "none",
        )

    return entries


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _hgvs_or_id(v: Optional[Variant], fallback: str) -> str:
    if v is None:
        return fallback
    return v.hgvs_c or fallback


def _fmt_af(af: Optional[float]) -> str:
    if af is None:
        return "absent"
    if af == 0.0:
        return "0"
    return f"{af:.6g}"


def _fmt_score(score: Optional[float], fmt: str) -> str:
    if score is None:
        return "—"
    return format(score, fmt)


def _loftee(v: Optional[Variant]) -> str:
    if v is None:
        return "—"
    lof = v.loftee_lof or "—"
    filt = f" ({v.loftee_filter})" if v.loftee_filter else ""
    return f"{lof}{filt}"
