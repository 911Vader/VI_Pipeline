"""
concordance.py — Compare pipeline classifications against the truth set.

Concordance categories (per validation SOP):
  Concordant       — pipeline == expected (exact category match)
  Borderline       — one is VUS and other is Likely Pathogenic or Likely Benign
  Direction Error  — pipeline and expected are on opposite sides of the
                     pathogenic/benign boundary (P/LP vs LB/B), but neither
                     is a critical failure
  Critical Failure — Pathogenic called as Benign, or Benign called as Pathogenic
  Not Found        — variant in truth set but absent from VCF

Metrics:
  Concordance rate = Concordant / (total matched)
  Sensitivity      = TP_path / (TP_path + FN_path)  [ability to detect P/LP]
  Specificity      = TN_ben  / (TN_ben  + FP_ben )  [ability to correctly call B/LB]
  VUS rate         = pipeline VUS calls / total matched
  Critical failure rate = critical failures / total matched
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from analysis.parser import Variant
from validation.truth_set import TruthEntry

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Classification tier helpers
# ---------------------------------------------------------------------------

_PATHOGENIC_TIER = {"Pathogenic", "Likely Pathogenic"}
_BENIGN_TIER     = {"Likely Benign", "Benign"}
_VUS_TIER        = {"VUS"}

# Ordered tier severity: 1 = most pathogenic, 5 = most benign
_TIER_RANK = {
    "Pathogenic":        1,
    "Likely Pathogenic": 2,
    "VUS":               3,
    "Likely Benign":     4,
    "Benign":            5,
}

CONCORDANT       = "Concordant"
BORDERLINE       = "Borderline"
DIRECTION_ERROR  = "Direction Error"
CRITICAL_FAILURE = "Critical Failure"
NOT_FOUND        = "Not Found"


def classify_concordance(
    expected: str,
    pipeline: str,
) -> str:
    """
    Return the concordance category for one variant.

    Parameters
    ----------
    expected:  Expected classification from truth set.
    pipeline:  Classification produced by the pipeline.
    """
    if expected == pipeline:
        return CONCORDANT

    # Critical failure: Pathogenic ↔ Benign direction swap (extreme ends)
    if (expected == "Pathogenic" and pipeline == "Benign") or \
       (expected == "Benign"     and pipeline == "Pathogenic"):
        return CRITICAL_FAILURE

    # Direction error: either side of the P/LP ↔ LB/B divide
    if (expected in _PATHOGENIC_TIER and pipeline in _BENIGN_TIER) or \
       (expected in _BENIGN_TIER     and pipeline in _PATHOGENIC_TIER):
        return DIRECTION_ERROR

    # Borderline: one VUS, adjacent on the tier ladder
    exp_rank  = _TIER_RANK.get(expected, 3)
    pipe_rank = _TIER_RANK.get(pipeline, 3)
    if abs(exp_rank - pipe_rank) == 1:
        return BORDERLINE

    # Larger rank gap (e.g. VUS vs Pathogenic/Benign, or LP vs LB)
    return DIRECTION_ERROR


# ---------------------------------------------------------------------------
# Result record
# ---------------------------------------------------------------------------

@dataclass
class VariantResult:
    """Outcome of comparing one truth set entry against the pipeline."""
    truth_entry: TruthEntry
    variant: Optional[Variant]          # None when not found in VCF
    pipeline_classification: str        # "Not Found" when variant is None
    pipeline_acmg_codes: List[str]      = field(default_factory=list)
    concordance_status: str             = NOT_FOUND
    is_critical_failure: bool           = False

    @property
    def variant_id(self) -> str:
        return self.truth_entry.variant_id

    @property
    def gene(self) -> str:
        return self.truth_entry.gene

    @property
    def expected(self) -> str:
        return self.truth_entry.expected_classification


# ---------------------------------------------------------------------------
# Metrics container
# ---------------------------------------------------------------------------

@dataclass
class ValidationMetrics:
    n_truth_total: int          = 0
    n_matched: int              = 0
    n_not_found: int            = 0
    n_concordant: int           = 0
    n_borderline: int           = 0
    n_direction_error: int      = 0
    n_critical_failure: int     = 0
    n_pipeline_vus: int         = 0

    # Sensitivity / specificity components
    tp_pathogenic: int          = 0   # truth P/LP, pipeline P/LP
    fn_pathogenic: int          = 0   # truth P/LP, pipeline not P/LP
    tn_benign: int              = 0   # truth B/LB, pipeline B/LB
    fp_benign: int              = 0   # truth B/LB, pipeline P/LP

    @property
    def concordance_rate(self) -> Optional[float]:
        if self.n_matched == 0:
            return None
        return self.n_concordant / self.n_matched

    @property
    def sensitivity(self) -> Optional[float]:
        denom = self.tp_pathogenic + self.fn_pathogenic
        return self.tp_pathogenic / denom if denom > 0 else None

    @property
    def specificity(self) -> Optional[float]:
        denom = self.tn_benign + self.fp_benign
        return self.tn_benign / denom if denom > 0 else None

    @property
    def vus_rate(self) -> Optional[float]:
        if self.n_matched == 0:
            return None
        return self.n_pipeline_vus / self.n_matched

    @property
    def critical_failure_rate(self) -> Optional[float]:
        if self.n_matched == 0:
            return None
        return self.n_critical_failure / self.n_matched

    def pass_fail_verdict(self) -> Tuple[str, List[str]]:
        """
        Return (verdict, [reason_strings]) per the SOP pass/fail criteria.

        Verdicts:
          PASS             — concordance ≥ 95%, zero critical failures
          CONDITIONAL PASS — concordance 90–94%, no critical failures
          FAIL             — any critical failure OR concordance < 90%
        """
        reasons: List[str] = []
        rate = self.concordance_rate or 0.0

        if self.n_critical_failure > 0:
            reasons.append(
                f"CRITICAL FAILURE: {self.n_critical_failure} Pathogenic ↔ Benign "
                f"direction error(s) — pipeline must be corrected before go-live"
            )

        if rate < 0.90:
            reasons.append(
                f"Concordance {rate:.1%} is below the 90% minimum threshold"
            )

        if reasons:
            return "FAIL", reasons

        if rate < 0.95:
            reasons.append(
                f"Concordance {rate:.1%} is in the conditional range (90–94%). "
                f"Each discordance must be individually documented and reviewed."
            )
            return "CONDITIONAL PASS", reasons

        return "PASS", [f"Concordance {rate:.1%} — all critical checks passed"]


# ---------------------------------------------------------------------------
# Concordance engine
# ---------------------------------------------------------------------------

class ConcordanceEngine:
    """
    Match pipeline-classified variants against the truth set and compute metrics.
    """

    def __init__(
        self,
        truth_index: Dict,
        classified_variants: List[Variant],
    ) -> None:
        self._truth_index    = truth_index
        self._variants_by_key: Dict = {}
        for v in classified_variants:
            key = (_normalise_chrom(v.chrom), v.pos, v.ref, v.alt)
            self._variants_by_key[key] = v

    def run(self) -> Tuple[List[VariantResult], ValidationMetrics]:
        """
        Compare all truth set entries against the classified variants.

        Returns
        -------
        (results, metrics) — results is one VariantResult per truth entry.
        """
        results: List[VariantResult] = []
        m = ValidationMetrics(n_truth_total=len(self._truth_index))

        for truth_key, truth_entry in self._truth_index.items():
            norm_key = (_normalise_chrom(truth_entry.chrom), truth_entry.pos,
                        truth_entry.ref, truth_entry.alt)
            variant = self._variants_by_key.get(norm_key)

            if variant is None:
                results.append(VariantResult(
                    truth_entry           = truth_entry,
                    variant               = None,
                    pipeline_classification = NOT_FOUND,
                    concordance_status    = NOT_FOUND,
                ))
                m.n_not_found += 1
                logger.warning("Truth entry %s not found in VCF", truth_entry.variant_id)
                continue

            m.n_matched += 1
            pipe_cls = variant.acmg_classification or "VUS"
            status   = classify_concordance(truth_entry.expected_classification, pipe_cls)

            result = VariantResult(
                truth_entry           = truth_entry,
                variant               = variant,
                pipeline_classification = pipe_cls,
                pipeline_acmg_codes   = list(variant.acmg_evidence_codes or []),
                concordance_status    = status,
                is_critical_failure   = (status == CRITICAL_FAILURE),
            )
            results.append(result)

            # Accumulate metrics
            if status == CONCORDANT:
                m.n_concordant += 1
            elif status == BORDERLINE:
                m.n_borderline += 1
            elif status == DIRECTION_ERROR:
                m.n_direction_error += 1
            elif status == CRITICAL_FAILURE:
                m.n_critical_failure += 1

            if pipe_cls == "VUS":
                m.n_pipeline_vus += 1

            expected = truth_entry.expected_classification
            if expected in _PATHOGENIC_TIER:
                if pipe_cls in _PATHOGENIC_TIER:
                    m.tp_pathogenic += 1
                else:
                    m.fn_pathogenic += 1
            elif expected in _BENIGN_TIER:
                if pipe_cls in _BENIGN_TIER:
                    m.tn_benign += 1
                elif pipe_cls in _PATHOGENIC_TIER:
                    m.fp_benign += 1

        logger.info(
            "Concordance: %d/%d matched | %d concordant | %d borderline | "
            "%d direction errors | %d critical failures | %d not found",
            m.n_matched, m.n_truth_total,
            m.n_concordant, m.n_borderline,
            m.n_direction_error, m.n_critical_failure,
            m.n_not_found,
        )
        return results, m


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise_chrom(chrom: str) -> str:
    """Strip leading 'chr' for consistent matching."""
    return chrom.lower().removeprefix("chr") if chrom else chrom
