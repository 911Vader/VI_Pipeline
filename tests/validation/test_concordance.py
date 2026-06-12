"""Tests for validation/concordance.py"""

import pytest
from dataclasses import dataclass
from typing import List, Optional

from validation.concordance import (
    classify_concordance,
    ConcordanceEngine,
    ValidationMetrics,
    CONCORDANT, BORDERLINE, DIRECTION_ERROR, CRITICAL_FAILURE, NOT_FOUND,
)
from validation.truth_set import TruthEntry, build_truth_index
from analysis.parser import Variant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_truth(vid, gene, expected):
    return TruthEntry(variant_id=vid, gene=gene, expected_classification=expected)


def _make_variant(chrom, pos, ref, alt, classification, codes=None):
    v = Variant(chrom=chrom, pos=pos, ref=ref, alt=alt)
    v.acmg_classification = classification
    v.acmg_evidence_codes = codes or []
    v.reviewer_flags = []
    return v


# ---------------------------------------------------------------------------
# classify_concordance
# ---------------------------------------------------------------------------

class TestClassifyConcordance:
    def test_exact_match_is_concordant(self):
        for cls in ("Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"):
            assert classify_concordance(cls, cls) == CONCORDANT

    def test_pathogenic_called_benign_is_critical(self):
        assert classify_concordance("Pathogenic", "Benign") == CRITICAL_FAILURE

    def test_benign_called_pathogenic_is_critical(self):
        assert classify_concordance("Benign", "Pathogenic") == CRITICAL_FAILURE

    def test_pathogenic_called_likely_benign_is_direction_error(self):
        assert classify_concordance("Pathogenic", "Likely Benign") == DIRECTION_ERROR

    def test_likely_pathogenic_called_benign_is_direction_error(self):
        assert classify_concordance("Likely Pathogenic", "Benign") == DIRECTION_ERROR

    def test_benign_called_likely_pathogenic_is_direction_error(self):
        assert classify_concordance("Benign", "Likely Pathogenic") == DIRECTION_ERROR

    def test_vus_called_likely_pathogenic_is_borderline(self):
        assert classify_concordance("VUS", "Likely Pathogenic") == BORDERLINE

    def test_vus_called_likely_benign_is_borderline(self):
        assert classify_concordance("VUS", "Likely Benign") == BORDERLINE

    def test_likely_pathogenic_called_vus_is_borderline(self):
        assert classify_concordance("Likely Pathogenic", "VUS") == BORDERLINE

    def test_likely_benign_called_vus_is_borderline(self):
        assert classify_concordance("Likely Benign", "VUS") == BORDERLINE

    def test_pathogenic_called_likely_pathogenic_is_borderline(self):
        assert classify_concordance("Pathogenic", "Likely Pathogenic") == BORDERLINE

    def test_benign_called_likely_benign_is_borderline(self):
        assert classify_concordance("Benign", "Likely Benign") == BORDERLINE

    def test_pathogenic_called_vus_is_direction_error(self):
        assert classify_concordance("Pathogenic", "VUS") == DIRECTION_ERROR

    def test_benign_called_vus_is_direction_error(self):
        assert classify_concordance("Benign", "VUS") == DIRECTION_ERROR

    def test_lp_called_lb_is_direction_error(self):
        assert classify_concordance("Likely Pathogenic", "Likely Benign") == DIRECTION_ERROR


# ---------------------------------------------------------------------------
# ConcordanceEngine
# ---------------------------------------------------------------------------

class TestConcordanceEngine:
    def _build(self, truth_entries, variants):
        index = build_truth_index(truth_entries)
        return ConcordanceEngine(index, variants)

    def test_concordant_variant(self):
        te = _make_truth("chr1:100:A:T", "GENE1", "Pathogenic")
        v  = _make_variant("chr1", 100, "A", "T", "Pathogenic")
        engine = self._build([te], [v])
        results, metrics = engine.run()
        assert len(results) == 1
        assert results[0].concordance_status == CONCORDANT
        assert metrics.n_concordant == 1
        assert metrics.n_critical_failure == 0

    def test_critical_failure(self):
        te = _make_truth("chr1:100:A:T", "GENE1", "Pathogenic")
        v  = _make_variant("chr1", 100, "A", "T", "Benign")
        engine = self._build([te], [v])
        results, metrics = engine.run()
        assert results[0].concordance_status == CRITICAL_FAILURE
        assert metrics.n_critical_failure == 1
        assert results[0].is_critical_failure is True

    def test_not_found_when_variant_absent(self):
        te = _make_truth("chr1:999:A:T", "GENE1", "Pathogenic")
        v  = _make_variant("chr1", 100, "A", "T", "Pathogenic")
        engine = self._build([te], [v])
        results, metrics = engine.run()
        assert results[0].concordance_status == NOT_FOUND
        assert metrics.n_not_found == 1
        assert metrics.n_matched == 0

    def test_chrom_normalisation_strips_chr_prefix(self):
        """Truth set 'chr1' should match VCF variant '1' and vice versa."""
        te = _make_truth("chr1:100:A:T", "GENE1", "Pathogenic")
        v  = _make_variant("1", 100, "A", "T", "Pathogenic")  # VCF without 'chr'
        engine = self._build([te], [v])
        results, metrics = engine.run()
        assert results[0].concordance_status == CONCORDANT

    def test_borderline_result(self):
        te = _make_truth("chr1:100:A:T", "GENE1", "Likely Pathogenic")
        v  = _make_variant("chr1", 100, "A", "T", "VUS")
        engine = self._build([te], [v])
        results, metrics = engine.run()
        assert results[0].concordance_status == BORDERLINE
        assert metrics.n_borderline == 1

    def test_multiple_results_metrics(self):
        entries = [
            _make_truth("chr1:100:A:T", "G1", "Pathogenic"),
            _make_truth("chr1:200:A:T", "G2", "Benign"),
            _make_truth("chr1:300:A:T", "G3", "VUS"),
        ]
        variants = [
            _make_variant("chr1", 100, "A", "T", "Pathogenic"),  # concordant
            _make_variant("chr1", 200, "A", "T", "Pathogenic"),  # critical failure
            _make_variant("chr1", 300, "A", "T", "VUS"),          # concordant
        ]
        engine = self._build(entries, variants)
        results, metrics = engine.run()
        assert metrics.n_concordant == 2
        assert metrics.n_critical_failure == 1
        assert metrics.n_matched == 3

    def test_sensitivity_calculation(self):
        entries = [
            _make_truth("chr1:100:A:T", "G1", "Pathogenic"),
            _make_truth("chr1:200:A:T", "G2", "Likely Pathogenic"),
        ]
        variants = [
            _make_variant("chr1", 100, "A", "T", "Pathogenic"),        # TP
            _make_variant("chr1", 200, "A", "T", "VUS"),               # FN
        ]
        engine = self._build(entries, variants)
        _, metrics = engine.run()
        assert metrics.tp_pathogenic == 1
        assert metrics.fn_pathogenic == 1
        assert metrics.sensitivity == pytest.approx(0.5)

    def test_specificity_calculation(self):
        entries = [
            _make_truth("chr1:100:A:T", "G1", "Benign"),
            _make_truth("chr1:200:A:T", "G2", "Likely Benign"),
        ]
        variants = [
            _make_variant("chr1", 100, "A", "T", "Benign"),            # TN
            _make_variant("chr1", 200, "A", "T", "Pathogenic"),        # FP
        ]
        engine = self._build(entries, variants)
        _, metrics = engine.run()
        assert metrics.tn_benign == 1
        assert metrics.fp_benign == 1
        assert metrics.specificity == pytest.approx(0.5)

    def test_concordance_rate_all_concordant(self):
        entries = [_make_truth(f"chr1:{i}:A:T", f"G{i}", "VUS") for i in range(1, 6)]
        variants = [_make_variant("chr1", i, "A", "T", "VUS") for i in range(1, 6)]
        _, metrics = self._build(entries, variants).run()
        assert metrics.concordance_rate == pytest.approx(1.0)

    def test_vus_rate(self):
        entries = [_make_truth(f"chr1:{i}:A:T", f"G{i}", "VUS") for i in range(1, 5)]
        variants = [_make_variant("chr1", i, "A", "T", "VUS") for i in range(1, 5)]
        _, metrics = self._build(entries, variants).run()
        assert metrics.vus_rate == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# ValidationMetrics.pass_fail_verdict
# ---------------------------------------------------------------------------

class TestPassFailVerdict:
    def _metrics(self, n_matched, n_concordant, n_critical):
        m = ValidationMetrics(
            n_truth_total=n_matched,
            n_matched=n_matched,
            n_concordant=n_concordant,
            n_critical_failure=n_critical,
        )
        return m

    def test_perfect_concordance_is_pass(self):
        m = self._metrics(10, 10, 0)
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "PASS"

    def test_ninety_five_pct_concordance_is_pass(self):
        m = self._metrics(20, 19, 0)  # 19/20 = 95%
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "PASS"

    def test_ninety_pct_is_conditional_pass(self):
        m = self._metrics(10, 9, 0)   # 9/10 = 90%
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "CONDITIONAL PASS"

    def test_below_ninety_pct_is_fail(self):
        m = self._metrics(10, 8, 0)   # 8/10 = 80%
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "FAIL"

    def test_any_critical_failure_is_fail(self):
        m = self._metrics(10, 10, 1)  # 100% concordance but 1 critical failure
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "FAIL"

    def test_fail_reason_mentions_critical_failure(self):
        m = self._metrics(10, 10, 1)
        _, reasons = m.pass_fail_verdict()
        assert any("CRITICAL" in r.upper() for r in reasons)

    def test_no_matched_variants_is_fail(self):
        """Zero matched variants means no VCF variants found — that is a FAIL."""
        m = ValidationMetrics(n_truth_total=5, n_matched=0, n_not_found=5)
        verdict, _ = m.pass_fail_verdict()
        assert verdict == "FAIL"
