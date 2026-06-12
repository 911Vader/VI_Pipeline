"""
Tests for analysis/acmg/classifier.py and evidence_codes.classify()

Strategy:
  - classify() (evidence_codes.py) is tested directly for all combination rules
    from Richards 2015 Table 5.
  - classify_variant() and classify_all() are tested with constructed Variants.
  - Each ACMG classification category is covered by at least one combination test.
  - Missing-data and edge-case behavior is verified.
  - Tests confirm that unimplemented rules never produce codes.

Run: pytest tests/analysis/acmg/test_classifier.py -v
"""

import pytest
from analysis.parser import Variant
from analysis.acmg.evidence_codes import (
    classify,
    PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN,
)
from analysis.acmg.classifier import (
    classify_variant,
    classify_all,
    get_classification_summary,
)
from analysis.acmg.rules import IMPLEMENTED_RULES, UNIMPLEMENTED_V1_RULES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    defaults = dict(
        chrom="chr1", pos=100, ref="A", alt="T",
        filter_status="PASS",
        consequence="",
        gene_symbol="BRCA1",
        gnomad_af=None,
        cadd_phred=None,
        revel_score=None,
        spliceai_max_score=None,
        clinvar_significance=None,
        clinvar_review_status=None,
        excluded=False,
        in_gene_list=False,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


SETTINGS = {
    "acmg": {
        "insilico_thresholds": {
            "cadd_phred_pathogenic": 20,
            "revel_pathogenic": 0.75,
            "revel_benign": 0.15,
            "spliceai_high": 0.5,
            "spliceai_low": 0.1,
        },
        "ba1_af_threshold": 0.05,
        "bs1_af_threshold": 0.01,
        "pm2_af_threshold": 0.0001,
    }
}


# ===========================================================================
# classify() — Richards 2015 Table 5 combination logic
# ===========================================================================

class TestClassifyPathogenic:
    def test_pvs1_plus_one_ps(self):
        assert classify(["PVS1", "PS3"]) == PATHOGENIC

    def test_pvs1_plus_two_pm(self):
        assert classify(["PVS1", "PM2", "PM4"]) == PATHOGENIC

    def test_pvs1_plus_one_pm_one_pp(self):
        assert classify(["PVS1", "PM2", "PP3"]) == PATHOGENIC

    def test_pvs1_plus_two_pp(self):
        assert classify(["PVS1", "PP3", "BP7"]) != PATHOGENIC  # BP7 is benign, PP3 only 1 PP
        assert classify(["PVS1", "PP3", "PP3"]) == PATHOGENIC   # 2 PP (same code twice)

    def test_two_strong_pathogenic(self):
        assert classify(["PS3", "PS3"]) == PATHOGENIC

    def test_one_ps_three_pm(self):
        assert classify(["PS3", "PM2", "PM4", "PM2"]) == PATHOGENIC

    def test_empty_codes_is_vus(self):
        assert classify([]) == VUS

    def test_unknown_code_ignored(self):
        """Unknown codes are silently ignored; result should still be VUS."""
        assert classify(["FAKE99"]) == VUS


class TestClassifyLikelyPathogenic:
    def test_pvs1_plus_one_pm(self):
        assert classify(["PVS1", "PM2"]) == LIKELY_PATHOGENIC

    def test_one_ps_one_pm(self):
        assert classify(["PS3", "PM2"]) == LIKELY_PATHOGENIC

    def test_one_ps_two_pm(self):
        assert classify(["PS3", "PM2", "PM4"]) == LIKELY_PATHOGENIC

    def test_one_ps_two_pp(self):
        assert classify(["PS3", "PP3", "PP3"]) == LIKELY_PATHOGENIC

    def test_three_pm(self):
        assert classify(["PM2", "PM4", "PM2"]) == LIKELY_PATHOGENIC

    def test_two_pm_two_pp(self):
        assert classify(["PM2", "PM4", "PP3", "PP3"]) == LIKELY_PATHOGENIC

    def test_one_pm_four_pp(self):
        codes = ["PM2"] + ["PP3"] * 4
        assert classify(codes) == LIKELY_PATHOGENIC


class TestClassifyBenign:
    def test_ba1_alone_is_benign(self):
        assert classify(["BA1"]) == BENIGN

    def test_ba1_overrides_pathogenic_evidence(self):
        """BA1 is stand-alone and overrides all pathogenic evidence."""
        assert classify(["BA1", "PVS1", "PS3"]) == BENIGN

    def test_two_strong_benign(self):
        assert classify(["BS1", "BS1"]) == BENIGN

    def test_two_different_benign_strong(self):
        # BS1 and BS1 (counted twice)
        assert classify(["BS1", "BS1"]) == BENIGN


class TestClassifyLikelyBenign:
    def test_one_bs_one_bp(self):
        assert classify(["BS1", "BP4"]) == LIKELY_BENIGN

    def test_one_bs_one_bp7(self):
        assert classify(["BS1", "BP7"]) == LIKELY_BENIGN

    def test_two_supporting_benign(self):
        assert classify(["BP4", "BP7"]) == LIKELY_BENIGN

    def test_three_supporting_benign(self):
        assert classify(["BP4", "BP7", "BP4"]) == LIKELY_BENIGN


class TestClassifyVUS:
    def test_single_moderate_is_vus(self):
        assert classify(["PM2"]) == VUS

    def test_single_supporting_pathogenic_is_vus(self):
        assert classify(["PP3"]) == VUS

    def test_single_supporting_benign_is_vus(self):
        assert classify(["BP4"]) == VUS

    def test_pvs1_alone_is_vus(self):
        """PVS1 alone is not sufficient for Pathogenic or Likely Pathogenic."""
        assert classify(["PVS1"]) == VUS

    def test_conflicting_lp_lb_is_vus(self):
        """Conflicting LP + LB evidence → VUS (conservative)."""
        # LP: PS3 + PM2; LB: BS1 + BP4
        assert classify(["PS3", "PM2", "BS1", "BP4"]) == VUS

    def test_pvs1_plus_one_pp_is_vus(self):
        """PVS1 + 1 PP is not in Table 5 (needs PVS1 + ≥2 PP for Pathogenic)."""
        assert classify(["PVS1", "PP3"]) == VUS


class TestClassifyConflictingEvidence:
    def test_ba1_overrides_pvs1(self):
        assert classify(["BA1", "PVS1"]) == BENIGN

    def test_pathogenic_with_lp_benign_still_pathogenic(self):
        """Pathogenic evidence strong enough to override likely-benign."""
        # PVS1 + 2 PM → Pathogenic; also BP4 + BP7 → Likely Benign
        result = classify(["PVS1", "PM2", "PM4", "BP4", "BP7"])
        assert result == PATHOGENIC


# ===========================================================================
# classify_variant() — integration
# ===========================================================================

class TestClassifyVariant:
    def test_returns_tuple(self):
        v = _make_variant()
        codes, cat, confidence = classify_variant(v, SETTINGS)
        assert isinstance(codes, list)
        assert isinstance(cat, str)
        assert confidence in ("High", "Moderate", "Low")

    def test_null_variant_not_in_db_gets_codes(self):
        """stop_gained + LOFTEE HC + absent from gnomAD → PVS1 + PM2 → Likely Pathogenic."""
        v = _make_variant(consequence="stop_gained", gnomad_af=None, loftee_lof="HC")
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PVS1" in codes
        assert "PM2" in codes
        assert cat == LIKELY_PATHOGENIC

    def test_pvs1_plus_ps3_gives_pathogenic(self):
        """stop_gained + LOFTEE HC + ClinVar P 2-star + absent AF → PVS1 + PS3 + PM2 → Pathogenic."""
        v = _make_variant(
            consequence="stop_gained",
            gnomad_af=None,
            loftee_lof="HC",
            clinvar_significance="Pathogenic",
            clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts",
        )
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PVS1" in codes
        assert "PS3" in codes
        assert "PM2" in codes
        assert cat == PATHOGENIC

    def test_common_variant_gets_ba1(self):
        """gnomad_af > 5% → BA1 → Benign."""
        v = _make_variant(consequence="missense_variant", gnomad_af=0.10)
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "BA1" in codes
        assert cat == BENIGN

    def test_synonymous_no_scores_gets_bp7(self):
        v = _make_variant(consequence="synonymous_variant")
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "BP7" in codes

    def test_two_independent_groups_gets_pp3(self):
        """CADD (Group 1) + SpliceAI (Group 2) = 2 independent votes → PP3 fires."""
        v = _make_variant(
            consequence="missense_variant",
            cadd_phred=25, spliceai_max_score=0.6,
        )
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PP3" in codes

    def test_cadd_revel_alone_does_not_trigger_pp3(self):
        """CADD + REVEL are in the same anti-correlation group → only 1 vote → PP3 does not fire."""
        v = _make_variant(
            consequence="missense_variant",
            cadd_phred=25, revel_score=0.80,
        )
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PP3" not in codes

    def test_no_data_variant_is_vus(self):
        """Variant with no usable annotation → no codes → VUS."""
        v = _make_variant(
            consequence="missense_variant",
            gnomad_af=0.001,          # rare enough to not trigger BS1/BA1/PM2
        )
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert cat == VUS

    def test_inframe_indel_novel_gets_pm4_pm2(self):
        v = _make_variant(consequence="inframe_deletion", gnomad_af=None)
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PM4" in codes
        assert "PM2" in codes

    def test_excluded_variant_still_classified(self):
        """classifier runs regardless of variant.excluded flag."""
        v = _make_variant(consequence="stop_gained", gnomad_af=None, excluded=True, loftee_lof="HC")
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PVS1" in codes

    def test_no_unimplemented_codes_triggered(self):
        """Unimplemented rules (PS1, PM1, PM5, PP2, BS2) must never appear in output."""
        # Use a variant that might naively trigger these rules
        v = _make_variant(
            consequence="missense_variant",
            gnomad_af=0.0001,
            clinvar_significance="Pathogenic",
            clinvar_review_status="reviewed_by_expert_panel",
            cadd_phred=30,
            revel_score=0.90,
        )
        codes, _, _conf = classify_variant(v, SETTINGS)
        for unimplemented in UNIMPLEMENTED_V1_RULES:
            assert unimplemented not in codes, (
                f"{unimplemented} (not implemented in v1) appeared in output: {codes}"
            )


# ===========================================================================
# classify_all() — batch classification
# ===========================================================================

class TestClassifyAll:
    def test_populates_fields_in_place(self):
        variants = [
            _make_variant(consequence="stop_gained", gnomad_af=None),
            _make_variant(consequence="missense_variant", gnomad_af=0.10),
        ]
        result = classify_all(variants, SETTINGS)
        assert result is variants  # same list, mutated in place
        for v in variants:
            assert isinstance(v.acmg_evidence_codes, list)
            assert v.acmg_classification in (
                PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN
            )

    def test_empty_list(self):
        result = classify_all([], SETTINGS)
        assert result == []

    def test_all_variants_get_classification(self):
        variants = [_make_variant() for _ in range(5)]
        classify_all(variants, SETTINGS)
        for v in variants:
            assert v.acmg_classification is not None

    def test_excluded_variants_still_classified(self):
        v = _make_variant(excluded=True)
        classify_all([v], SETTINGS)
        assert v.acmg_classification is not None


# ===========================================================================
# No fabricated evidence
# ===========================================================================

class TestNoFabricatedEvidence:
    """
    Ensure no evidence codes are generated from missing data.
    A completely un-annotated variant (no scores, no ClinVar, no consequence)
    must produce only codes that legitimately fire for those conditions.
    """

    def test_variant_with_no_scores_has_no_pp3_or_bp4(self):
        v = _make_variant(consequence="missense_variant")
        codes, _, _conf = classify_variant(v, SETTINGS)
        assert "PP3" not in codes
        assert "BP4" not in codes

    def test_variant_with_no_clinvar_has_no_ps3(self):
        v = _make_variant(consequence="missense_variant")
        codes, _, _conf = classify_variant(v, SETTINGS)
        assert "PS3" not in codes

    def test_variant_with_common_af_has_no_pm2(self):
        v = _make_variant(consequence="missense_variant", gnomad_af=0.005)
        codes, _, _conf = classify_variant(v, SETTINGS)
        assert "PM2" not in codes

    def test_only_pvs1_no_extras_when_just_null_variant(self):
        """A null variant with LOFTEE HC + common AF should get PVS1 but also BA1 (if > 5%)."""
        v = _make_variant(consequence="stop_gained", gnomad_af=0.10, loftee_lof="HC")
        codes, cat, _ = classify_variant(v, SETTINGS)
        assert "PVS1" in codes
        assert "BA1" in codes
        assert cat == BENIGN   # BA1 overrides PVS1

    def test_one_in_silico_tool_does_not_give_pp3(self):
        """Only CADD available → cannot confirm 'multiple tools' → no PP3."""
        v = _make_variant(consequence="missense_variant", cadd_phred=35)
        codes, _, _conf = classify_variant(v, SETTINGS)
        assert "PP3" not in codes

    def test_absent_clinvar_review_status_does_not_give_ps3(self):
        """P/LP significance but no review status → 0 stars → no PS3."""
        v = _make_variant(
            consequence="missense_variant",
            clinvar_significance="Pathogenic",
            clinvar_review_status=None,
        )
        codes, _, _conf = classify_variant(v, SETTINGS)
        assert "PS3" not in codes


# ===========================================================================
# get_classification_summary
# ===========================================================================

class TestGetClassificationSummary:
    def test_summary_counts(self):
        variants = [
            _make_variant(consequence="stop_gained", gnomad_af=None),  # → LP+
            _make_variant(consequence="missense_variant", gnomad_af=0.10),  # → Benign (BA1)
            _make_variant(consequence="missense_variant", gnomad_af=0.001),  # → VUS
        ]
        classify_all(variants, SETTINGS)
        summary = get_classification_summary(variants)
        assert summary["total"] == 3
        assert summary["unclassified"] == 0
        assert sum(summary["by_category"].values()) == 3

    def test_unclassified_count(self):
        v = _make_variant()
        # Don't call classify_all — acmg_classification remains None
        summary = get_classification_summary([v])
        assert summary["unclassified"] == 1

    def test_empty_list(self):
        summary = get_classification_summary([])
        assert summary["total"] == 0

    def test_all_vus(self):
        variants = [_make_variant(consequence="missense_variant", gnomad_af=0.001)
                    for _ in range(3)]
        classify_all(variants, SETTINGS)
        summary = get_classification_summary(variants)
        assert summary["by_category"].get(VUS, 0) == 3
