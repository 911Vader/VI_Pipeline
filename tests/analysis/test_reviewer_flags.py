"""
Tests for analysis/reviewer_flags.py

Strategy:
  - Tests operate on minimal Variant fixtures.
  - Each test verifies that specific flags are (or are not) generated.
  - Settings are built inline to match config/settings.yaml structure.

Run: pytest tests/analysis/test_reviewer_flags.py -v
"""

import pytest
from analysis.parser import Variant
from analysis.reviewer_flags import annotate_reviewer_flags, _build_flags, _thresholds


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    """Minimal Variant with safe defaults."""
    defaults = dict(
        chrom="chr1", pos=100, ref="A", alt="T",
        filter_status="PASS",
        consequence="missense_variant",
        gene_symbol="BRCA1",
        gnomad_af=None,
        cadd_phred=None,
        revel_score=None,
        spliceai_max_score=None,
        alphamissense_score=None,
        clinvar_significance=None,
        clinvar_review_status=None,
        loftee_lof=None,
        loftee_filter=None,
        excluded=False,
        deprioritized=False,
        deprioritization_reason=None,
        in_gene_list=False,
        acmg_classification=None,
        acmg_evidence_codes=[],
        depth=None,
        genotype_quality=None,
        allele_depth=None,
        vaf=None,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def _make_settings(**overrides) -> dict:
    """Minimal settings dict."""
    base = {
        "acmg": {
            "insilico_thresholds": {
                "cadd_phred_pathogenic": 20,
                "revel_pathogenic": 0.75,
                "spliceai_high": 0.5,
            },
            "ba1_af_threshold": 0.05,
        },
        "filtering": {
            "max_af": 0.01,
        },
        "reviewer_flags": {
            "cadd_high": 20,
            "cadd_very_high": 30,
            "revel_high": 0.75,
            "spliceai_high": 0.5,
            "spliceai_moderate": 0.2,
            "low_depth": 20,
            "low_gq": 20,
            "low_vaf": 0.20,
        },
    }
    base.update(overrides)
    return base


SETTINGS = _make_settings()


# ===========================================================================
# annotate_reviewer_flags — entry point
# ===========================================================================

class TestAnnotateReviewerFlags:
    def test_populates_reviewer_flags_for_all_variants(self):
        """annotate_reviewer_flags populates reviewer_flags for every variant."""
        variants = [_make_variant(), _make_variant(gnomad_af=0.001)]
        annotate_reviewer_flags(variants, SETTINGS)
        for v in variants:
            assert isinstance(v.reviewer_flags, list)

    def test_returns_same_list(self):
        """Function returns the same list object (in-place mutation)."""
        variants = [_make_variant()]
        result = annotate_reviewer_flags(variants, SETTINGS)
        assert result is variants

    def test_empty_variant_list(self):
        result = annotate_reviewer_flags([], SETTINGS)
        assert result == []

    def test_gene_list_passed_through(self):
        """Gene list parameter is used to set gene-list flags."""
        v = _make_variant(gene_symbol="BRCA1", in_gene_list=True)
        annotate_reviewer_flags([v], SETTINGS, gene_list=["BRCA1"])
        assert any("gene list" in f.lower() for f in v.reviewer_flags)


# ===========================================================================
# Population frequency flags
# ===========================================================================

class TestGnomADFlags:
    def test_absent_from_gnomad(self):
        v = _make_variant(gnomad_af=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Absent from gnomAD" in f for f in flags)

    def test_very_rare_gnomad(self):
        v = _make_variant(gnomad_af=0.000005)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Very rare in gnomAD" in f for f in flags)

    def test_rare_gnomad(self):
        v = _make_variant(gnomad_af=0.005)   # < 0.01 max_af
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Rare in gnomAD" in f for f in flags)

    def test_intermediate_af(self):
        v = _make_variant(gnomad_af=0.02)   # between max_af (0.01) and BA1 (0.05)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Intermediate AF" in f for f in flags)

    def test_common_variant(self):
        v = _make_variant(gnomad_af=0.10)   # > BA1 threshold
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Common variant in gnomAD" in f for f in flags)


# ===========================================================================
# Gene list flags
# ===========================================================================

class TestGeneListFlags:
    def test_in_gene_list_flag(self):
        v = _make_variant(in_gene_list=True)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, {"BRCA1"})
        assert "In supplied gene list" in flags

    def test_outside_gene_list_flag(self):
        v = _make_variant(in_gene_list=False)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, {"BRCA1"})
        assert "Outside supplied gene list" in flags

    def test_no_gene_set_no_gene_flag(self):
        v = _make_variant(in_gene_list=False)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert "In supplied gene list" not in flags
        assert "Outside supplied gene list" not in flags


# ===========================================================================
# CADD flags
# ===========================================================================

class TestCADDFlags:
    def test_very_high_cadd_flag(self):
        v = _make_variant(cadd_phred=35.0)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Very high CADD" in f for f in flags)

    def test_high_cadd_flag(self):
        v = _make_variant(cadd_phred=25.0)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("High CADD" in f for f in flags)
        assert not any("Very high CADD" in f for f in flags)

    def test_below_cadd_threshold_no_flag(self):
        v = _make_variant(cadd_phred=10.0)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("CADD" in f for f in flags)

    def test_none_cadd_no_flag(self):
        v = _make_variant(cadd_phred=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("CADD" in f for f in flags)


# ===========================================================================
# LOFTEE flags
# ===========================================================================

class TestLOFTEEFlags:
    def test_loftee_hc_flag(self):
        v = _make_variant(loftee_lof="HC")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert "LOFTEE high-confidence LoF" in flags

    def test_loftee_lc_flag(self):
        v = _make_variant(loftee_lof="LC")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("LOFTEE low-confidence LoF" in f for f in flags)

    def test_loftee_lc_with_filter_reason(self):
        v = _make_variant(loftee_lof="LC", loftee_filter="END_TRUNC")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("END_TRUNC" in f for f in flags)

    def test_no_loftee_no_flag(self):
        v = _make_variant(loftee_lof=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert "LOFTEE high-confidence LoF" not in flags

    def test_lof_consequence_without_loftee_hc(self):
        """LoF consequence without LOFTEE HC → explicit warning flag."""
        v = _make_variant(consequence="stop_gained", loftee_lof=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("LoF consequence" in f and "LOFTEE HC not confirmed" in f for f in flags)

    def test_lof_consequence_with_loftee_hc_no_warning(self):
        """LoF consequence with LOFTEE HC → no warning flag (HC confirmed)."""
        v = _make_variant(consequence="stop_gained", loftee_lof="HC")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("LOFTEE HC not confirmed" in f for f in flags)


# ===========================================================================
# QC flags
# ===========================================================================

class TestQCFlags:
    def test_low_depth_flag(self):
        v = _make_variant(depth=10)   # < 20 threshold
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Low read depth" in f for f in flags)

    def test_depth_at_threshold_no_flag(self):
        v = _make_variant(depth=20)   # == threshold, not below
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("Low read depth" in f for f in flags)

    def test_none_depth_no_flag(self):
        v = _make_variant(depth=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("Low read depth" in f for f in flags)

    def test_low_gq_flag(self):
        v = _make_variant(genotype_quality=10)   # < 20 threshold
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Low genotype quality" in f for f in flags)

    def test_low_vaf_flag(self):
        v = _make_variant(vaf=0.10)   # < 0.20 threshold
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Low VAF" in f for f in flags)

    def test_sufficient_vaf_no_flag(self):
        v = _make_variant(vaf=0.45)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("Low VAF" in f for f in flags)

    def test_none_vaf_no_flag(self):
        v = _make_variant(vaf=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("Low VAF" in f for f in flags)


# ===========================================================================
# ClinVar flags
# ===========================================================================

class TestClinVarFlags:
    def test_pathogenic_clinvar_flag(self):
        v = _make_variant(clinvar_significance="Pathogenic", clinvar_review_status="reviewed_by_expert_panel")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("ClinVar: Pathogenic" in f for f in flags)

    def test_benign_clinvar_flag(self):
        v = _make_variant(clinvar_significance="Benign")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("ClinVar: Benign" in f for f in flags)

    def test_conflicting_clinvar_flag(self):
        v = _make_variant(clinvar_significance="Conflicting_interpretations_of_pathogenicity")
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("conflicting interpretations" in f.lower() for f in flags)

    def test_no_clinvar_no_flag(self):
        v = _make_variant(clinvar_significance=None)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("ClinVar:" in f for f in flags)


# ===========================================================================
# Deprioritization flag
# ===========================================================================

class TestDeprioritizationFlag:
    def test_deprioritized_variant_has_flag(self):
        v = _make_variant(
            deprioritized=True,
            deprioritization_reason="gnomAD AF 0.03 > threshold 0.01 (above common-variant cutoff; deprioritized for review)",
        )
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert any("Deprioritized" in f for f in flags)

    def test_not_deprioritized_no_flag(self):
        v = _make_variant(deprioritized=False)
        t = _thresholds(SETTINGS)
        flags = _build_flags(v, t, set())
        assert not any("Deprioritized:" in f for f in flags)


# ===========================================================================
# Integration — annotate_reviewer_flags with full workflow
# ===========================================================================

class TestIntegration:
    def test_high_cadd_flag_set_after_annotate(self):
        v = _make_variant(cadd_phred=35.0, gnomad_af=None)
        annotate_reviewer_flags([v], SETTINGS)
        assert any("Very high CADD" in f for f in v.reviewer_flags)

    def test_absent_gnomad_flag_set_after_annotate(self):
        v = _make_variant(gnomad_af=None)
        annotate_reviewer_flags([v], SETTINGS)
        assert any("Absent from gnomAD" in f for f in v.reviewer_flags)

    def test_all_none_variant_no_error(self):
        """A variant with all optional fields=None must not raise."""
        v = _make_variant()
        try:
            annotate_reviewer_flags([v], SETTINGS)
        except Exception as exc:
            pytest.fail(f"annotate_reviewer_flags raised on all-None variant: {exc}")
        assert isinstance(v.reviewer_flags, list)
