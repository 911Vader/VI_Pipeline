"""
Tests for analysis/filter_engine.py

Strategy:
  - All tests operate on Variant objects constructed directly (no VCF files).
  - Settings dicts are built inline to match config/settings.yaml structure.
  - Each test covers exactly one filtering behavior for clarity.

Coverage:
  - AF threshold: gnomad_af above/below/equal to threshold
  - MAX_AF primary / gnomAD_AF fallback (via variant.gnomad_af field)
  - Consequence-tier filtering: tier1, tier2, tier3, excluded, compound
  - Gene-list matching: case-insensitive exact, partial no-match, absent gene
  - ClinVar-based deprioritization: benign, likely_benign, pathogenic retained, None retained
  - PASS-only behavior: PASS retained, non-PASS excluded, "." excluded
  - Exclusion reason tracking: populated and human-readable
  - Deprioritization: AF between max_af and BA1, tier3, ClinVar benign

Run: pytest tests/analysis/test_filter_engine.py -v
"""

import pytest
from dataclasses import replace

from analysis.parser import Variant
from analysis.filter_engine import (
    apply_filters,
    get_filter_summary,
    _check_pass_filter,
    _check_af_exclusion,
    _check_consequence_excluded,
    _check_tier3_deprioritize,
    _check_clinvar_filter,
    _check_clinvar_deprioritize,
    _get_consequence_tier,
    _is_clinvar_pathogenic,
    _mark_gene_list,
    _set_gene_relevance,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    """
    Create a Variant with safe defaults for filtering tests.

    Defaults:
      filter_status="PASS", gnomad_af=None, consequence="missense_variant",
      gene_symbol="BRCA1", clinvar_significance=None, excluded=False
    """
    defaults = dict(
        chrom="chr1",
        pos=100,
        ref="A",
        alt="T",
        filter_status="PASS",
        gnomad_af=None,
        consequence="missense_variant",
        gene_symbol="BRCA1",
        clinvar_significance=None,
        excluded=False,
        exclusion_reason=None,
        in_gene_list=False,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def _make_settings(
    require_pass: bool = True,
    max_af: float = 0.01,
    ba1_af_threshold: float = 0.05,
    tier1=None,
    tier2=None,
    tier3=None,
    excluded_tier=None,
    exclude_clinvar_sigs=None,
) -> dict:
    """Build a minimal settings dict matching the config/settings.yaml structure."""
    return {
        "filtering": {
            "require_pass": require_pass,
            "max_af": max_af,
            "consequence_tiers": {
                "tier1": tier1 or [
                    "stop_gained", "frameshift_variant",
                    "splice_acceptor_variant", "splice_donor_variant",
                    "start_lost", "stop_lost",
                ],
                "tier2": tier2 or [
                    "missense_variant", "splice_region_variant",
                    "inframe_insertion", "inframe_deletion",
                ],
                "tier3": tier3 or [
                    "synonymous_variant", "5_prime_UTR_variant",
                    "3_prime_UTR_variant", "intron_variant",
                ],
                "excluded": excluded_tier or [
                    "intergenic_variant", "downstream_gene_variant",
                    "upstream_gene_variant",
                ],
            },
            "exclude_clinvar_significances": (
                exclude_clinvar_sigs
                if exclude_clinvar_sigs is not None
                else ["benign", "likely benign"]
            ),
        },
        "acmg": {
            "ba1_af_threshold": ba1_af_threshold,
        },
    }


SETTINGS = _make_settings()


# ===========================================================================
# PASS-only filter
# ===========================================================================

class TestPassFilter:
    def test_pass_retained(self):
        v = _make_variant(filter_status="PASS")
        result = _check_pass_filter(v, require_pass=True)
        assert result is None

    def test_non_pass_excluded(self):
        v = _make_variant(filter_status="LowQual")
        result = _check_pass_filter(v, require_pass=True)
        assert result is not None
        assert "LowQual" in result

    def test_dot_filter_status_excluded(self):
        """filter_status='.' must be excluded when require_pass=True (v1 strict behavior)."""
        v = _make_variant(filter_status=".")
        result = _check_pass_filter(v, require_pass=True)
        assert result is not None
        assert "." in result

    def test_require_pass_false_allows_any_status(self):
        """When require_pass=False, all filter statuses pass."""
        for status in ("PASS", ".", "LowQual", "FAIL"):
            v = _make_variant(filter_status=status)
            assert _check_pass_filter(v, require_pass=False) is None

    def test_exclusion_reason_recorded(self):
        """Exclusion reason must be non-empty and informative."""
        v = _make_variant(filter_status="SoftFilter")
        reason = _check_pass_filter(v, require_pass=True)
        assert reason
        assert len(reason) > 10  # not empty/trivial

    def test_apply_filters_pass_excluded(self):
        v = _make_variant(filter_status="LowQual")
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.exclusion_reason is not None

    def test_apply_filters_dot_excluded(self):
        v = _make_variant(filter_status=".")
        apply_filters([v], SETTINGS)
        assert v.excluded is True

    def test_apply_filters_pass_retained(self):
        v = _make_variant(filter_status="PASS")
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.exclusion_reason is None

    def test_require_pass_false_in_settings(self):
        settings = _make_settings(require_pass=False)
        v = _make_variant(filter_status="SoftFilter")
        apply_filters([v], settings)
        assert v.excluded is False


# ===========================================================================
# Allele-frequency filter — new deprioritize/exclude split
# ===========================================================================

class TestAFFilter:
    def test_af_below_threshold_retained(self):
        """gnomad_af < max_af (1%) → retained."""
        v = _make_variant(gnomad_af=0.005)
        result = _check_af_exclusion(v, max_af=0.01)
        assert result is None

    def test_af_equal_threshold_retained(self):
        """gnomad_af == max_af → retained (threshold is exclusive upper bound)."""
        v = _make_variant(gnomad_af=0.01)
        result = _check_af_exclusion(v, max_af=0.01)
        assert result is None

    def test_af_above_max_af_excluded(self):
        """gnomad_af > max_af (1%) → hard excluded per SOP."""
        v = _make_variant(gnomad_af=0.03)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.deprioritized is False

    def test_af_well_above_max_af_excluded(self):
        """gnomad_af > 5% → also excluded (same code path)."""
        v = _make_variant(gnomad_af=0.06)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.deprioritized is False

    def test_af_none_retained(self):
        """Variant not in population databases (gnomad_af=None) → retained."""
        v = _make_variant(gnomad_af=None)
        result = _check_af_exclusion(v, max_af=0.01)
        assert result is None

    def test_check_af_exclusion_above_max_af(self):
        v = _make_variant(gnomad_af=0.03)
        result = _check_af_exclusion(v, max_af=0.01)
        assert result is not None
        assert "0.03" in result or "3" in result

    def test_check_af_exclusion_below_max_af(self):
        v = _make_variant(gnomad_af=0.005)
        result = _check_af_exclusion(v, max_af=0.01)
        assert result is None

    def test_apply_filters_af_none_retained(self):
        v = _make_variant(gnomad_af=None)
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is False


class TestAFFieldSource:
    """
    The parser sets variant.gnomad_af = MAX_AF (primary) or gnomAD_AF (fallback).
    The filter engine reads variant.gnomad_af directly — it does not look at
    primary_transcript fields.  These tests verify that the contract holds:
    filter_engine excludes/deprioritizes based on whatever value gnomad_af holds.

    This mirrors the v1 design decision: MAX_AF is primary, gnomAD_AF is fallback.
    That logic lives in parser.py; filter_engine.py trusts gnomad_af as-is.
    """

    def test_high_gnomad_af_excluded(self):
        """Simulate gnomAD_AF fallback having a value above BA1 threshold → excluded."""
        v = _make_variant(gnomad_af=0.12)
        apply_filters([v], SETTINGS)
        assert v.excluded is True

    def test_high_max_af_excluded(self):
        """Simulate MAX_AF (primary) having a value above BA1 threshold → excluded."""
        v = _make_variant(gnomad_af=0.20)
        apply_filters([v], SETTINGS)
        assert v.excluded is True

    def test_intermediate_af_excluded(self):
        """AF > max_af (1%) → excluded per SOP, not deprioritized."""
        v = _make_variant(gnomad_af=0.03)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.deprioritized is False

    def test_low_max_af_retained(self):
        """MAX_AF below threshold → retained (not deprioritized)."""
        v = _make_variant(gnomad_af=0.001)
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is False

    def test_af_none_means_not_in_database(self):
        """
        gnomad_af=None means neither MAX_AF nor gnomAD_AF was present.
        Variant is treated as novel / not in population databases → retained.
        """
        v = _make_variant(gnomad_af=None)
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is False


# ===========================================================================
# Consequence-tier filtering
# ===========================================================================

class TestConsequenceTier:
    def test_tier1_always_retained(self):
        for csq in ["stop_gained", "frameshift_variant", "splice_acceptor_variant"]:
            v = _make_variant(consequence=csq)
            result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
            assert result is None, f"tier1 consequence {csq!r} should be retained"

    def test_tier2_retained(self):
        for csq in ["missense_variant", "splice_region_variant"]:
            v = _make_variant(consequence=csq)
            result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
            assert result is None, f"tier2 consequence {csq!r} should be retained"

    def test_excluded_tier_excluded(self):
        for csq in ["intergenic_variant", "downstream_gene_variant", "upstream_gene_variant"]:
            v = _make_variant(consequence=csq)
            result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
            assert result is not None, f"excluded consequence {csq!r} should be excluded"
            assert csq in result

    def test_tier3_deprioritized_not_excluded(self):
        """Tier3 consequence → deprioritized, not excluded."""
        v = _make_variant(consequence="intron_variant", in_gene_list=False)
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is True

    def test_tier3_deprioritize_function(self):
        v = _make_variant(consequence="synonymous_variant")
        result = _check_tier3_deprioritize(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is not None
        assert "tier3" in result

    def test_tier3_in_gene_list_also_deprioritized(self):
        """Tier3 in gene list → still deprioritized (gene list is just a tag now)."""
        v = _make_variant(consequence="synonymous_variant", in_gene_list=True)
        result = _check_tier3_deprioritize(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is not None

    def test_tier3_clinvar_pathogenic_deprioritized(self):
        """Tier3 with ClinVar P/LP → deprioritized (not excluded, but still tier3)."""
        v = _make_variant(
            consequence="intron_variant",
            in_gene_list=False,
            clinvar_significance="Pathogenic",
        )
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is True

    def test_unrecognised_consequence_retained(self):
        """Unknown consequences are not in any tier → pass through (conservative)."""
        v = _make_variant(consequence="some_novel_consequence")
        result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is None

    def test_empty_consequence_retained(self):
        v = _make_variant(consequence="")
        result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is None

    def test_compound_consequence_most_favorable_tier_wins(self):
        """
        'missense_variant&intron_variant' → tier2 wins over tier3.
        Variant should not be excluded or deprioritized by consequence filter.
        """
        v = _make_variant(
            consequence="missense_variant&intron_variant",
            in_gene_list=False,
        )
        result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is None
        result2 = _check_tier3_deprioritize(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result2 is None

    def test_compound_consequence_excluded_plus_tier2(self):
        """
        'missense_variant&intergenic_variant' → tier2 wins over excluded.
        Variant not excluded by consequence filter.
        """
        v = _make_variant(consequence="missense_variant&intergenic_variant")
        result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is None

    def test_compound_all_excluded(self):
        """All parts of compound consequence in excluded tier → excluded."""
        v = _make_variant(consequence="intergenic_variant&upstream_gene_variant")
        result = _check_consequence_excluded(v, SETTINGS["filtering"]["consequence_tiers"])
        assert result is not None


class TestGetConsequenceTier:
    def test_known_tier1(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("stop_gained", tier_cfg) == "tier1"

    def test_known_tier2(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("missense_variant", tier_cfg) == "tier2"

    def test_known_tier3(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("intron_variant", tier_cfg) == "tier3"

    def test_known_excluded(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("intergenic_variant", tier_cfg) == "excluded"

    def test_unknown_returns_none(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("novel_variant", tier_cfg) is None

    def test_empty_returns_none(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        assert _get_consequence_tier("", tier_cfg) is None

    def test_compound_most_favorable(self):
        tier_cfg = SETTINGS["filtering"]["consequence_tiers"]
        # tier2 wins over tier3
        assert _get_consequence_tier("missense_variant&intron_variant", tier_cfg) == "tier2"
        # tier1 wins over excluded
        assert _get_consequence_tier("stop_gained&intergenic_variant", tier_cfg) == "tier1"


# ===========================================================================
# Gene-list matching and gene relevance
# ===========================================================================

class TestGeneListMatching:
    def test_exact_match_case_insensitive(self):
        v = _make_variant(gene_symbol="BRCA1")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is True

    def test_lowercase_gene_symbol_matches(self):
        v = _make_variant(gene_symbol="brca1")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is True

    def test_mixed_case_symbol_matches(self):
        v = _make_variant(gene_symbol="Brca1")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is True

    def test_partial_match_not_matched(self):
        """Partial gene name must NOT match — only exact match is allowed."""
        v = _make_variant(gene_symbol="BRCA")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is False

    def test_prefix_not_matched(self):
        v = _make_variant(gene_symbol="BRCA1X")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is False

    def test_no_gene_list_not_marked(self):
        v = _make_variant(gene_symbol="BRCA1")
        _mark_gene_list(v, set())
        assert v.in_gene_list is False

    def test_empty_gene_symbol_not_marked(self):
        v = _make_variant(gene_symbol="")
        _mark_gene_list(v, {"BRCA1"})
        assert v.in_gene_list is False

    def test_apply_filters_gene_list_set_before_filtering(self):
        """in_gene_list must be set even on excluded variants."""
        v = _make_variant(
            gene_symbol="BRCA1",
            filter_status="LowQual",  # will be excluded by PASS filter
        )
        apply_filters([v], SETTINGS, gene_list=["BRCA1"])
        assert v.in_gene_list is True
        assert v.excluded is True  # still excluded by PASS filter

    def test_gene_list_saves_tier3_variant_as_deprioritized(self):
        """Tier3 gene in gene list → deprioritized (not excluded)."""
        v = _make_variant(consequence="intron_variant", gene_symbol="BRCA1")
        apply_filters([v], SETTINGS, gene_list=["BRCA1"])
        assert v.in_gene_list is True
        assert v.excluded is False
        assert v.deprioritized is True

    def test_gene_list_none_does_not_mark(self):
        v = _make_variant(gene_symbol="BRCA1")
        apply_filters([v], SETTINGS, gene_list=None)
        assert v.in_gene_list is False


class TestSetGeneRelevance:
    def test_no_gene_set_empty_relevance(self):
        v = _make_variant()
        _set_gene_relevance(v, set())
        assert v.gene_relevance == ""

    def test_in_gene_list_relevance(self):
        v = _make_variant(gene_symbol="BRCA1")
        _mark_gene_list(v, {"BRCA1"})
        _set_gene_relevance(v, {"BRCA1"})
        assert v.gene_relevance == "In supplied gene list"

    def test_outside_gene_list_relevance(self):
        v = _make_variant(gene_symbol="TP53")
        _mark_gene_list(v, {"BRCA1"})
        _set_gene_relevance(v, {"BRCA1"})
        assert v.gene_relevance == "Outside supplied gene list"


# ===========================================================================
# ClinVar-based deprioritization (was exclusion)
# ===========================================================================

class TestClinVarFilter:
    def test_benign_deprioritized_not_excluded(self):
        """ClinVar Benign → deprioritized, not excluded."""
        v = _make_variant(clinvar_significance="Benign")
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is True

    def test_likely_benign_deprioritized(self):
        v = _make_variant(clinvar_significance="Likely_benign")
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is not None

    def test_likely_benign_underscore_deprioritized(self):
        """ClinVar uses both 'Likely_benign' and 'Likely benign' — both must match."""
        v = _make_variant(clinvar_significance="Likely_benign")
        result = _check_clinvar_deprioritize(v, {"likely benign"})
        assert result is not None

    def test_pathogenic_retained(self):
        v = _make_variant(clinvar_significance="Pathogenic")
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is None

    def test_likely_pathogenic_retained(self):
        v = _make_variant(clinvar_significance="Likely_pathogenic")
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is None

    def test_vus_retained(self):
        v = _make_variant(clinvar_significance="Uncertain_significance")
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is None

    def test_none_significance_retained(self):
        v = _make_variant(clinvar_significance=None)
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is None

    def test_empty_significance_retained(self):
        v = _make_variant(clinvar_significance="")
        result = _check_clinvar_deprioritize(v, {"benign", "likely benign"})
        assert result is None

    def test_deprioritization_reason_contains_significance(self):
        v = _make_variant(clinvar_significance="Benign")
        reason = _check_clinvar_deprioritize(v, {"benign"})
        assert reason is not None
        assert "Benign" in reason

    def test_empty_exclude_sigs_retains_all(self):
        """If no ClinVar sigs are configured, nothing is deprioritized."""
        v = _make_variant(clinvar_significance="Benign")
        result = _check_clinvar_deprioritize(v, set())
        assert result is None

    def test_apply_filters_benign_deprioritized(self):
        v = _make_variant(clinvar_significance="Benign")
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is True
        assert v.deprioritization_reason is not None

    def test_apply_filters_pathogenic_retained(self):
        v = _make_variant(clinvar_significance="Pathogenic")
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.deprioritized is False

    def test_check_clinvar_filter_still_works(self):
        """_check_clinvar_filter is retained for backward compatibility."""
        v = _make_variant(clinvar_significance="Benign")
        result = _check_clinvar_filter(v, {"benign", "likely benign"})
        assert result is not None
        assert "Benign" in result


class TestIsClinVarPathogenic:
    def test_pathogenic(self):
        assert _is_clinvar_pathogenic("Pathogenic") is True

    def test_likely_pathogenic(self):
        assert _is_clinvar_pathogenic("Likely_pathogenic") is True

    def test_pathogenic_likely_pathogenic(self):
        assert _is_clinvar_pathogenic("Pathogenic/Likely_pathogenic") is True

    def test_benign_not_pathogenic(self):
        assert _is_clinvar_pathogenic("Benign") is False

    def test_likely_benign_not_pathogenic(self):
        assert _is_clinvar_pathogenic("Likely_benign") is False

    def test_vus_not_pathogenic(self):
        assert _is_clinvar_pathogenic("Uncertain_significance") is False

    def test_none_not_pathogenic(self):
        assert _is_clinvar_pathogenic(None) is False

    def test_conflicting_not_pathogenic(self):
        assert _is_clinvar_pathogenic("Conflicting_interpretations_of_pathogenicity") is False


# ===========================================================================
# Exclusion reason tracking
# ===========================================================================

class TestExclusionReasonTracking:
    def test_pass_filter_reason_set(self):
        v = _make_variant(filter_status="LowQ")
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.exclusion_reason is not None
        assert len(v.exclusion_reason) > 5

    def test_af_exclusion_reason_set(self):
        """AF above BA1 threshold → excluded with reason."""
        v = _make_variant(gnomad_af=0.99)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.exclusion_reason is not None
        assert "AF" in v.exclusion_reason or "0.99" in v.exclusion_reason

    def test_consequence_filter_reason_set(self):
        v = _make_variant(consequence="intergenic_variant")
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.exclusion_reason is not None
        assert "intergenic_variant" in v.exclusion_reason

    def test_clinvar_deprioritization_reason_set(self):
        v = _make_variant(clinvar_significance="Benign")
        apply_filters([v], SETTINGS)
        assert v.deprioritized is True
        assert "Benign" in v.deprioritization_reason

    def test_first_failing_filter_reason_wins(self):
        """When PASS and AF both fail, the PASS filter reason is recorded (first)."""
        v = _make_variant(filter_status="LowQ", gnomad_af=0.99)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        # PASS filter runs before AF filter
        assert "LowQ" in v.exclusion_reason or "FILTER" in v.exclusion_reason

    def test_retained_variant_has_no_reason(self):
        v = _make_variant()
        apply_filters([v], SETTINGS)
        assert v.excluded is False
        assert v.exclusion_reason is None

    def test_exclusion_reason_is_human_readable(self):
        """Reason should be a full sentence, not a code."""
        v = _make_variant(filter_status=".")
        apply_filters([v], SETTINGS)
        assert v.exclusion_reason is not None
        # Must not be just a single word or character
        assert len(v.exclusion_reason) > 20
        # Must contain some meaningful text
        assert any(c.isalpha() for c in v.exclusion_reason)


# ===========================================================================
# apply_filters: integration scenarios
# ===========================================================================

class TestApplyFilters:
    def test_returns_same_list(self):
        variants = [_make_variant(), _make_variant(gnomad_af=0.5)]
        result = apply_filters(variants, SETTINGS)
        assert result is variants  # same object, mutated in place

    def test_mixed_variants(self):
        """Test that the new filter logic correctly categorizes variants."""
        variants = [
            _make_variant(gnomad_af=0.001),                         # retained
            _make_variant(gnomad_af=0.5),                           # excluded AF (> BA1)
            _make_variant(filter_status="."),                        # excluded PASS
            _make_variant(consequence="intergenic_variant"),         # excluded tier
            _make_variant(clinvar_significance="Likely_benign"),     # deprioritized ClinVar
        ]
        apply_filters(variants, SETTINGS)
        retained = [v for v in variants if not v.excluded and not v.deprioritized]
        deprioritized = [v for v in variants if v.deprioritized]
        excluded = [v for v in variants if v.excluded]
        assert len(retained) == 1
        assert len(deprioritized) == 1
        assert len(excluded) == 3

    def test_empty_variant_list(self):
        result = apply_filters([], SETTINGS)
        assert result == []

    def test_all_retained(self):
        variants = [
            _make_variant(gnomad_af=0.001),
            _make_variant(gnomad_af=None),
        ]
        apply_filters(variants, SETTINGS)
        assert all(not v.excluded and not v.deprioritized for v in variants)

    def test_gene_list_passed_through(self):
        variants = [
            _make_variant(gene_symbol="BRCA1"),
            _make_variant(gene_symbol="TP53"),
        ]
        apply_filters(variants, SETTINGS, gene_list=["BRCA1"])
        brca1 = next(v for v in variants if v.gene_symbol == "BRCA1")
        tp53  = next(v for v in variants if v.gene_symbol == "TP53")
        assert brca1.in_gene_list is True
        assert tp53.in_gene_list is False

    def test_idempotent_on_already_excluded(self):
        """Calling apply_filters twice does not double-set fields or change reason."""
        v = _make_variant(filter_status="LowQ")
        apply_filters([v], SETTINGS)
        first_reason = v.exclusion_reason
        # Reset and re-apply
        v.excluded = False
        v.exclusion_reason = None
        apply_filters([v], SETTINGS)
        assert v.exclusion_reason == first_reason

    def test_af_above_max_af_excluded(self):
        """AF > max_af (0.01) → excluded per SOP."""
        v = _make_variant(gnomad_af=0.03)
        apply_filters([v], SETTINGS)
        assert v.excluded is True
        assert v.deprioritized is False
        assert v.exclusion_reason is not None


# ===========================================================================
# get_filter_summary
# ===========================================================================

class TestGetFilterSummary:
    def test_summary_counts(self):
        variants = [
            _make_variant(),
            _make_variant(gnomad_af=0.5),   # excluded (> BA1)
            _make_variant(filter_status="."),  # excluded
        ]
        apply_filters(variants, SETTINGS)
        summary = get_filter_summary(variants)
        assert summary["total"] == 3
        assert summary["retained"] == 1
        assert summary["excluded"] == 2

    def test_summary_af_above_threshold_excluded(self):
        """AF > 1% is now hard-excluded per SOP, not deprioritized."""
        variants = [
            _make_variant(),
            _make_variant(gnomad_af=0.03),   # excluded (> max_af per SOP)
            _make_variant(filter_status="."),  # excluded
        ]
        apply_filters(variants, SETTINGS)
        summary = get_filter_summary(variants)
        assert summary["total"] == 3
        assert summary["retained"] == 1
        assert summary["deprioritized"] == 0
        assert summary["excluded"] == 2

    def test_summary_all_retained(self):
        variants = [_make_variant(), _make_variant()]
        apply_filters(variants, SETTINGS)
        summary = get_filter_summary(variants)
        assert summary["retained"] == 2
        assert summary["excluded"] == 0
        assert summary["deprioritized"] == 0

    def test_summary_exclusion_reasons_present(self):
        v = _make_variant(gnomad_af=0.5)  # > BA1 threshold → excluded
        apply_filters([v], SETTINGS)
        summary = get_filter_summary([v])
        assert len(summary["exclusion_reasons"]) > 0

    def test_summary_empty_list(self):
        summary = get_filter_summary([])
        assert summary["total"] == 0
        assert summary["retained"] == 0
        assert summary["excluded"] == 0
        assert summary["deprioritized"] == 0
