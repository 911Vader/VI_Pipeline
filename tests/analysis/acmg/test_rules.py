"""
Tests for analysis/acmg/rules.py

Strategy:
  - Each rule is tested independently with minimal Variant fixtures.
  - Boundary conditions and missing-data (None) cases are always covered.
  - Tests verify that rules do NOT fire when data is absent or insufficient.
  - No unsupported evidence is fabricated (all None-field tests confirm no code).

Run: pytest tests/analysis/acmg/test_rules.py -v
"""

import pytest
from analysis.parser import Variant
from analysis.acmg.rules import (
    eval_PVS1,
    eval_PS3,
    eval_PM2,
    eval_PM4,
    eval_PP3,
    eval_BA1,
    eval_BS1,
    eval_BP4,
    eval_BP7,
    _is_clinvar_p_or_lp,
    _clinvar_review_stars,
    IMPLEMENTED_RULES,
    UNIMPLEMENTED_V1_RULES,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    """Minimal Variant with safe defaults (all None for optional score fields)."""
    defaults = dict(
        chrom="chr1", pos=100, ref="A", alt="T",
        filter_status="PASS",
        consequence="",
        gene_symbol="BRCA1",
        gnomad_af=None,
        cadd_phred=None,
        revel_score=None,
        spliceai_max_score=None,
        alphamissense_score=None,
        clinvar_significance=None,
        clinvar_review_status=None,
        loftee_lof=None,
        excluded=False,
        in_gene_list=False,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def _settings(**overrides) -> dict:
    """Minimal settings dict with standard ACMG thresholds."""
    base = {
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
    base["acmg"].update(overrides)
    return base


SETTINGS = _settings()


# ===========================================================================
# Catalogue helpers
# ===========================================================================

class TestCatalogueHelpers:
    def test_implemented_rules_keys_match_catalogue(self):
        """Every key in IMPLEMENTED_RULES must be a known ACMG code."""
        from analysis.acmg.evidence_codes import EVIDENCE_CATALOGUE
        for code in IMPLEMENTED_RULES:
            assert code in EVIDENCE_CATALOGUE, f"{code} not in evidence catalogue"

    def test_unimplemented_rules_have_reasons(self):
        for code, reason in UNIMPLEMENTED_V1_RULES.items():
            assert reason, f"{code} unimplemented reason is empty"

    def test_implemented_and_unimplemented_do_not_overlap(self):
        overlap = set(IMPLEMENTED_RULES) & set(UNIMPLEMENTED_V1_RULES)
        assert not overlap, f"Codes in both dicts: {overlap}"


# ===========================================================================
# PVS1 — null variant types
# ===========================================================================

class TestEvalPVS1:
    def test_stop_gained_triggers(self):
        v = _make_variant(consequence="stop_gained", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_frameshift_triggers(self):
        v = _make_variant(consequence="frameshift_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_splice_acceptor_triggers(self):
        v = _make_variant(consequence="splice_acceptor_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_splice_donor_triggers(self):
        v = _make_variant(consequence="splice_donor_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_start_lost_triggers(self):
        v = _make_variant(consequence="start_lost", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_missense_does_not_trigger(self):
        v = _make_variant(consequence="missense_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_synonymous_does_not_trigger(self):
        v = _make_variant(consequence="synonymous_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_inframe_indel_does_not_trigger(self):
        """in-frame indels → PM4, not PVS1."""
        v = _make_variant(consequence="inframe_deletion", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_stop_lost_does_not_trigger(self):
        """stop_lost → PM4 (protein elongation), not PVS1."""
        v = _make_variant(consequence="stop_lost", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_empty_consequence_does_not_trigger(self):
        v = _make_variant(consequence="", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_compound_consequence_with_null_triggers(self):
        """A compound consequence containing a null type + LOFTEE HC should trigger PVS1."""
        v = _make_variant(consequence="stop_gained&splice_region_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_compound_consequence_without_null_does_not_trigger(self):
        v = _make_variant(consequence="missense_variant&splice_region_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is False

    # New tests for LOFTEE HC requirement
    def test_stop_gained_without_loftee_hc_does_not_trigger(self):
        """consequence=stop_gained, loftee_lof=None → False (LOFTEE HC required)."""
        v = _make_variant(consequence="stop_gained", loftee_lof=None)
        assert eval_PVS1(v, SETTINGS) is False

    def test_frameshift_with_loftee_lc_does_not_trigger(self):
        """consequence=frameshift_variant, loftee_lof=LC → False."""
        v = _make_variant(consequence="frameshift_variant", loftee_lof="LC")
        assert eval_PVS1(v, SETTINGS) is False

    def test_frameshift_with_loftee_hc_triggers(self):
        """consequence=frameshift_variant, loftee_lof=HC → True."""
        v = _make_variant(consequence="frameshift_variant", loftee_lof="HC")
        assert eval_PVS1(v, SETTINGS) is True

    def test_compound_with_null_and_loftee_hc_triggers(self):
        """Compound consequence with null type and LOFTEE HC → True."""
        v = _make_variant(
            consequence="frameshift_variant&splice_region_variant",
            loftee_lof="HC",
        )
        assert eval_PVS1(v, SETTINGS) is True


# ===========================================================================
# PS3 — ClinVar P/LP proxy
# ===========================================================================

class TestEvalPS3:
    def test_pathogenic_two_star_triggers(self):
        v = _make_variant(
            clinvar_significance="Pathogenic",
            clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts",
        )
        assert eval_PS3(v, SETTINGS) is True

    def test_likely_pathogenic_expert_panel_triggers(self):
        v = _make_variant(
            clinvar_significance="Likely_pathogenic",
            clinvar_review_status="reviewed_by_expert_panel",
        )
        assert eval_PS3(v, SETTINGS) is True

    def test_pathogenic_practice_guideline_triggers(self):
        v = _make_variant(
            clinvar_significance="Pathogenic",
            clinvar_review_status="practice_guideline",
        )
        assert eval_PS3(v, SETTINGS) is True

    def test_pathogenic_single_submitter_does_not_trigger(self):
        """1-star review status is insufficient for PS3 proxy."""
        v = _make_variant(
            clinvar_significance="Pathogenic",
            clinvar_review_status="criteria_provided,_single_submitter",
        )
        assert eval_PS3(v, SETTINGS) is False

    def test_benign_does_not_trigger(self):
        v = _make_variant(
            clinvar_significance="Benign",
            clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts",
        )
        assert eval_PS3(v, SETTINGS) is False

    def test_vus_does_not_trigger(self):
        v = _make_variant(
            clinvar_significance="Uncertain_significance",
            clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts",
        )
        assert eval_PS3(v, SETTINGS) is False

    def test_conflicting_does_not_trigger(self):
        v = _make_variant(
            clinvar_significance="Conflicting_interpretations_of_pathogenicity",
            clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts",
        )
        assert eval_PS3(v, SETTINGS) is False

    def test_none_significance_does_not_trigger(self):
        v = _make_variant(clinvar_significance=None, clinvar_review_status=None)
        assert eval_PS3(v, SETTINGS) is False

    def test_pathogenic_no_review_status_does_not_trigger(self):
        """P/LP with missing review status → 0 stars → does not trigger."""
        v = _make_variant(
            clinvar_significance="Pathogenic",
            clinvar_review_status=None,
        )
        assert eval_PS3(v, SETTINGS) is False


class TestClinVarHelpers:
    def test_is_p_or_lp_pathogenic(self):
        assert _is_clinvar_p_or_lp("Pathogenic") is True

    def test_is_p_or_lp_likely_pathogenic(self):
        assert _is_clinvar_p_or_lp("Likely_pathogenic") is True

    def test_is_p_or_lp_pathogenic_slash_lp(self):
        assert _is_clinvar_p_or_lp("Pathogenic/Likely_pathogenic") is True

    def test_is_p_or_lp_benign(self):
        assert _is_clinvar_p_or_lp("Benign") is False

    def test_is_p_or_lp_conflicting(self):
        assert _is_clinvar_p_or_lp("Conflicting_interpretations_of_pathogenicity") is False

    def test_is_p_or_lp_none(self):
        assert _is_clinvar_p_or_lp(None) is False

    def test_review_stars_practice_guideline(self):
        assert _clinvar_review_stars("practice_guideline") == 4

    def test_review_stars_expert_panel(self):
        assert _clinvar_review_stars("reviewed_by_expert_panel") == 3

    def test_review_stars_multiple_submitters(self):
        assert _clinvar_review_stars(
            "criteria_provided,_multiple_submitters,_no_conflicts"
        ) == 2

    def test_review_stars_single_submitter(self):
        assert _clinvar_review_stars("criteria_provided,_single_submitter") == 1

    def test_review_stars_none(self):
        assert _clinvar_review_stars(None) == 0

    def test_review_stars_no_assertion(self):
        assert _clinvar_review_stars("no_assertion_criteria_provided") == 0


# ===========================================================================
# PM2 — absent / extremely rare
# ===========================================================================

class TestEvalPM2:
    def test_none_af_triggers(self):
        """Absent from population databases → PM2."""
        v = _make_variant(gnomad_af=None)
        assert eval_PM2(v, SETTINGS) is True

    def test_very_rare_af_triggers(self):
        v = _make_variant(gnomad_af=0.00005)   # < 0.0001
        assert eval_PM2(v, SETTINGS) is True

    def test_at_threshold_does_not_trigger(self):
        """gnomad_af == pm2 threshold: not below → does not trigger."""
        v = _make_variant(gnomad_af=0.0001)
        assert eval_PM2(v, SETTINGS) is False

    def test_common_variant_does_not_trigger(self):
        v = _make_variant(gnomad_af=0.001)
        assert eval_PM2(v, SETTINGS) is False

    def test_custom_threshold(self):
        """Custom pm2 threshold is respected."""
        settings = _settings(pm2_af_threshold=0.001)
        v = _make_variant(gnomad_af=0.0005)
        assert eval_PM2(v, settings) is True
        v2 = _make_variant(gnomad_af=0.002)
        assert eval_PM2(v2, settings) is False


# ===========================================================================
# PM4 — protein length change
# ===========================================================================

class TestEvalPM4:
    def test_inframe_insertion_triggers(self):
        v = _make_variant(consequence="inframe_insertion")
        assert eval_PM4(v, SETTINGS) is True

    def test_inframe_deletion_triggers(self):
        v = _make_variant(consequence="inframe_deletion")
        assert eval_PM4(v, SETTINGS) is True

    def test_stop_lost_triggers(self):
        v = _make_variant(consequence="stop_lost")
        assert eval_PM4(v, SETTINGS) is True

    def test_missense_does_not_trigger(self):
        v = _make_variant(consequence="missense_variant")
        assert eval_PM4(v, SETTINGS) is False

    def test_frameshift_does_not_trigger(self):
        """Frameshifts → PVS1, not PM4."""
        v = _make_variant(consequence="frameshift_variant")
        assert eval_PM4(v, SETTINGS) is False

    def test_stop_gained_does_not_trigger(self):
        v = _make_variant(consequence="stop_gained")
        assert eval_PM4(v, SETTINGS) is False

    def test_empty_consequence_does_not_trigger(self):
        v = _make_variant(consequence="")
        assert eval_PM4(v, SETTINGS) is False

    def test_compound_with_inframe_triggers(self):
        v = _make_variant(consequence="inframe_deletion&splice_region_variant")
        assert eval_PM4(v, SETTINGS) is True


# ===========================================================================
# PP3 — in-silico damaging (multiple tools)
# ===========================================================================

class TestEvalPP3:
    def test_group1_plus_group2_triggers(self):
        """CADD (Group 1) + SpliceAI (Group 2) = 2 independent groups → PP3 fires."""
        v = _make_variant(cadd_phred=30, spliceai_max_score=0.6)
        assert eval_PP3(v, SETTINGS) is True

    def test_cadd_revel_alone_does_not_trigger(self):
        """CADD + REVEL are correlated (Group 1) — only 1 vote — PP3 does NOT fire."""
        v = _make_variant(cadd_phred=30, revel_score=0.85)
        assert eval_PP3(v, SETTINGS) is False

    def test_three_groups_all_damaging_triggers(self):
        """All three groups vote damaging → PP3 fires (≥2 required)."""
        v = _make_variant(cadd_phred=30, revel_score=0.85, spliceai_max_score=0.6)
        assert eval_PP3(v, SETTINGS) is True

    def test_one_tool_only_does_not_trigger(self):
        """Single group (CADD only) = 1 vote — insufficient for PP3."""
        v = _make_variant(cadd_phred=30)
        assert eval_PP3(v, SETTINGS) is False

    def test_no_tools_does_not_trigger(self):
        v = _make_variant()  # all score fields None
        assert eval_PP3(v, SETTINGS) is False

    def test_intragroup_conflict_cancels_group1_vote(self):
        """CADD damaging + REVEL clearly benign (≤ revel_benign) → conflict cancels Group 1 vote."""
        v = _make_variant(cadd_phred=30, revel_score=0.10)
        assert eval_PP3(v, SETTINGS) is False

    def test_revel_gray_zone_does_not_cancel_cadd_damaging_vote(self):
        """CADD damaging + REVEL gray zone (0.15 < score < 0.75) → REVEL is neutral, not conflicting.
        Group 1 still votes damaging. Combined with SpliceAI, PP3 fires."""
        v = _make_variant(cadd_phred=30, revel_score=0.50, spliceai_max_score=0.6)
        assert eval_PP3(v, SETTINGS) is True

    def test_revel_gray_zone_alone_insufficient(self):
        """CADD damaging + REVEL gray zone, no other group → still only 1 group voting → no PP3."""
        v = _make_variant(cadd_phred=30, revel_score=0.50)
        assert eval_PP3(v, SETTINGS) is False

    def test_two_tools_both_below_threshold_does_not_trigger(self):
        v = _make_variant(cadd_phred=15, revel_score=0.50)
        assert eval_PP3(v, SETTINGS) is False

    def test_cadd_at_threshold_plus_spliceai_triggers(self):
        """CADD == threshold (≥ is inclusive) + SpliceAI = 2 independent groups."""
        v = _make_variant(cadd_phred=20.0, spliceai_max_score=0.6)
        assert eval_PP3(v, SETTINGS) is True

    def test_revel_at_threshold_plus_spliceai_triggers(self):
        """REVEL == pathogenic threshold (≥ inclusive) + SpliceAI = 2 groups."""
        v = _make_variant(revel_score=0.75, spliceai_max_score=0.6)
        assert eval_PP3(v, SETTINGS) is True

    def test_spliceai_damaging_with_cadd_triggers(self):
        v = _make_variant(cadd_phred=22, spliceai_max_score=0.55)
        assert eval_PP3(v, SETTINGS) is True

    def test_cadd_only_one_group_does_not_trigger(self):
        """CADD only (REVEL=None, SpliceAI=None) → 1 group → no PP3."""
        v = _make_variant(cadd_phred=35, revel_score=None, spliceai_max_score=None)
        assert eval_PP3(v, SETTINGS) is False


# ===========================================================================
# BA1 — AF > 5%
# ===========================================================================

class TestEvalBA1:
    def test_high_af_triggers(self):
        v = _make_variant(gnomad_af=0.10)
        assert eval_BA1(v, SETTINGS) is True

    def test_af_just_above_threshold_triggers(self):
        v = _make_variant(gnomad_af=0.051)
        assert eval_BA1(v, SETTINGS) is True

    def test_af_at_threshold_does_not_trigger(self):
        """gnomad_af == 0.05: not above → does not trigger."""
        v = _make_variant(gnomad_af=0.05)
        assert eval_BA1(v, SETTINGS) is False

    def test_rare_variant_does_not_trigger(self):
        v = _make_variant(gnomad_af=0.001)
        assert eval_BA1(v, SETTINGS) is False

    def test_none_af_does_not_trigger(self):
        v = _make_variant(gnomad_af=None)
        assert eval_BA1(v, SETTINGS) is False

    def test_custom_threshold(self):
        settings = _settings(ba1_af_threshold=0.01)
        v = _make_variant(gnomad_af=0.02)
        assert eval_BA1(v, settings) is True


# ===========================================================================
# BS1 — AF > disorder threshold
# ===========================================================================

class TestEvalBS1:
    def test_af_above_bs1_threshold_triggers(self):
        v = _make_variant(gnomad_af=0.02)   # > 0.01 default
        assert eval_BS1(v, SETTINGS) is True

    def test_af_at_bs1_threshold_does_not_trigger(self):
        v = _make_variant(gnomad_af=0.01)
        assert eval_BS1(v, SETTINGS) is False

    def test_af_below_bs1_threshold_does_not_trigger(self):
        v = _make_variant(gnomad_af=0.005)
        assert eval_BS1(v, SETTINGS) is False

    def test_none_af_does_not_trigger(self):
        v = _make_variant(gnomad_af=None)
        assert eval_BS1(v, SETTINGS) is False

    def test_ba1_level_also_triggers_bs1(self):
        """AF > BA1 threshold also exceeds BS1 threshold."""
        v = _make_variant(gnomad_af=0.10)
        assert eval_BS1(v, SETTINGS) is True


# ===========================================================================
# BP4 — in-silico benign (multiple tools)
# ===========================================================================

class TestEvalBP4:
    def test_group1_plus_group2_triggers(self):
        """REVEL benign (Group 1) + SpliceAI low (Group 2) = 2 independent groups → BP4 fires."""
        v = _make_variant(revel_score=0.05, spliceai_max_score=0.02)
        assert eval_BP4(v, SETTINGS) is True

    def test_revel_and_cadd_alone_does_not_trigger(self):
        """REVEL + CADD are correlated (Group 1) — only 1 vote — BP4 does NOT fire."""
        v = _make_variant(revel_score=0.05, cadd_phred=10)
        assert eval_BP4(v, SETTINGS) is False

    def test_all_three_groups_benign_triggers(self):
        """Group 1 (REVEL+CADD benign) + Group 2 (SpliceAI low) = 2 groups → BP4 fires."""
        v = _make_variant(revel_score=0.05, spliceai_max_score=0.02, cadd_phred=12)
        assert eval_BP4(v, SETTINGS) is True

    def test_one_group_only_does_not_trigger(self):
        """REVEL only (Group 1 only) = 1 vote — insufficient for BP4."""
        v = _make_variant(revel_score=0.05)
        assert eval_BP4(v, SETTINGS) is False

    def test_no_tools_does_not_trigger(self):
        v = _make_variant()
        assert eval_BP4(v, SETTINGS) is False

    def test_high_spliceai_cancels_group2_vote(self):
        """REVEL benign (Group 1 vote) + SpliceAI high (Group 2 NOT benign) → only 1 vote."""
        v = _make_variant(revel_score=0.05, spliceai_max_score=0.6)
        assert eval_BP4(v, SETTINGS) is False

    def test_revel_at_benign_threshold_plus_spliceai_triggers(self):
        """REVEL == 0.15 (≤ inclusive) counts as benign in Group 1; SpliceAI provides Group 2."""
        v = _make_variant(revel_score=0.15, spliceai_max_score=0.05)
        assert eval_BP4(v, SETTINGS) is True

    def test_revel_in_gray_zone_no_group1_vote(self):
        """REVEL 0.15–0.75 is the gray zone — no signal from Group 1."""
        v = _make_variant(revel_score=0.50, spliceai_max_score=0.05)
        assert eval_BP4(v, SETTINGS) is False

    def test_intragroup_conflict_cancels_group1_vote(self):
        """CADD damaging + REVEL benign inside Group 1 → conflict → group1 vote cancels."""
        v = _make_variant(cadd_phred=25, revel_score=0.05, spliceai_max_score=0.02)
        assert eval_BP4(v, SETTINGS) is False


# ===========================================================================
# BP7 — synonymous + no splice impact
# ===========================================================================

class TestEvalBP7:
    def test_synonymous_no_spliceai_triggers(self):
        """Synonymous with no SpliceAI data → no positive splice evidence → BP7."""
        v = _make_variant(consequence="synonymous_variant", spliceai_max_score=None)
        assert eval_BP7(v, SETTINGS) is True

    def test_synonymous_low_spliceai_triggers(self):
        v = _make_variant(consequence="synonymous_variant", spliceai_max_score=0.02)
        assert eval_BP7(v, SETTINGS) is True

    def test_synonymous_high_spliceai_does_not_trigger(self):
        """SpliceAI shows potential splice impact → BP7 not appropriate."""
        v = _make_variant(consequence="synonymous_variant", spliceai_max_score=0.7)
        assert eval_BP7(v, SETTINGS) is False

    def test_synonymous_spliceai_at_low_threshold_triggers(self):
        """spliceai_max_score == spliceai_low threshold (≤ inclusive) → triggers."""
        v = _make_variant(consequence="synonymous_variant", spliceai_max_score=0.1)
        assert eval_BP7(v, SETTINGS) is True

    def test_missense_does_not_trigger(self):
        v = _make_variant(consequence="missense_variant", spliceai_max_score=None)
        assert eval_BP7(v, SETTINGS) is False

    def test_frameshift_does_not_trigger(self):
        v = _make_variant(consequence="frameshift_variant", spliceai_max_score=None)
        assert eval_BP7(v, SETTINGS) is False

    def test_empty_consequence_does_not_trigger(self):
        v = _make_variant(consequence="", spliceai_max_score=None)
        assert eval_BP7(v, SETTINGS) is False

    def test_compound_containing_synonymous_triggers(self):
        """Compound consequence including synonymous → BP7 if no splice impact."""
        v = _make_variant(
            consequence="synonymous_variant&intron_variant",
            spliceai_max_score=0.05,
        )
        assert eval_BP7(v, SETTINGS) is True

    def test_compound_containing_synonymous_high_spliceai_does_not_trigger(self):
        v = _make_variant(
            consequence="synonymous_variant&intron_variant",
            spliceai_max_score=0.8,
        )
        assert eval_BP7(v, SETTINGS) is False


# ===========================================================================
# Missing-data safety — no rule should raise on all-None variant
# ===========================================================================

class TestMissingDataSafety:
    """
    All rule functions must handle a variant with all optional fields=None
    without raising an exception.  Missing data must result in no code, not an error.
    """
    def _all_none_variant(self) -> Variant:
        return _make_variant(
            consequence="",
            gnomad_af=None,
            cadd_phred=None,
            revel_score=None,
            spliceai_max_score=None,
            clinvar_significance=None,
            clinvar_review_status=None,
        )

    def test_pvs1_all_none_no_error(self):
        assert eval_PVS1(self._all_none_variant(), SETTINGS) is False

    def test_ps3_all_none_no_error(self):
        assert eval_PS3(self._all_none_variant(), SETTINGS) is False

    def test_pm2_all_none_triggers_absent_from_db(self):
        """gnomad_af=None means absent from databases → PM2 fires."""
        assert eval_PM2(self._all_none_variant(), SETTINGS) is True

    def test_pm4_all_none_no_error(self):
        assert eval_PM4(self._all_none_variant(), SETTINGS) is False

    def test_pp3_all_none_no_error(self):
        """No score data → cannot confirm multiple tools → no PP3."""
        assert eval_PP3(self._all_none_variant(), SETTINGS) is False

    def test_ba1_all_none_no_error(self):
        assert eval_BA1(self._all_none_variant(), SETTINGS) is False

    def test_bs1_all_none_no_error(self):
        assert eval_BS1(self._all_none_variant(), SETTINGS) is False

    def test_bp4_all_none_no_error(self):
        assert eval_BP4(self._all_none_variant(), SETTINGS) is False

    def test_bp7_all_none_no_error(self):
        assert eval_BP7(self._all_none_variant(), SETTINGS) is False

    def test_all_rules_safe_on_minimal_variant(self):
        """Exhaustive safety: every implemented rule must not raise."""
        v = self._all_none_variant()
        for code, rule_fn in IMPLEMENTED_RULES.items():
            try:
                rule_fn(v, SETTINGS)
            except Exception as exc:
                pytest.fail(f"Rule {code} raised on all-None variant: {exc}")
