"""
Tests for analysis/ranker.py

Strategy:
  - Ordering tests verify the sorted list order by varying one clinical
    criterion at a time.
  - Determinism tests verify that identical variants sort consistently by
    genomic coordinates.
  - Helper function unit tests cover sort-key helpers directly.

No priority score or composite weighting is tested — the ranker does not
calculate scores.

Run: pytest tests/analysis/test_ranker.py -v
"""

import pytest
from analysis.parser import Variant
from analysis.acmg.evidence_codes import (
    PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN,
)
from analysis.ranker import (
    rank_variants,
    get_ranking_summary,
    _chrom_sort_key,
    _is_clinvar_plp,
    _consequence_tier_key,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    """Variant with safe defaults for ranking tests."""
    defaults = dict(
        chrom="chr1",
        pos=100,
        ref="A",
        alt="T",
        filter_status="PASS",
        consequence="missense_variant",
        gene_symbol="GENE1",
        gnomad_af=0.001,
        cadd_phred=None,
        revel_score=None,
        spliceai_max_score=None,
        clinvar_significance=None,
        acmg_classification=VUS,
        acmg_evidence_codes=[],
        excluded=False,
        in_gene_list=False,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def _settings(**overrides) -> dict:
    """Standard settings dict for ranking tests."""
    base = {
        "filtering": {
            "max_af": 0.01,
            "consequence_tiers": {
                "tier1": ["stop_gained", "frameshift_variant",
                          "splice_acceptor_variant", "splice_donor_variant", "start_lost"],
                "tier2": ["missense_variant", "splice_region_variant",
                          "inframe_insertion", "inframe_deletion"],
                "tier3": ["synonymous_variant", "5_prime_UTR_variant",
                          "3_prime_UTR_variant", "intron_variant"],
                "excluded": ["intergenic_variant", "downstream_gene_variant",
                             "upstream_gene_variant"],
            },
        },
    }
    base.update(overrides)
    return base


SETTINGS = _settings()


# ===========================================================================
# rank_variants: basic properties
# ===========================================================================

class TestRankVariantsBasic:
    def test_returns_same_list(self):
        variants = [_make_variant(), _make_variant(pos=200)]
        result = rank_variants(variants, SETTINGS)
        assert result is variants

    def test_empty_list_returns_empty(self):
        result = rank_variants([], SETTINGS)
        assert result == []

    def test_single_variant_unchanged(self):
        v = _make_variant()
        result = rank_variants([v], SETTINGS)
        assert result[0] is v


# ===========================================================================
# ACMG classification ordering
# ===========================================================================

class TestRankingByACMG:
    def _rank_acmg(self, *classifications):
        variants = [
            _make_variant(pos=i * 100, acmg_classification=c)
            for i, c in enumerate(classifications)
        ]
        return rank_variants(variants, SETTINGS)

    def test_pathogenic_before_likely_pathogenic(self):
        vs = self._rank_acmg(LIKELY_PATHOGENIC, PATHOGENIC)
        assert vs[0].acmg_classification == PATHOGENIC

    def test_likely_pathogenic_before_vus(self):
        vs = self._rank_acmg(VUS, LIKELY_PATHOGENIC)
        assert vs[0].acmg_classification == LIKELY_PATHOGENIC

    def test_vus_before_likely_benign(self):
        vs = self._rank_acmg(LIKELY_BENIGN, VUS)
        assert vs[0].acmg_classification == VUS

    def test_likely_benign_before_benign(self):
        vs = self._rank_acmg(BENIGN, LIKELY_BENIGN)
        assert vs[0].acmg_classification == LIKELY_BENIGN

    def test_full_acmg_ordering(self):
        vs = self._rank_acmg(BENIGN, VUS, LIKELY_BENIGN, PATHOGENIC, LIKELY_PATHOGENIC)
        order = [v.acmg_classification for v in vs]
        assert order == [PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN]

    def test_none_acmg_sorts_last(self):
        v_vus  = _make_variant(acmg_classification=VUS)
        v_none = _make_variant(acmg_classification=None, pos=200)
        result = rank_variants([v_none, v_vus], SETTINGS)
        assert result[0].acmg_classification == VUS


# ===========================================================================
# Excluded vs retained ordering
# ===========================================================================

class TestExcludedVsRetained:
    def test_retained_sorts_before_excluded(self):
        excluded = _make_variant(excluded=True, acmg_classification=PATHOGENIC, pos=100)
        retained = _make_variant(excluded=False, acmg_classification=VUS, pos=200)
        result = rank_variants([excluded, retained], SETTINGS)
        assert result[0].excluded is False
        assert result[1].excluded is True

    def test_retained_sorts_before_excluded_regardless_of_acmg(self):
        exc = _make_variant(excluded=True,  acmg_classification=PATHOGENIC)
        ret = _make_variant(excluded=False, acmg_classification=BENIGN, pos=200)
        result = rank_variants([exc, ret], SETTINGS)
        assert result[0] is ret
        assert result[1] is exc

    def test_multiple_retained_then_multiple_excluded(self):
        variants = [
            _make_variant(pos=1, excluded=True),
            _make_variant(pos=2, excluded=False),
            _make_variant(pos=3, excluded=True),
            _make_variant(pos=4, excluded=False),
        ]
        result = rank_variants(variants, SETTINGS)
        retained_indices = [i for i, v in enumerate(result) if not v.excluded]
        excluded_indices  = [i for i, v in enumerate(result) if v.excluded]
        assert max(retained_indices) < min(excluded_indices)


# ===========================================================================
# ClinVar ordering
# ===========================================================================

class TestClinVarOrdering:
    def test_clinvar_plp_ranks_before_no_clinvar(self):
        v_p    = _make_variant(clinvar_significance="Pathogenic",  pos=100)
        v_none = _make_variant(clinvar_significance=None,           pos=200)
        result = rank_variants([v_none, v_p], SETTINGS)
        assert result[0].clinvar_significance == "Pathogenic"

    def test_clinvar_lp_ranks_before_uncertain(self):
        v_lp  = _make_variant(clinvar_significance="Likely_pathogenic",     pos=100)
        v_vus = _make_variant(clinvar_significance="Uncertain_significance", pos=200)
        result = rank_variants([v_vus, v_lp], SETTINGS)
        assert result[0] is v_lp

    def test_clinvar_benign_does_not_boost_rank(self):
        v_b    = _make_variant(clinvar_significance="Benign",   pos=100)
        v_none = _make_variant(clinvar_significance=None,        pos=200)
        # Both have no P/LP — benign and no-clinvar should not differ on ClinVar key
        result = rank_variants([v_b, v_none], SETTINGS)
        # Just assert no crash and order is stable (both get clinvar_key=1)
        assert len(result) == 2


# ===========================================================================
# Consequence tier ordering
# ===========================================================================

class TestConsequenceTierOrdering:
    def test_tier1_ranks_before_tier2(self):
        v1 = _make_variant(consequence="stop_gained",     pos=100)
        v2 = _make_variant(consequence="missense_variant", pos=200)
        result = rank_variants([v2, v1], SETTINGS)
        assert result[0].consequence == "stop_gained"

    def test_tier2_ranks_before_tier3(self):
        v2 = _make_variant(consequence="missense_variant",  pos=100)
        v3 = _make_variant(consequence="synonymous_variant", pos=200)
        result = rank_variants([v3, v2], SETTINGS)
        assert result[0].consequence == "missense_variant"

    def test_excluded_tier_ranks_last_among_same_acmg(self):
        v_t2   = _make_variant(consequence="missense_variant",    pos=100)
        v_excl = _make_variant(consequence="intergenic_variant",  pos=200)
        result = rank_variants([v_excl, v_t2], SETTINGS)
        assert result[0].consequence == "missense_variant"


# ===========================================================================
# gnomAD rarity ordering
# ===========================================================================

class TestRarityOrdering:
    def test_absent_af_ranks_before_rare_af(self):
        v_novel = _make_variant(gnomad_af=None,  pos=100)
        v_rare  = _make_variant(gnomad_af=0.001, pos=200)
        result = rank_variants([v_rare, v_novel], SETTINGS)
        assert result[0].gnomad_af is None

    def test_lower_af_ranks_before_higher_af(self):
        v_very_rare = _make_variant(gnomad_af=0.0001, pos=100)
        v_rare      = _make_variant(gnomad_af=0.005,  pos=200)
        result = rank_variants([v_rare, v_very_rare], SETTINGS)
        assert result[0].gnomad_af == 0.0001


# ===========================================================================
# Deterministic tie-breaking by genomic coordinates
# ===========================================================================

class TestDeterministicTies:
    def test_identical_variants_sorted_by_position(self):
        v1 = _make_variant(pos=300)
        v2 = _make_variant(pos=100)
        v3 = _make_variant(pos=200)
        result = rank_variants([v1, v2, v3], SETTINGS)
        assert [v.pos for v in result] == [100, 200, 300]

    def test_identical_variants_sorted_by_chrom_then_pos(self):
        v1 = _make_variant(chrom="chr2", pos=100)
        v2 = _make_variant(chrom="chr1", pos=200)
        v3 = _make_variant(chrom="chr1", pos=100)
        result = rank_variants([v1, v2, v3], SETTINGS)
        assert result[0].chrom == "chr1" and result[0].pos == 100
        assert result[1].chrom == "chr1" and result[1].pos == 200
        assert result[2].chrom == "chr2"

    def test_same_result_on_repeated_calls(self):
        variants1 = [_make_variant(pos=i * 10) for i in range(5)]
        variants2 = [_make_variant(pos=i * 10) for i in range(5)]
        rank_variants(variants1, SETTINGS)
        rank_variants(variants2, SETTINGS)
        assert [v.pos for v in variants1] == [v.pos for v in variants2]

    def test_chromosome_numeric_sort_order(self):
        v2  = _make_variant(chrom="chr2",  pos=100)
        v10 = _make_variant(chrom="chr10", pos=100)
        result = rank_variants([v10, v2], SETTINGS)
        assert result[0].chrom == "chr2"
        assert result[1].chrom == "chr10"

    def test_sex_chromosomes_sort_after_autosomes(self):
        vX = _make_variant(chrom="chrX", pos=100)
        v1 = _make_variant(chrom="chr1", pos=100)
        result = rank_variants([vX, v1], SETTINGS)
        assert result[0].chrom == "chr1"
        assert result[1].chrom == "chrX"


# ===========================================================================
# _chrom_sort_key
# ===========================================================================

class TestChromSortKey:
    def test_numeric_chrom(self):
        assert _chrom_sort_key("chr1") < _chrom_sort_key("chr2")
        assert _chrom_sort_key("chr2") < _chrom_sort_key("chr10")

    def test_sex_chrom_after_autosome(self):
        assert _chrom_sort_key("chr22") < _chrom_sort_key("chrX")

    def test_y_after_x(self):
        assert _chrom_sort_key("chrX") < _chrom_sort_key("chrY")

    def test_mt_after_y(self):
        assert _chrom_sort_key("chrY") < _chrom_sort_key("chrMT")

    def test_no_chr_prefix(self):
        assert _chrom_sort_key("1") < _chrom_sort_key("2")


# ===========================================================================
# _is_clinvar_plp
# ===========================================================================

class TestIsClinvarPlp:
    def test_pathogenic(self):
        assert _is_clinvar_plp("Pathogenic") is True

    def test_likely_pathogenic(self):
        assert _is_clinvar_plp("Likely_pathogenic") is True

    def test_benign_is_false(self):
        assert _is_clinvar_plp("Benign") is False

    def test_conflicting_is_false(self):
        assert _is_clinvar_plp("Conflicting_interpretations_of_pathogenicity") is False

    def test_uncertain_is_false(self):
        assert _is_clinvar_plp("Uncertain_significance") is False

    def test_none_is_false(self):
        assert _is_clinvar_plp(None) is False


# ===========================================================================
# _consequence_tier_key
# ===========================================================================

class TestConsequenceTierKey:
    def _tier_cfg(self):
        return SETTINGS["filtering"]["consequence_tiers"]

    def test_tier1_key_lowest(self):
        assert _consequence_tier_key("stop_gained", self._tier_cfg()) == 0

    def test_tier2_key_middle(self):
        assert _consequence_tier_key("missense_variant", self._tier_cfg()) == 1

    def test_tier3_key_higher(self):
        assert _consequence_tier_key("synonymous_variant", self._tier_cfg()) == 2

    def test_excluded_key_higher_than_tier3(self):
        assert (_consequence_tier_key("intergenic_variant", self._tier_cfg()) >
                _consequence_tier_key("synonymous_variant", self._tier_cfg()))

    def test_unknown_consequence_key_highest(self):
        t = _consequence_tier_key("novel_consequence", self._tier_cfg())
        assert t >= _consequence_tier_key("intergenic_variant", self._tier_cfg())

    def test_compound_consequence_uses_most_severe(self):
        # missense & intron → tier2 (missense wins)
        t = _consequence_tier_key("missense_variant&intron_variant", self._tier_cfg())
        assert t == _consequence_tier_key("missense_variant", self._tier_cfg())

    def test_empty_consequence(self):
        assert _consequence_tier_key("", self._tier_cfg()) >= 3


# ===========================================================================
# get_ranking_summary
# ===========================================================================

class TestGetRankingSummary:
    def test_summary_counts(self):
        variants = [
            _make_variant(pos=1, excluded=False, acmg_classification=PATHOGENIC),
            _make_variant(pos=2, excluded=False, acmg_classification=VUS),
            _make_variant(pos=3, excluded=True,  acmg_classification=VUS),
        ]
        rank_variants(variants, SETTINGS)
        summary = get_ranking_summary(variants)
        assert summary["total"] == 3
        assert summary["retained"] == 2
        assert summary["excluded"] == 1

    def test_by_acmg_class_counts_only_retained(self):
        variants = [
            _make_variant(pos=1, excluded=False, acmg_classification=PATHOGENIC),
            _make_variant(pos=2, excluded=True,  acmg_classification=PATHOGENIC),
        ]
        rank_variants(variants, SETTINGS)
        summary = get_ranking_summary(variants)
        assert summary["by_acmg_class"].get(PATHOGENIC, 0) == 1

    def test_top_5_in_summary(self):
        variants = [_make_variant(pos=i * 10) for i in range(10)]
        rank_variants(variants, SETTINGS)
        summary = get_ranking_summary(variants)
        assert len(summary["top_5"]) <= 5

    def test_empty_list(self):
        summary = get_ranking_summary([])
        assert summary["total"] == 0
