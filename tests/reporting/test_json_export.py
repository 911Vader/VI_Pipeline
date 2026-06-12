"""Tests for reporting/json_export.py"""
import json
import os
import tempfile
from pathlib import Path

import pytest

from analysis.parser import Variant
from reporting.json_export import (
    MANDATORY_DISCLAIMER,
    PIPELINE_VERSION,
    write_json,
    _variant_to_dict,
    _build_summary,
)


def _make_variant(
    chrom="chr1", pos=100, ref="A", alt="T",
    gene_symbol="GENE1", acmg_classification="VUS",
    excluded=False, exclusion_reason=None,
    acmg_evidence_codes=None,
):
    v = Variant(
        chrom=chrom, pos=pos, ref=ref, alt=alt,
        gene_symbol=gene_symbol,
        acmg_classification=acmg_classification,
        excluded=excluded,
        exclusion_reason=exclusion_reason,
        acmg_evidence_codes=acmg_evidence_codes or [],
    )
    return v


# ---------------------------------------------------------------------------
# _variant_to_dict
# ---------------------------------------------------------------------------

class TestVariantToDict:
    def test_retained_variant_has_rank(self):
        v = _make_variant()
        d = _variant_to_dict(v, rank=3)
        assert d["rank"] == 3

    def test_excluded_variant_rank_none(self):
        v = _make_variant(excluded=True, exclusion_reason="gnomAD AF > 0.01")
        d = _variant_to_dict(v, rank=None)
        assert d["rank"] is None
        assert d["exclusion_reason"] == "gnomAD AF > 0.01"

    def test_all_expected_keys_present(self):
        v = _make_variant()
        d = _variant_to_dict(v, rank=1)
        expected_keys = {
            "rank", "gene_symbol", "chromosome", "position",
            "ref", "alt", "variant_id", "transcript_id", "hgvs_c", "hgvs_p",
            "consequence", "zygosity", "acmg_classification", "acmg_evidence_codes",
            "clinvar_significance", "clinvar_review_status", "clinvar_disease",
            "gnomad_af", "cadd_phred", "revel_score", "spliceai_max_score",
            "alphamissense_score", "in_gene_list", "excluded", "exclusion_reason",
        }
        assert expected_keys.issubset(d.keys())

    def test_gene_symbol_none_when_empty(self):
        v = _make_variant(gene_symbol="")
        d = _variant_to_dict(v, rank=1)
        assert d["gene_symbol"] is None


# ---------------------------------------------------------------------------
# _build_summary
# ---------------------------------------------------------------------------

class TestBuildSummary:
    def test_counts_retained_and_excluded(self):
        variants = [
            _make_variant(acmg_classification="Pathogenic"),
            _make_variant(acmg_classification="VUS"),
            _make_variant(excluded=True),
        ]
        s = _build_summary(variants)
        assert s["total_variants"] == 3
        assert s["retained"] == 2
        assert s["excluded"] == 1

    def test_acmg_counts(self):
        variants = [
            _make_variant(acmg_classification="Pathogenic"),
            _make_variant(acmg_classification="Pathogenic"),
            _make_variant(acmg_classification="VUS"),
        ]
        s = _build_summary(variants)
        assert s["acmg_classification_counts"]["Pathogenic"] == 2
        assert s["acmg_classification_counts"]["VUS"] == 1

    def test_top_gene_and_acmg_class(self):
        variants = [
            _make_variant(gene_symbol="BRCA1", acmg_classification="Pathogenic"),
            _make_variant(gene_symbol="TP53", acmg_classification="VUS"),
        ]
        s = _build_summary(variants)
        assert s["top_gene"] == "BRCA1"
        assert s["top_acmg_class"] == "Pathogenic"

    def test_empty_variants(self):
        s = _build_summary([])
        assert s["retained"] == 0
        assert s["excluded"] == 0
        assert s["top_gene"] is None
        assert s["top_acmg_class"] is None

    def test_all_excluded(self):
        variants = [_make_variant(excluded=True), _make_variant(excluded=True)]
        s = _build_summary(variants)
        assert s["retained"] == 0
        assert s["top_gene"] is None


# ---------------------------------------------------------------------------
# write_json
# ---------------------------------------------------------------------------

class TestWriteJson:
    def test_creates_file(self, tmp_path):
        variants = [_make_variant()]
        out = str(tmp_path / "out.json")
        result = write_json(variants, out, {"sample_id": "S1"})
        assert Path(result).exists()

    def test_returns_path(self, tmp_path):
        out = str(tmp_path / "out.json")
        result = write_json([_make_variant()], out, {"sample_id": "S1"})
        assert result == out

    def test_json_structure(self, tmp_path):
        retained = _make_variant()
        excluded = _make_variant(excluded=True, exclusion_reason="benign")
        out = str(tmp_path / "out.json")
        write_json([retained, excluded], out, {"sample_id": "S1", "case_id": "C1"})
        with open(out) as f:
            payload = json.load(f)
        assert payload["pipeline_version"] == PIPELINE_VERSION
        assert payload["disclaimer"] == MANDATORY_DISCLAIMER
        assert "generated_at" in payload
        assert len(payload["retained_variants"]) == 1
        assert len(payload["excluded_variants"]) == 1
        assert payload["retained_variants"][0]["rank"] == 1
        assert payload["excluded_variants"][0]["rank"] is None

    def test_disclaimer_exact_text(self, tmp_path):
        out = str(tmp_path / "out.json")
        write_json([_make_variant()], out, {"sample_id": "S1"})
        with open(out) as f:
            payload = json.load(f)
        assert payload["disclaimer"] == (
            "Top genomic candidate(s) requiring expert review. "
            "Variant pathogenicity alone does not prove cause of death."
        )

    def test_creates_parent_dirs(self, tmp_path):
        out = str(tmp_path / "deep" / "nested" / "out.json")
        write_json([_make_variant()], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_metadata_passed_through(self, tmp_path):
        out = str(tmp_path / "out.json")
        meta = {"sample_id": "S1", "case_id": "C99", "input_vcf": "test.vcf"}
        write_json([_make_variant()], out, meta)
        with open(out) as f:
            payload = json.load(f)
        assert payload["metadata"]["case_id"] == "C99"

    def test_empty_variants_list(self, tmp_path):
        out = str(tmp_path / "out.json")
        write_json([], out, {"sample_id": "S1"})
        with open(out) as f:
            payload = json.load(f)
        assert payload["retained_variants"] == []
        assert payload["excluded_variants"] == []
        assert payload["summary"]["total_variants"] == 0

    def test_ranking_order(self, tmp_path):
        v1 = _make_variant(gene_symbol="TOP")
        v2 = _make_variant(gene_symbol="SECOND")
        out = str(tmp_path / "out.json")
        write_json([v1, v2], out, {"sample_id": "S1"})
        with open(out) as f:
            payload = json.load(f)
        assert payload["retained_variants"][0]["rank"] == 1
        assert payload["retained_variants"][1]["rank"] == 2
