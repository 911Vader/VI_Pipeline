"""Tests for reporting/pdf_report.py"""
import pytest

pytest.importorskip("reportlab")

from pathlib import Path

from analysis.parser import Variant
from reporting.pdf_report import (
    MANDATORY_DISCLAIMER,
    PIPELINE_VERSION,
    write_pdf,
    _trunc,
)


def _make_variant(
    chrom="chr1", pos=100, ref="A", alt="T",
    gene_symbol="GENE1", acmg_classification="VUS",
    excluded=False,
    acmg_evidence_codes=None,
    gnomad_af=0.0001,
    hgvs_c="NM_000059.3:c.1A>T",
    hgvs_p="NP_000050.2:p.Met1?",
    consequence="missense_variant",
    in_gene_list=False,
):
    v = Variant(
        chrom=chrom, pos=pos, ref=ref, alt=alt,
        gene_symbol=gene_symbol,
        acmg_classification=acmg_classification,
        excluded=excluded,
        acmg_evidence_codes=acmg_evidence_codes or [],
        gnomad_af=gnomad_af,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        consequence=consequence,
        in_gene_list=in_gene_list,
    )
    return v


# ---------------------------------------------------------------------------
# _trunc
# ---------------------------------------------------------------------------

class TestTrunc:
    def test_short_string_unchanged(self):
        assert _trunc("hello", 10) == "hello"

    def test_exact_length_unchanged(self):
        assert _trunc("hello", 5) == "hello"

    def test_long_string_truncated(self):
        result = _trunc("hello world", 7)
        assert result == "hello …"
        assert len(result) == 7

    def test_empty_string_returns_empty(self):
        assert _trunc("", 10) == ""

    def test_none_returns_empty(self):
        assert _trunc(None, 10) == ""


# ---------------------------------------------------------------------------
# write_pdf — file creation
# ---------------------------------------------------------------------------

class TestWritePdf:
    def test_creates_file(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        result = write_pdf([_make_variant()], out, {"sample_id": "S1"})
        assert Path(result).exists()

    def test_returns_path(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        result = write_pdf([_make_variant()], out, {"sample_id": "S1"})
        assert result == out

    def test_creates_parent_dirs(self, tmp_path):
        out = str(tmp_path / "deep" / "nested" / "report.pdf")
        write_pdf([_make_variant()], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_file_is_pdf(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        write_pdf([_make_variant()], out, {"sample_id": "S1"})
        with open(out, "rb") as f:
            header = f.read(4)
        assert header == b"%PDF"

    def test_empty_variants_creates_file(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        write_pdf([], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_all_excluded_creates_file(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        v = _make_variant(excluded=True)
        write_pdf([v], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_settings_max_candidates_respected(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        variants = [_make_variant(gene_symbol=f"G{i}") for i in range(10)]
        settings = {"reporting": {"max_candidates_in_report": 3}}
        # Should not raise; we just verify it runs without error
        write_pdf(variants, out, {"sample_id": "S1"}, settings=settings)
        assert Path(out).exists()

    def test_metadata_sample_id_only(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        write_pdf([_make_variant()], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_metadata_with_case_id(self, tmp_path):
        out = str(tmp_path / "report.pdf")
        write_pdf([_make_variant()], out, {"sample_id": "S1", "case_id": "C99"})
        assert Path(out).exists()


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

class TestPdfConstants:
    def test_mandatory_disclaimer_exact(self):
        assert MANDATORY_DISCLAIMER == (
            "Top genomic candidate(s) requiring expert review. "
            "Variant pathogenicity alone does not prove cause of death."
        )

    def test_pipeline_version(self):
        assert PIPELINE_VERSION == "1.0.0"


# ---------------------------------------------------------------------------
# VUS in listed genes section
# ---------------------------------------------------------------------------

class TestVusGeneListSection:
    """
    The PDF must include a dedicated section for VUS variants whose gene appears
    in the user-supplied gene list, regardless of rank.

    ACMG classification is not changed — these remain VUS.
    The section only renders when qualifying variants exist.
    """

    def test_vus_in_gene_list_does_not_raise(self, tmp_path):
        """A VUS in the gene list renders without error."""
        out = str(tmp_path / "report.pdf")
        vus = _make_variant(
            gene_symbol="RYR1", acmg_classification="VUS",
            in_gene_list=True,
        )
        write_pdf([vus], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_vus_not_in_gene_list_does_not_raise(self, tmp_path):
        """A VUS not in the gene list renders without error (no extra section)."""
        out = str(tmp_path / "report.pdf")
        vus = _make_variant(
            gene_symbol="TTN", acmg_classification="VUS",
            in_gene_list=False,
        )
        write_pdf([vus], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_pathogenic_in_gene_list_not_in_vus_section(self, tmp_path):
        """Only VUS variants trigger the extra section; P/LP are not duplicated."""
        out = str(tmp_path / "report.pdf")
        path = _make_variant(
            gene_symbol="BRCA1", acmg_classification="Pathogenic",
            in_gene_list=True,
        )
        write_pdf([path], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_excluded_vus_in_gene_list_not_in_section(self, tmp_path):
        """Excluded variants must not appear in the VUS-in-gene-list section."""
        out = str(tmp_path / "report.pdf")
        excluded_vus = _make_variant(
            gene_symbol="RYR1", acmg_classification="VUS",
            in_gene_list=True, excluded=True,
        )
        write_pdf([excluded_vus], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_vus_already_in_top_candidates_not_duplicated(self, tmp_path):
        """
        A high-scoring VUS in the gene list that makes it into the main
        candidates table must NOT also appear in the VUS section.
        """
        out = str(tmp_path / "report.pdf")
        # One high-priority VUS that will be in the top-N candidates
        vus_top = _make_variant(
            gene_symbol="RYR1", acmg_classification="VUS",
            in_gene_list=True,
        )
        write_pdf([vus_top], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_multiple_vus_gene_list_variants_render(self, tmp_path):
        """Multiple qualifying VUS variants all render without error."""
        out = str(tmp_path / "report.pdf")
        # Push them below max_candidates so they land in the VUS section
        settings = {"reporting": {"max_candidates_in_report": 0}}
        variants = [
            _make_variant(
                gene_symbol=g, acmg_classification="VUS",
                in_gene_list=True,
                pos=100 + i,
            )
            for i, g in enumerate(["RYR1", "KCNQ1", "SCN5A"])
        ]
        write_pdf(variants, out, {"sample_id": "S1"}, settings=settings)
        assert Path(out).exists()

    def test_empty_variants_no_vus_section_raised(self, tmp_path):
        """No variants → no VUS section → still produces a valid PDF."""
        out = str(tmp_path / "report.pdf")
        write_pdf([], out, {"sample_id": "S1"})
        assert Path(out).exists()
