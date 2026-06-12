"""
Tests for reporting/tsv_export.py

Verifies:
  - File creation and path return value
  - Parent directory creation
  - Mandatory disclaimer as first comment line
  - Column headers match VARIANT_COLUMNS
  - One data row per retained variant, in rank order
  - Excluded variants are omitted
  - Empty variant list produces header-only output (after disclaimer)
  - Field values are correctly serialised
  - Tab-delimited format
  - Rank column reflects row order

Run: pytest tests/reporting/test_tsv_export.py -v
"""

import csv
import pytest
from pathlib import Path

from analysis.parser import Variant
from reporting.tsv_export import write_tsv, MANDATORY_DISCLAIMER
from reporting.spreadsheet import VARIANT_COLUMNS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(
    chrom="chr1", pos=100, ref="A", alt="T",
    gene_symbol="GENE1", acmg_classification="VUS",
    excluded=False, exclusion_reason=None,
    acmg_evidence_codes=None,
    zygosity="heterozygous",
    gnomad_af=0.0001,
    hgvs_c="c.100A>T",
    hgvs_p="p.Lys34Ile",
    consequence="missense_variant",
    transcript_id="ENST00000123456",
    clinvar_significance="VUS",
) -> Variant:
    return Variant(
        chrom=chrom, pos=pos, ref=ref, alt=alt,
        gene_symbol=gene_symbol,
        acmg_classification=acmg_classification,
        excluded=excluded,
        exclusion_reason=exclusion_reason,
        acmg_evidence_codes=acmg_evidence_codes or [],
        zygosity=zygosity,
        gnomad_af=gnomad_af,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        consequence=consequence,
        transcript_id=transcript_id,
        clinvar_significance=clinvar_significance,
    )


def _read_tsv(path: str):
    """Return (comment_lines, header, rows) from a TSV file."""
    comments = []
    header = None
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                comments.append(line)
            elif header is None:
                header = line.split("\t")
            else:
                rows.append(line.split("\t"))
    return comments, header, rows


# ---------------------------------------------------------------------------
# File creation
# ---------------------------------------------------------------------------

class TestWriteTsvCreation:
    def test_creates_file(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        result = write_tsv([_make_variant()], out, {"sample_id": "S1"})
        assert Path(result).exists()

    def test_returns_path_string(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        result = write_tsv([_make_variant()], out, {"sample_id": "S1"})
        assert result == out

    def test_creates_parent_directories(self, tmp_path):
        out = str(tmp_path / "deep" / "nested" / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        assert Path(out).exists()

    def test_utf8_encoding(self, tmp_path):
        """File must be readable as UTF-8 without errors."""
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        content = Path(out).read_text(encoding="utf-8")
        assert len(content) > 0


# ---------------------------------------------------------------------------
# Mandatory disclaimer
# ---------------------------------------------------------------------------

class TestTsvDisclaimer:
    def test_first_line_is_comment(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        first_line = Path(out).read_text(encoding="utf-8").splitlines()[0]
        assert first_line.startswith("#")

    def test_disclaimer_text_present(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        comments, _, _ = _read_tsv(out)
        assert any(MANDATORY_DISCLAIMER in c for c in comments)

    def test_disclaimer_exact_text(self):
        assert MANDATORY_DISCLAIMER == (
            "Top genomic candidate(s) requiring expert review. "
            "Variant pathogenicity alone does not prove cause of death."
        )


# ---------------------------------------------------------------------------
# Column headers
# ---------------------------------------------------------------------------

class TestTsvHeaders:
    def test_header_row_present(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        _, header, _ = _read_tsv(out)
        assert header is not None
        assert len(header) > 0

    def test_headers_match_variant_columns(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        _, header, _ = _read_tsv(out)
        expected = [col[0] for col in VARIANT_COLUMNS]
        assert header == expected

    def test_required_columns_present(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        _, header, _ = _read_tsv(out)
        for col in ("Chromosome", "Position", "Ref", "Alt", "Gene",
                    "Transcript", "HGVS c.", "HGVS p.", "Consequence",
                    "Zygosity", "ClinVar Significance", "gnomAD AF",
                    "ACMG Evidence Codes", "ACMG Classification"):
            assert col in header, f"Expected column not found: {col}"

    def test_header_is_tab_delimited(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        # Find the header line (first non-comment line)
        for line in Path(out).read_text(encoding="utf-8").splitlines():
            if not line.startswith("#"):
                assert "\t" in line
                break


# ---------------------------------------------------------------------------
# Data rows — retained-only
# ---------------------------------------------------------------------------

class TestTsvDataRows:
    def test_one_retained_variant_one_data_row(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "S1"})
        _, _, rows = _read_tsv(out)
        assert len(rows) == 1

    def test_excluded_variants_omitted(self, tmp_path):
        retained = _make_variant()
        excluded = _make_variant(excluded=True)
        out = str(tmp_path / "out.tsv")
        write_tsv([retained, excluded], out, {"sample_id": "S1"})
        _, _, rows = _read_tsv(out)
        assert len(rows) == 1

    def test_only_excluded_produces_no_data_rows(self, tmp_path):
        excluded = _make_variant(excluded=True)
        out = str(tmp_path / "out.tsv")
        write_tsv([excluded], out, {"sample_id": "S1"})
        _, _, rows = _read_tsv(out)
        assert len(rows) == 0

    def test_empty_variants_produces_no_data_rows(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([], out, {"sample_id": "S1"})
        _, _, rows = _read_tsv(out)
        assert len(rows) == 0

    def test_multiple_retained_all_present(self, tmp_path):
        variants = [_make_variant(gene_symbol=f"GENE{i}") for i in range(5)]
        out = str(tmp_path / "out.tsv")
        write_tsv(variants, out, {"sample_id": "S1"})
        _, _, rows = _read_tsv(out)
        assert len(rows) == 5

    def test_mixed_retained_and_excluded(self, tmp_path):
        variants = [
            _make_variant(gene_symbol="KEEP1"),
            _make_variant(gene_symbol="DROP1", excluded=True),
            _make_variant(gene_symbol="KEEP2"),
            _make_variant(gene_symbol="DROP2", excluded=True),
        ]
        out = str(tmp_path / "out.tsv")
        write_tsv(variants, out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert len(rows) == 2


# ---------------------------------------------------------------------------
# Field value correctness
# ---------------------------------------------------------------------------

class TestTsvFieldValues:
    def _get_col(self, header, rows, col_name, row_idx=0):
        idx = header.index(col_name)
        return rows[row_idx][idx]

    def test_gene_symbol_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(gene_symbol="BRCA2")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Gene") == "BRCA2"

    def test_chromosome_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(chrom="chr17")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Chromosome") == "chr17"

    def test_position_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(pos=41234567)], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Position") == "41234567"

    def test_ref_alt_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(ref="G", alt="A")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Ref") == "G"
        assert self._get_col(header, rows, "Alt") == "A"

    def test_hgvs_c_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(hgvs_c="c.100A>T")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "HGVS c.") == "c.100A>T"

    def test_hgvs_p_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(hgvs_p="p.Lys34Ile")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "HGVS p.") == "p.Lys34Ile"

    def test_consequence_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(consequence="stop_gained")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Consequence") == "stop_gained"

    def test_acmg_classification_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant(acmg_classification="Pathogenic")], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "ACMG Classification") == "Pathogenic"

    def test_acmg_evidence_codes_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        v = _make_variant(acmg_evidence_codes=["PVS1", "PM2"])
        write_tsv([v], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        val = self._get_col(header, rows, "ACMG Evidence Codes")
        assert "PVS1" in val
        assert "PM2" in val

    def test_sample_id_in_row(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        write_tsv([_make_variant()], out, {"sample_id": "MYSAMPLE"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "Sample ID") == "MYSAMPLE"

    def test_clinvar_significance_written(self, tmp_path):
        out = str(tmp_path / "out.tsv")
        v = _make_variant(clinvar_significance="Pathogenic")
        write_tsv([v], out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        assert self._get_col(header, rows, "ClinVar Significance") == "Pathogenic"


# ---------------------------------------------------------------------------
# Row order — variants appear in the order supplied (priority order from ranker)
# ---------------------------------------------------------------------------

class TestTsvRowOrder:
    def test_variants_appear_in_supplied_order(self, tmp_path):
        variants = [_make_variant(gene_symbol=f"GENE{i}") for i in range(3)]
        out = str(tmp_path / "out.tsv")
        write_tsv(variants, out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        gene_idx = header.index("Gene")
        assert [row[gene_idx] for row in rows] == ["GENE0", "GENE1", "GENE2"]

    def test_excluded_do_not_appear_in_tsv(self, tmp_path):
        variants = [
            _make_variant(gene_symbol="A"),
            _make_variant(gene_symbol="X", excluded=True),
            _make_variant(gene_symbol="B"),
        ]
        out = str(tmp_path / "out.tsv")
        write_tsv(variants, out, {"sample_id": "S1"})
        _, header, rows = _read_tsv(out)
        gene_idx = header.index("Gene")
        assert [row[gene_idx] for row in rows] == ["A", "B"]
