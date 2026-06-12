"""Tests for validation/truth_set.py"""

import io
import pytest
from unittest.mock import patch, mock_open

from validation.truth_set import (
    TruthEntry,
    load_truth_set,
    build_truth_index,
)

# ---------------------------------------------------------------------------
# TruthEntry
# ---------------------------------------------------------------------------

class TestTruthEntry:
    def test_parses_variant_id(self):
        e = TruthEntry(
            variant_id="chr17:43071077:A:AAGAA",
            gene="BRCA1",
            expected_classification="Pathogenic",
        )
        assert e.chrom == "chr17"
        assert e.pos == 43071077
        assert e.ref == "A"
        assert e.alt == "AAGAA"

    def test_match_key(self):
        e = TruthEntry(
            variant_id="chr13:32339659:TTTACA:T",
            gene="BRCA2",
            expected_classification="Pathogenic",
        )
        assert e.match_key == ("chr13", 32339659, "TTTACA", "T")

    def test_invalid_variant_id_does_not_crash(self):
        e = TruthEntry(
            variant_id="NOTVALID",
            gene="BRCA1",
            expected_classification="VUS",
        )
        assert e.chrom == ""
        assert e.pos == 0

    def test_expected_acmg_codes_defaults_empty(self):
        e = TruthEntry(
            variant_id="chr1:100:A:T",
            gene="GENE1",
            expected_classification="VUS",
        )
        assert e.expected_acmg_codes == []


# ---------------------------------------------------------------------------
# load_truth_set
# ---------------------------------------------------------------------------

_VALID_TSV = (
    "variant_id\tgene\thgvs_c\thgvs_p\texpected_classification\texpected_acmg_codes\tcurator_notes\n"
    "chr17:43071077:A:AAGAA\tBRCA1\tc.5266dupC\tp.Gln1756fs\tPathogenic\tPVS1 PS3\tClassic BRCA1\n"
    "chr13:32339659:TTTACA:T\tBRCA2\t\t\tVUS\t\t\n"
    "# comment line — must be ignored\n"
    "chr1:100:A:T\tGENE1\t\t\tBenign\tBA1\tnote\n"
)


def _write_tmp_tsv(tmp_path, content):
    p = tmp_path / "truth.tsv"
    p.write_text(content, encoding="utf-8")
    return str(p)


class TestLoadTruthSet:
    def test_loads_valid_tsv(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        assert len(entries) == 3

    def test_comment_lines_skipped(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        assert all(e.gene != "comment" for e in entries)

    def test_gene_and_classification_parsed(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        assert entries[0].gene == "BRCA1"
        assert entries[0].expected_classification == "Pathogenic"
        assert entries[0].expected_acmg_codes == ["PVS1", "PS3"]

    def test_optional_fields_populated(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        assert entries[0].hgvs_c == "c.5266dupC"
        assert entries[0].curator_notes == "Classic BRCA1"

    def test_empty_optional_fields_ok(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        assert entries[1].hgvs_c == ""
        assert entries[1].expected_acmg_codes == []

    def test_invalid_classification_skipped(self, tmp_path):
        tsv = (
            "variant_id\tgene\texpected_classification\n"
            "chr1:100:A:T\tGENE1\tINVALID_CATEGORY\n"
            "chr1:200:A:T\tGENE2\tPathogenic\n"
        )
        path = _write_tmp_tsv(tmp_path, tsv)
        entries = load_truth_set(path)
        assert len(entries) == 1
        assert entries[0].gene == "GENE2"

    def test_missing_required_column_raises(self, tmp_path):
        tsv = "variant_id\tgene\n" "chr1:100:A:T\tGENE1\n"
        path = _write_tmp_tsv(tmp_path, tsv)
        with pytest.raises(ValueError, match="expected_classification"):
            load_truth_set(path)

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            load_truth_set("/nonexistent/file.tsv")

    def test_all_five_categories_accepted(self, tmp_path):
        rows = "\n".join(
            f"chr1:{i}:A:T\tGENE{i}\t{cls}"
            for i, cls in enumerate([
                "Pathogenic", "Likely Pathogenic", "VUS",
                "Likely Benign", "Benign",
            ], start=100)
        )
        tsv = "variant_id\tgene\texpected_classification\n" + rows + "\n"
        path = _write_tmp_tsv(tmp_path, tsv)
        entries = load_truth_set(path)
        assert len(entries) == 5


# ---------------------------------------------------------------------------
# build_truth_index
# ---------------------------------------------------------------------------

class TestBuildTruthIndex:
    def test_keys_are_match_keys(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        index = build_truth_index(entries)
        assert ("chr17", 43071077, "A", "AAGAA") in index

    def test_lookup_returns_correct_entry(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        index = build_truth_index(entries)
        e = index[("chr17", 43071077, "A", "AAGAA")]
        assert e.gene == "BRCA1"
        assert e.expected_classification == "Pathogenic"

    def test_size_matches_entry_count(self, tmp_path):
        path = _write_tmp_tsv(tmp_path, _VALID_TSV)
        entries = load_truth_set(path)
        index = build_truth_index(entries)
        assert len(index) == len(entries)
