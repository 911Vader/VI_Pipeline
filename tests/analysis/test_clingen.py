"""
Tests for analysis/clingen.py

Verifies:
  - load_clingen parses the CSV correctly
  - Strongest classification is kept when a gene has multiple disease curations
  - Unknown classifications are silently skipped
  - Missing file returns empty dict (no exception)
  - annotate_clingen sets clingen_validity on matched variants
  - Unmatched variants keep clingen_validity=None
  - Return value of annotate_clingen is the count of matched variants

Run: pytest tests/analysis/test_clingen.py -v
"""

import csv
import pytest
from pathlib import Path

from analysis.parser import Variant
from analysis.clingen import load_clingen, annotate_clingen, _CLASSIFICATION_RANK


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_clingen_csv(tmp_path, rows: list[tuple]) -> str:
    """
    Write a minimal ClinGen gene validity CSV with the 4 metadata header lines,
    the column header row, then the supplied data rows.

    Each tuple in `rows` is (gene_symbol, classification).
    The CSV has at least 7 columns; column 0 = gene, column 6 = classification.
    """
    path = tmp_path / "ClinGen_gene_validity.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        # 4 metadata lines
        writer.writerow(["ClinGen Gene Disease Validity Curations"])
        writer.writerow(["Date generated: 2026-04-29"])
        writer.writerow([])
        writer.writerow([])
        # column header
        writer.writerow(["GENE SYMBOL", "GENE ID (HGNC)", "DISEASE LABEL",
                          "DISEASE ID (MONDO)", "MOI", "SOP", "CLASSIFICATION",
                          "ONLINE REPORT", "CLASSIFICATION DATE", "GCEP"])
        for gene, classification in rows:
            writer.writerow([gene, "HGNC:1234", "Disease X",
                             "MONDO:0000001", "AD", "SOP8", classification,
                             "https://example.com", "2024-01-01", "GCEP1"])
    return str(path)


def _make_variant(gene_symbol="BRCA1") -> Variant:
    return Variant(
        chrom="chr17", pos=41234567, ref="A", alt="T",
        gene_symbol=gene_symbol,
    )


# ---------------------------------------------------------------------------
# load_clingen
# ---------------------------------------------------------------------------

class TestLoadClingen:
    def test_returns_empty_dict_for_missing_file(self, tmp_path):
        result = load_clingen(str(tmp_path / "does_not_exist.csv"))
        assert result == {}

    def test_parses_single_row(self, tmp_path):
        csv_path = _write_clingen_csv(tmp_path, [("BRCA1", "Definitive")])
        result = load_clingen(csv_path)
        assert result == {"BRCA1": "Definitive"}

    def test_parses_multiple_distinct_genes(self, tmp_path):
        rows = [("BRCA1", "Definitive"), ("BRCA2", "Strong"), ("TP53", "Moderate")]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert result["BRCA1"] == "Definitive"
        assert result["BRCA2"] == "Strong"
        assert result["TP53"] == "Moderate"

    def test_keeps_strongest_for_gene_with_multiple_diseases(self, tmp_path):
        # BRCA1 appears twice: Moderate then Definitive — Definitive should win
        rows = [("BRCA1", "Moderate"), ("BRCA1", "Definitive")]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert result["BRCA1"] == "Definitive"

    def test_keeps_strongest_regardless_of_row_order(self, tmp_path):
        # Stronger entry comes first
        rows = [("BRCA1", "Definitive"), ("BRCA1", "Limited")]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert result["BRCA1"] == "Definitive"

    def test_skips_unknown_classification(self, tmp_path):
        rows = [("BRCA1", "Definitive"), ("UNKNOWNGENE", "NotARealClass")]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert "UNKNOWNGENE" not in result
        assert result["BRCA1"] == "Definitive"

    def test_all_valid_classifications_accepted(self, tmp_path):
        valid = list(_CLASSIFICATION_RANK.keys())
        rows = [(f"GENE{i}", cls) for i, cls in enumerate(valid)]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert len(result) == len(valid)

    def test_classification_rank_order(self):
        # Definitive (0) < Strong (1) < Moderate (2) < Limited (3) <
        # Disputed (4) < Refuted (5) < No Known Disease Relationship (6)
        assert _CLASSIFICATION_RANK["Definitive"] < _CLASSIFICATION_RANK["Strong"]
        assert _CLASSIFICATION_RANK["Strong"] < _CLASSIFICATION_RANK["Moderate"]
        assert _CLASSIFICATION_RANK["Moderate"] < _CLASSIFICATION_RANK["Limited"]
        assert _CLASSIFICATION_RANK["Limited"] < _CLASSIFICATION_RANK["Disputed"]
        assert _CLASSIFICATION_RANK["Disputed"] < _CLASSIFICATION_RANK["Refuted"]
        assert _CLASSIFICATION_RANK["Refuted"] < _CLASSIFICATION_RANK["No Known Disease Relationship"]

    def test_empty_gene_rows_skipped(self, tmp_path):
        rows = [("", "Definitive"), ("BRCA1", "Strong")]
        csv_path = _write_clingen_csv(tmp_path, rows)
        result = load_clingen(csv_path)
        assert "" not in result
        assert result["BRCA1"] == "Strong"

    def test_returns_dict_not_none(self, tmp_path):
        csv_path = _write_clingen_csv(tmp_path, [])
        result = load_clingen(csv_path)
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# annotate_clingen
# ---------------------------------------------------------------------------

class TestAnnotateClingen:
    def test_matched_variant_gets_classification(self, tmp_path):
        gene_dict = {"BRCA1": "Definitive"}
        v = _make_variant("BRCA1")
        annotate_clingen([v], gene_dict)
        assert v.clingen_validity == "Definitive"

    def test_unmatched_variant_keeps_none(self):
        gene_dict = {"BRCA1": "Definitive"}
        v = _make_variant("TP53")
        annotate_clingen([v], gene_dict)
        assert v.clingen_validity is None

    def test_returns_count_of_matched_variants(self):
        gene_dict = {"BRCA1": "Definitive", "BRCA2": "Strong"}
        variants = [
            _make_variant("BRCA1"),
            _make_variant("BRCA2"),
            _make_variant("TP53"),
        ]
        count = annotate_clingen(variants, gene_dict)
        assert count == 2

    def test_empty_gene_dict_returns_zero(self):
        v = _make_variant("BRCA1")
        count = annotate_clingen([v], {})
        assert count == 0
        assert v.clingen_validity is None

    def test_empty_variant_list_returns_zero(self):
        gene_dict = {"BRCA1": "Definitive"}
        count = annotate_clingen([], gene_dict)
        assert count == 0

    def test_variant_with_no_gene_symbol_unmatched(self):
        gene_dict = {"BRCA1": "Definitive"}
        v = _make_variant("")
        v.gene_symbol = ""
        annotate_clingen([v], gene_dict)
        assert v.clingen_validity is None

    def test_modifies_variants_in_place(self):
        gene_dict = {"BRCA1": "Definitive"}
        v = _make_variant("BRCA1")
        original_id = id(v)
        annotate_clingen([v], gene_dict)
        assert id(v) == original_id
        assert v.clingen_validity == "Definitive"

    def test_multiple_variants_same_gene(self):
        gene_dict = {"BRCA1": "Definitive"}
        variants = [_make_variant("BRCA1"), _make_variant("BRCA1")]
        count = annotate_clingen(variants, gene_dict)
        assert count == 2
        assert all(v.clingen_validity == "Definitive" for v in variants)

    @pytest.mark.parametrize("classification", [
        "Definitive", "Strong", "Moderate", "Limited",
        "Disputed", "Refuted", "No Known Disease Relationship",
    ])
    def test_all_classifications_annotated(self, classification):
        gene_dict = {"GENE1": classification}
        v = _make_variant("GENE1")
        annotate_clingen([v], gene_dict)
        assert v.clingen_validity == classification
