"""
Tests for gui/validation.py

No display required — all functions are pure.
Run: pytest tests/gui/test_validation.py -v
"""

import pytest
from pathlib import Path
from gui.validation import validate_inputs


def _base(tmp_path):
    """Return a minimal valid inputs dict using tmp_path for real files (VCF input)."""
    f = tmp_path / "input.vcf.gz"
    f.write_text("")
    g = tmp_path / "genes.txt"
    g.write_text("BRCA1\n")
    return {
        "input":      str(f),
        "input2":     "",
        "gene_list":  str(g),
        "sample_id":  "S001",
        "case_id":    "",
    }


def _base_fastq(tmp_path):
    """Return a minimal valid inputs dict with a FASTQ input file."""
    f = tmp_path / "input.fastq.gz"
    f.write_text("")
    g = tmp_path / "genes.txt"
    g.write_text("BRCA1\n")
    return {
        "input":      str(f),
        "input2":     "",
        "gene_list":  str(g),
        "sample_id":  "S001",
        "case_id":    "",
    }


class TestValidInputs:
    def test_all_valid_returns_ok(self, tmp_path):
        ok, msg = validate_inputs(_base(tmp_path))
        assert ok is True
        assert msg == ""

    def test_optional_fields_absent_is_ok(self, tmp_path):
        d = _base(tmp_path)
        d["case_id"] = ""
        d["input2"] = ""
        ok, _ = validate_inputs(d)
        assert ok is True

    def test_with_case_id_is_ok(self, tmp_path):
        d = _base(tmp_path)
        d["case_id"] = "CASE-42"
        ok, _ = validate_inputs(d)
        assert ok is True

    def test_with_input2_present_and_existing_is_ok(self, tmp_path):
        d = _base_fastq(tmp_path)
        r2 = tmp_path / "r2.fastq.gz"
        r2.write_text("")
        d["input2"] = str(r2)
        ok, _ = validate_inputs(d)
        assert ok is True

    def test_fastq_input_without_r2_is_ok(self, tmp_path):
        ok, _ = validate_inputs(_base_fastq(tmp_path))
        assert ok is True


class TestMissingRequiredFields:
    def test_missing_input_file(self, tmp_path):
        d = _base(tmp_path)
        d["input"] = ""
        ok, msg = validate_inputs(d)
        assert not ok
        assert "Input file" in msg

    def test_missing_gene_list(self, tmp_path):
        d = _base(tmp_path)
        d["gene_list"] = ""
        ok, msg = validate_inputs(d)
        assert not ok
        assert "Gene list" in msg

    def test_missing_sample_id(self, tmp_path):
        d = _base(tmp_path)
        d["sample_id"] = ""
        ok, msg = validate_inputs(d)
        assert not ok
        assert "Sample ID" in msg

    def test_whitespace_only_sample_id_is_missing(self, tmp_path):
        d = _base(tmp_path)
        d["sample_id"] = "   "
        ok, msg = validate_inputs(d)
        assert not ok
        assert "Sample ID" in msg


class TestNonexistentPaths:
    def test_nonexistent_input_file(self, tmp_path):
        d = _base(tmp_path)
        d["input"] = str(tmp_path / "no_such.vcf")
        ok, msg = validate_inputs(d)
        assert not ok
        assert "not found" in msg.lower()

    def test_nonexistent_gene_list(self, tmp_path):
        d = _base(tmp_path)
        d["gene_list"] = str(tmp_path / "no_such.txt")
        ok, msg = validate_inputs(d)
        assert not ok
        assert "not found" in msg.lower()

    def test_nonexistent_input2(self, tmp_path):
        d = _base_fastq(tmp_path)
        d["input2"] = str(tmp_path / "no_r2.fastq.gz")
        ok, msg = validate_inputs(d)
        assert not ok
        assert "R2" in msg
        assert "not found" in msg.lower()

    def test_input2_empty_string_does_not_check_existence(self, tmp_path):
        """Empty input2 must not trigger a path-not-found error."""
        d = _base(tmp_path)
        d["input2"] = ""
        ok, _ = validate_inputs(d)
        assert ok is True


class TestErrorMessageContent:
    def test_missing_input_error_names_field(self, tmp_path):
        d = _base(tmp_path)
        d["input"] = ""
        _, msg = validate_inputs(d)
        assert "Input file" in msg

    def test_nonexistent_input_error_includes_path(self, tmp_path):
        d = _base(tmp_path)
        bad = str(tmp_path / "bad.vcf")
        d["input"] = bad
        _, msg = validate_inputs(d)
        assert bad in msg

    def test_nonexistent_input2_error_includes_path(self, tmp_path):
        d = _base_fastq(tmp_path)
        bad = str(tmp_path / "bad_r2.fastq")
        d["input2"] = bad
        _, msg = validate_inputs(d)
        assert bad in msg

    def test_r2_with_vcf_input_is_rejected(self, tmp_path):
        d = _base(tmp_path)
        r2 = tmp_path / "r2.fastq.gz"
        r2.write_text("")
        d["input2"] = str(r2)
        ok, msg = validate_inputs(d)
        assert not ok
        assert "R2" in msg

    def test_unrecognised_extension_is_rejected(self, tmp_path):
        f = tmp_path / "input.txt"
        f.write_text("")
        d = _base(tmp_path)
        d["input"] = str(f)
        ok, msg = validate_inputs(d)
        assert not ok
        assert "Unrecognised" in msg
