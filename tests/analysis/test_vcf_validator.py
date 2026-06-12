"""
Tests for analysis/vcf_validator.py

Uses pytest tmp_path for all file fixtures so no real genomic data is needed.
Compressed fixtures use Python's gzip module, which is sufficient for header
validation tests (bgzip vs gzip distinction is not tested here — that is a
runtime bcftools concern, not a Python-level concern).

Run: pytest tests/analysis/test_vcf_validator.py -v
"""

import gzip
import os
import stat
import pytest
from pathlib import Path

from analysis.vcf_validator import (
    validate_vcf,
    has_csq_annotation,
    has_multiallelic_records,
    _is_compressed,
    _check_exists,
    _check_is_file,
    _check_readable,
    _check_not_empty,
    _check_tabix_index,
    _check_vcf_header,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

VALID_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
)

INVALID_HEADER_CONTENT = (
    "This is not a VCF file\n"
    "chr1\t100\t.\tA\tT\n"
)


def write_plain_vcf(path: Path, content: str = VALID_VCF_HEADER) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def write_gz_vcf(path: Path, content: str = VALID_VCF_HEADER) -> Path:
    with gzip.open(str(path), "wt", encoding="utf-8") as fh:
        fh.write(content)
    return path


def write_tbi(vcf_gz_path: Path) -> Path:
    tbi = Path(str(vcf_gz_path) + ".tbi")
    tbi.write_bytes(b"\x1f\x8b")  # minimal placeholder — just needs to exist
    return tbi


# ---------------------------------------------------------------------------
# _is_compressed
# ---------------------------------------------------------------------------

class TestIsCompressed:
    def test_vcf_gz_is_compressed(self):
        assert _is_compressed("sample.vcf.gz") is True

    def test_vcf_bgz_is_compressed(self):
        assert _is_compressed("sample.vcf.bgz") is True

    def test_plain_vcf_is_not_compressed(self):
        assert _is_compressed("sample.vcf") is False

    def test_case_insensitive(self):
        assert _is_compressed("SAMPLE.VCF.GZ") is True

    def test_bam_is_not_compressed(self):
        assert _is_compressed("sample.bam") is False


# ---------------------------------------------------------------------------
# _check_exists
# ---------------------------------------------------------------------------

class TestCheckExists:
    def test_existing_file_passes(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_text("data")
        _check_exists(str(f))  # should not raise

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            _check_exists(str(tmp_path / "nonexistent.vcf"))


# ---------------------------------------------------------------------------
# _check_is_file
# ---------------------------------------------------------------------------

class TestCheckIsFile:
    def test_regular_file_passes(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_text("data")
        _check_is_file(str(f))  # should not raise

    def test_directory_raises(self, tmp_path):
        with pytest.raises(ValueError, match="not a file"):
            _check_is_file(str(tmp_path))


# ---------------------------------------------------------------------------
# _check_readable
# ---------------------------------------------------------------------------

class TestCheckReadable:
    def test_readable_file_passes(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_text("data")
        _check_readable(str(f))  # should not raise

    @pytest.mark.skipif(os.getuid() == 0, reason="root ignores file permissions")
    def test_unreadable_file_raises(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_text("data")
        f.chmod(0o000)
        try:
            with pytest.raises(ValueError, match="cannot be read"):
                _check_readable(str(f))
        finally:
            f.chmod(0o644)  # restore so tmp_path cleanup works


# ---------------------------------------------------------------------------
# _check_not_empty
# ---------------------------------------------------------------------------

class TestCheckNotEmpty:
    def test_non_empty_file_passes(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_text("data")
        _check_not_empty(str(f))  # should not raise

    def test_empty_file_raises(self, tmp_path):
        f = tmp_path / "sample.vcf"
        f.write_bytes(b"")
        with pytest.raises(ValueError, match="empty"):
            _check_not_empty(str(f))


# ---------------------------------------------------------------------------
# _check_tabix_index
# ---------------------------------------------------------------------------

class TestCheckTabixIndex:
    def test_tbi_present_passes(self, tmp_path):
        vcf_gz = tmp_path / "sample.vcf.gz"
        write_gz_vcf(vcf_gz)
        write_tbi(vcf_gz)
        _check_tabix_index(str(vcf_gz))  # should not raise

    def test_tbi_missing_raises(self, tmp_path):
        vcf_gz = tmp_path / "sample.vcf.gz"
        write_gz_vcf(vcf_gz)
        # No .tbi created
        with pytest.raises(ValueError, match="Tabix index not found"):
            _check_tabix_index(str(vcf_gz))

    def test_error_message_contains_fix_command(self, tmp_path):
        vcf_gz = tmp_path / "sample.vcf.gz"
        write_gz_vcf(vcf_gz)
        with pytest.raises(ValueError, match="tabix"):
            _check_tabix_index(str(vcf_gz))


# ---------------------------------------------------------------------------
# _check_vcf_header
# ---------------------------------------------------------------------------

class TestCheckVcfHeader:
    def test_valid_plain_vcf_header_passes(self, tmp_path):
        f = write_plain_vcf(tmp_path / "sample.vcf")
        _check_vcf_header(str(f), is_gz=False)  # should not raise

    def test_valid_gz_vcf_header_passes(self, tmp_path):
        f = write_gz_vcf(tmp_path / "sample.vcf.gz")
        _check_vcf_header(str(f), is_gz=True)  # should not raise

    def test_invalid_plain_vcf_header_raises(self, tmp_path):
        f = write_plain_vcf(tmp_path / "sample.vcf", content=INVALID_HEADER_CONTENT)
        with pytest.raises(ValueError, match="##fileformat=VCF"):
            _check_vcf_header(str(f), is_gz=False)

    def test_invalid_gz_vcf_header_raises(self, tmp_path):
        f = write_gz_vcf(tmp_path / "sample.vcf.gz", content=INVALID_HEADER_CONTENT)
        with pytest.raises(ValueError, match="##fileformat=VCF"):
            _check_vcf_header(str(f), is_gz=True)

    def test_plain_file_opened_as_gz_raises(self, tmp_path):
        """A plain-text file opened as gzip should raise ValueError."""
        f = tmp_path / "sample.vcf.gz"
        f.write_text(VALID_VCF_HEADER, encoding="utf-8")  # not actually gzipped
        with pytest.raises(ValueError):
            _check_vcf_header(str(f), is_gz=True)


# ---------------------------------------------------------------------------
# validate_vcf — integration of all checks
# ---------------------------------------------------------------------------

class TestValidateVcf:
    def test_valid_plain_vcf_returns_path(self, tmp_path):
        f = write_plain_vcf(tmp_path / "sample.vcf")
        result = validate_vcf(str(f))
        assert result == str(f)

    def test_valid_vcf_gz_returns_path(self, tmp_path):
        f = write_gz_vcf(tmp_path / "sample.vcf.gz")
        write_tbi(f)
        result = validate_vcf(str(f))
        assert result == str(f)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            validate_vcf(str(tmp_path / "missing.vcf"))

    def test_vcf_gz_without_tbi_raises(self, tmp_path):
        f = write_gz_vcf(tmp_path / "sample.vcf.gz")
        with pytest.raises(ValueError, match="Tabix index not found"):
            validate_vcf(str(f))

    def test_empty_vcf_raises(self, tmp_path):
        f = tmp_path / "empty.vcf"
        f.write_bytes(b"")
        with pytest.raises(ValueError, match="empty"):
            validate_vcf(str(f))

    def test_non_vcf_content_raises(self, tmp_path):
        f = write_plain_vcf(tmp_path / "sample.vcf", content=INVALID_HEADER_CONTENT)
        with pytest.raises(ValueError, match="##fileformat=VCF"):
            validate_vcf(str(f))

    def test_paths_config_accepted_and_ignored_in_v1(self, tmp_path):
        """paths_config parameter must be accepted without error in version 1."""
        f = write_plain_vcf(tmp_path / "sample.vcf")
        result = validate_vcf(str(f), paths_config={"vep": {"path": "/usr/bin/vep"}})
        assert result == str(f)


# ---------------------------------------------------------------------------
# has_csq_annotation
# ---------------------------------------------------------------------------

_CSQ_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
)

_NO_CSQ_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
)


class TestHasCsqAnnotation:
    def test_plain_vcf_with_csq_returns_true(self, tmp_path):
        f = tmp_path / "annotated.vcf"
        f.write_text(_CSQ_HEADER, encoding="utf-8")
        assert has_csq_annotation(str(f)) is True

    def test_plain_vcf_without_csq_returns_false(self, tmp_path):
        f = tmp_path / "unannotated.vcf"
        f.write_text(_NO_CSQ_HEADER, encoding="utf-8")
        assert has_csq_annotation(str(f)) is False

    def test_gz_vcf_with_csq_returns_true(self, tmp_path):
        f = tmp_path / "annotated.vcf.gz"
        with gzip.open(str(f), "wt", encoding="utf-8") as fh:
            fh.write(_CSQ_HEADER)
        assert has_csq_annotation(str(f)) is True

    def test_gz_vcf_without_csq_returns_false(self, tmp_path):
        f = tmp_path / "unannotated.vcf.gz"
        with gzip.open(str(f), "wt", encoding="utf-8") as fh:
            fh.write(_NO_CSQ_HEADER)
        assert has_csq_annotation(str(f)) is False

    def test_returns_false_on_unreadable_file(self, tmp_path):
        # nonexistent file — should return False, not raise
        result = has_csq_annotation(str(tmp_path / "nonexistent.vcf"))
        assert result is False

    def test_csq_line_anywhere_in_header(self, tmp_path):
        content = (
            "##fileformat=VCFv4.2\n"
            "##FILTER=<ID=PASS,Description=\"pass\">\n"
            "##contig=<ID=chr1,length=248956422>\n"
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        f = tmp_path / "annotated.vcf"
        f.write_text(content, encoding="utf-8")
        assert has_csq_annotation(str(f)) is True

    def test_stops_reading_at_data_lines(self, tmp_path):
        # CSQ appears after the #CHROM line (not in header) — should return False
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="late">\n'
            "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
        )
        f = tmp_path / "late_csq.vcf"
        f.write_text(content, encoding="utf-8")
        assert has_csq_annotation(str(f)) is False

    def test_empty_file_returns_false(self, tmp_path):
        f = tmp_path / "empty.vcf"
        f.write_bytes(b"")
        assert has_csq_annotation(str(f)) is False


# ---------------------------------------------------------------------------
# has_multiallelic_records
# ---------------------------------------------------------------------------

_BIALLELIC_ONLY = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
    "chr1\t200\t.\tG\tC\t60\tPASS\t.\n"
)

_WITH_MULTIALLELIC = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
    "chr1\t200\t.\tG\tC,T\t60\tPASS\t.\n"   # multi-allelic
    "chr1\t300\t.\tA\tG\t70\tPASS\t.\n"
)

_ALL_MULTIALLELIC = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT,G\t50\tPASS\t.\n"
    "chr1\t200\t.\tC\tA,T,G\t60\tPASS\t.\n"
)


class TestHasMultiallelicRecords:
    def test_returns_false_for_biallelic_only(self, tmp_path):
        f = tmp_path / "biallelic.vcf"
        f.write_text(_BIALLELIC_ONLY, encoding="utf-8")
        assert has_multiallelic_records(str(f)) is False

    def test_returns_true_when_multiallelic_present(self, tmp_path):
        f = tmp_path / "multi.vcf"
        f.write_text(_WITH_MULTIALLELIC, encoding="utf-8")
        assert has_multiallelic_records(str(f)) is True

    def test_returns_true_for_all_multiallelic(self, tmp_path):
        f = tmp_path / "all_multi.vcf"
        f.write_text(_ALL_MULTIALLELIC, encoding="utf-8")
        assert has_multiallelic_records(str(f)) is True

    def test_returns_false_for_empty_file(self, tmp_path):
        f = tmp_path / "empty.vcf"
        f.write_bytes(b"")
        assert has_multiallelic_records(str(f)) is False

    def test_returns_false_for_header_only(self, tmp_path):
        f = tmp_path / "header_only.vcf"
        f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                     encoding="utf-8")
        assert has_multiallelic_records(str(f)) is False

    def test_returns_false_on_missing_file(self, tmp_path):
        assert has_multiallelic_records(str(tmp_path / "missing.vcf")) is False

    def test_gzipped_multiallelic_returns_true(self, tmp_path):
        f = tmp_path / "multi.vcf.gz"
        with gzip.open(str(f), "wt", encoding="utf-8") as fh:
            fh.write(_WITH_MULTIALLELIC)
        assert has_multiallelic_records(str(f)) is True

    def test_gzipped_biallelic_returns_false(self, tmp_path):
        f = tmp_path / "biallelic.vcf.gz"
        with gzip.open(str(f), "wt", encoding="utf-8") as fh:
            fh.write(_BIALLELIC_ONLY)
        assert has_multiallelic_records(str(f)) is False

    def test_max_records_limits_scan(self, tmp_path):
        """max_records=1 should only check the first data line."""
        # First line is biallelic, second is multi-allelic
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"      # biallelic — first data line
            "chr1\t200\t.\tG\tC,T\t60\tPASS\t.\n"    # multi-allelic — second
        )
        f = tmp_path / "mixed.vcf"
        f.write_text(content, encoding="utf-8")
        # With max_records=1 only the first data line is checked → False
        assert has_multiallelic_records(str(f), max_records=1) is False

    def test_three_alt_alleles_detected(self, tmp_path):
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tT,G,C\t50\tPASS\t.\n"
        )
        f = tmp_path / "tri.vcf"
        f.write_text(content, encoding="utf-8")
        assert has_multiallelic_records(str(f)) is True

    def test_indel_alt_without_comma_is_biallelic(self, tmp_path):
        """An indel ALT like 'ATCG' should not trigger multi-allelic detection."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tATCG\t50\tPASS\t.\n"
        )
        f = tmp_path / "indel.vcf"
        f.write_text(content, encoding="utf-8")
        assert has_multiallelic_records(str(f)) is False
