"""
Tests for analysis/input_router.py

Covers:
  - InputType detection by extension (no filesystem access)
  - route_input() with valid VCF and VCF.GZ
  - route_input() with FASTQ → NotImplementedError
  - route_input() with BAM/CRAM → NotImplementedError
  - route_input() with unknown extensions → ValueError
  - Error messages are informative

Run: pytest tests/analysis/test_input_router.py -v
"""

import gzip
import pytest
from pathlib import Path

from analysis.input_router import (
    InputType,
    detect_input_type,
    route_input,
)

# ---------------------------------------------------------------------------
# Helpers (shared with test_vcf_validator — kept local for test independence)
# ---------------------------------------------------------------------------

VALID_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
)


def make_plain_vcf(path: Path) -> Path:
    path.write_text(VALID_VCF_HEADER, encoding="utf-8")
    return path


def make_gz_vcf(path: Path) -> Path:
    with gzip.open(str(path), "wt", encoding="utf-8") as fh:
        fh.write(VALID_VCF_HEADER)
    return path


def make_tbi(vcf_gz_path: Path) -> Path:
    tbi = Path(str(vcf_gz_path) + ".tbi")
    tbi.write_bytes(b"\x1f\x8b")
    return tbi


# ---------------------------------------------------------------------------
# detect_input_type — extension detection only, no filesystem access
# ---------------------------------------------------------------------------

class TestDetectInputType:
    # VCF
    def test_plain_vcf(self):
        assert detect_input_type("sample.vcf") == InputType.VCF

    def test_plain_vcf_uppercase(self):
        assert detect_input_type("SAMPLE.VCF") == InputType.VCF

    def test_plain_vcf_with_path(self):
        assert detect_input_type("/data/cases/sample.vcf") == InputType.VCF

    # VCF.GZ
    def test_vcf_gz(self):
        assert detect_input_type("sample.vcf.gz") == InputType.VCF_GZ

    def test_vcf_bgz(self):
        assert detect_input_type("sample.vcf.bgz") == InputType.VCF_GZ

    def test_bcf(self):
        assert detect_input_type("sample.bcf") == InputType.VCF_GZ

    def test_vcf_gz_uppercase(self):
        assert detect_input_type("SAMPLE.VCF.GZ") == InputType.VCF_GZ

    def test_vcf_gz_with_path(self):
        assert detect_input_type("/data/cases/sample.vcf.gz") == InputType.VCF_GZ

    # FASTQ
    def test_fastq(self):
        assert detect_input_type("reads.fastq") == InputType.FASTQ

    def test_fq(self):
        assert detect_input_type("reads.fq") == InputType.FASTQ

    def test_fastq_gz(self):
        assert detect_input_type("reads.fastq.gz") == InputType.FASTQ

    def test_fq_gz(self):
        assert detect_input_type("reads.fq.gz") == InputType.FASTQ

    # BAM / CRAM
    def test_bam(self):
        assert detect_input_type("sample.bam") == InputType.BAM_CRAM

    def test_cram(self):
        assert detect_input_type("sample.cram") == InputType.BAM_CRAM

    # Unknown
    def test_txt_is_unknown(self):
        assert detect_input_type("sample.txt") == InputType.UNKNOWN

    def test_no_extension_is_unknown(self):
        assert detect_input_type("samplefile") == InputType.UNKNOWN

    def test_gz_without_vcf_is_unknown(self):
        # A bare .gz file with no .vcf prefix is unrecognised
        assert detect_input_type("archive.gz") == InputType.UNKNOWN

    def test_bed_is_unknown(self):
        assert detect_input_type("regions.bed") == InputType.UNKNOWN


# ---------------------------------------------------------------------------
# route_input — VCF paths (require real files via tmp_path)
# ---------------------------------------------------------------------------

class TestRouteInputVcf:
    def test_valid_plain_vcf_returns_path(self, tmp_path):
        f = make_plain_vcf(tmp_path / "sample.vcf")
        result = route_input(str(f))
        assert result == str(f)

    def test_valid_vcf_gz_returns_path(self, tmp_path):
        f = make_gz_vcf(tmp_path / "sample.vcf.gz")
        make_tbi(f)
        result = route_input(str(f))
        assert result == str(f)

    def test_missing_vcf_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            route_input(str(tmp_path / "missing.vcf"))

    def test_vcf_gz_without_tbi_raises_value_error(self, tmp_path):
        f = make_gz_vcf(tmp_path / "sample.vcf.gz")
        with pytest.raises(ValueError, match="Tabix index not found"):
            route_input(str(f))

    def test_paths_config_forwarded_without_error(self, tmp_path):
        f = make_plain_vcf(tmp_path / "sample.vcf")
        result = route_input(str(f), paths_config={"vep": {"path": "/usr/bin/vep"}})
        assert result == str(f)


# ---------------------------------------------------------------------------
# route_input — FASTQ: scaffolded, must raise NotImplementedError
# ---------------------------------------------------------------------------

class TestRouteInputFastq:
    """
    FASTQ input is not supported in version 1.
    These tests confirm the correct scaffolding is in place:
    the code exists, recognises FASTQ, and raises NotImplementedError
    with a beginner-friendly message — not a generic crash.
    """

    def test_fastq_raises_not_implemented(self, tmp_path):
        fake = tmp_path / "reads.fastq"
        fake.write_text("@read1\nACGT\n+\nIIII\n")
        with pytest.raises(NotImplementedError):
            route_input(str(fake))

    def test_fastq_gz_raises_not_implemented(self, tmp_path):
        fake = tmp_path / "reads.fastq.gz"
        fake.write_bytes(b"placeholder")
        with pytest.raises(NotImplementedError):
            route_input(str(fake))

    def test_fq_raises_not_implemented(self, tmp_path):
        fake = tmp_path / "reads.fq"
        fake.write_text("@read1\nACGT\n+\nIIII\n")
        with pytest.raises(NotImplementedError):
            route_input(str(fake))

    def test_fastq_error_message_mentions_version_1(self, tmp_path):
        fake = tmp_path / "reads.fastq"
        fake.write_text("@read1\nACGT\n+\nIIII\n")
        with pytest.raises(NotImplementedError, match="version 1"):
            route_input(str(fake))

    def test_fastq_error_message_suggests_vcf(self, tmp_path):
        fake = tmp_path / "reads.fastq"
        fake.write_text("@read1\nACGT\n+\nIIII\n")
        with pytest.raises(NotImplementedError, match="VCF"):
            route_input(str(fake))


# ---------------------------------------------------------------------------
# route_input — BAM/CRAM: scaffolded, must raise NotImplementedError
# ---------------------------------------------------------------------------

class TestRouteInputBamCram:
    """
    BAM/CRAM input is not supported in version 1.
    Same scaffolding contract as FASTQ tests above.
    """

    def test_bam_raises_not_implemented(self, tmp_path):
        fake = tmp_path / "sample.bam"
        fake.write_bytes(b"BAM\x01")
        with pytest.raises(NotImplementedError):
            route_input(str(fake))

    def test_cram_raises_not_implemented(self, tmp_path):
        fake = tmp_path / "sample.cram"
        fake.write_bytes(b"CRAM")
        with pytest.raises(NotImplementedError):
            route_input(str(fake))

    def test_bam_error_message_mentions_version_1(self, tmp_path):
        fake = tmp_path / "sample.bam"
        fake.write_bytes(b"BAM\x01")
        with pytest.raises(NotImplementedError, match="version 1"):
            route_input(str(fake))

    def test_bam_error_message_suggests_vcf(self, tmp_path):
        fake = tmp_path / "sample.bam"
        fake.write_bytes(b"BAM\x01")
        with pytest.raises(NotImplementedError, match="VCF"):
            route_input(str(fake))


# ---------------------------------------------------------------------------
# route_input — unknown / unsupported extensions
# ---------------------------------------------------------------------------

class TestRouteInputUnknown:
    def test_txt_extension_raises_value_error(self, tmp_path):
        fake = tmp_path / "sample.txt"
        fake.write_text("data")
        with pytest.raises(ValueError, match="Unrecognised"):
            route_input(str(fake))

    def test_no_extension_raises_value_error(self, tmp_path):
        fake = tmp_path / "samplefile"
        fake.write_text("data")
        with pytest.raises(ValueError, match="Unrecognised"):
            route_input(str(fake))

    def test_error_message_lists_supported_types(self, tmp_path):
        fake = tmp_path / "sample.txt"
        fake.write_text("data")
        with pytest.raises(ValueError, match=r"\.vcf"):
            route_input(str(fake))
