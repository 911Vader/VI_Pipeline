"""
Tests for setup/resource_checker.py

Uses tmp_path fixtures so no real files need to exist.
Run: pytest tests/setup/test_resource_checker.py -v
"""

import pytest
from pathlib import Path
from setup.resource_checker import (
    check_resource,
    check_all_resources,
    RESOURCE_SPEC,
    _check_prefix_resource,
    _BWA_INDEX_SUFFIXES,
    _BWA_BWT_MIN_BYTES,
)


class TestCheckResource:
    def test_none_path_required_resource_is_not_provided(self):
        result = check_resource("reference_fasta", None)
        assert result.status == "not_provided"
        assert result.required is True
        assert any("required" in m.lower() or "not provided" in m.lower()
                   for m in result.messages)

    def test_none_path_optional_resource_is_not_provided(self):
        result = check_resource("cadd_snv_file", None)
        assert result.status == "not_provided"
        assert result.required is False

    def test_missing_file_returns_missing(self, tmp_path):
        path = str(tmp_path / "nonexistent.fa")
        result = check_resource("reference_fasta", path)
        assert result.status == "missing"

    def test_found_file_returns_found(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        fai = tmp_path / "genome.fa.fai"
        fai.write_text("chr1\t4\t6\t4\t5\n")
        result = check_resource("reference_fasta", str(fasta))
        assert result.status == "found"

    def test_missing_companion_index_warns(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        # No .fai created
        result = check_resource("reference_fasta", str(fasta))
        assert result.status == "found"
        assert any("tbi" in m.lower() or "fai" in m.lower() or "index" in m.lower()
                   for m in result.messages)

    def test_missing_dir_resource_returns_missing(self, tmp_path):
        result = check_resource("vep_cache_dir", str(tmp_path / "nonexistent_dir"))
        assert result.status == "missing"

    def test_found_dir_resource(self, tmp_path):
        result = check_resource("vep_cache_dir", str(tmp_path))
        assert result.status == "found"

    def test_path_is_file_not_dir_returns_bad_format(self, tmp_path):
        f = tmp_path / "not_a_dir.txt"
        f.write_text("data")
        result = check_resource("vep_cache_dir", str(f))
        assert result.status == "bad_format"


class TestCheckAllResources:
    def test_returns_all_resource_keys(self):
        paths = {k: None for k in RESOURCE_SPEC}
        results = check_all_resources(paths)
        assert set(results.keys()) == set(RESOURCE_SPEC.keys())

    def test_empty_paths_all_not_provided(self):
        results = check_all_resources({})
        for key, result in results.items():
            assert result.status == "not_provided"


# ---------------------------------------------------------------------------
# Optional resource not_provided messages — file format guidance
# ---------------------------------------------------------------------------

class TestOptionalNotProvidedMessages:
    """
    Optional resources that are not provided must emit a message that includes:
    - the paths.yaml key so the user knows exactly where to add the path
    - a file format / download hint from RESOURCE_SPEC.user_hint
    - a re-run instruction
    """

    def _msg(self, resource: str) -> str:
        result = check_resource(resource, None)
        assert result.status == "not_provided"
        assert len(result.messages) == 1
        return result.messages[0]

    # --- paths.yaml key appears in every optional not_provided message ---

    def test_revel_message_contains_paths_yaml_key(self):
        assert "revel_file" in self._msg("revel_file")

    def test_spliceai_snv_message_contains_paths_yaml_key(self):
        assert "spliceai_snv_file" in self._msg("spliceai_snv_file")

    def test_spliceai_indel_message_contains_paths_yaml_key(self):
        assert "spliceai_indel_file" in self._msg("spliceai_indel_file")

    def test_alphamissense_message_contains_paths_yaml_key(self):
        assert "alphamissense_file" in self._msg("alphamissense_file")

    def test_clinvar_message_contains_paths_yaml_key(self):
        assert "clinvar_vcf" in self._msg("clinvar_vcf")

    # --- re-run instruction appears in every optional not_provided message ---

    def test_revel_message_contains_setup_instruction(self):
        assert "python cli.py --setup" in self._msg("revel_file")

    # --- file format hints (from user_hint) appear in the message ---

    def test_revel_message_mentions_file_name(self):
        msg = self._msg("revel_file")
        assert "revel_all_chromosomes" in msg.lower()

    def test_revel_message_mentions_tabix(self):
        msg = self._msg("revel_file")
        assert "tabix" in msg.lower()

    def test_spliceai_snv_message_mentions_hg38(self):
        msg = self._msg("spliceai_snv_file")
        assert "hg38" in msg.lower()

    def test_alphamissense_message_mentions_file_name(self):
        msg = self._msg("alphamissense_file")
        assert "alphamissense" in msg.lower()

    def test_clinvar_message_mentions_ncbi(self):
        msg = self._msg("clinvar_vcf")
        assert "ncbi" in msg.lower() or "clinvar" in msg.lower()

    # --- required resources still have their original concise message ---

    def test_required_resource_message_unchanged(self):
        result = check_resource("reference_fasta", None)
        assert result.status == "not_provided"
        assert result.required is True
        assert any("required" in m.lower() or "not provided" in m.lower()
                   for m in result.messages)

    # --- user_hint fields exist in RESOURCE_SPEC for the five key resources ---

    def test_revel_has_user_hint(self):
        assert RESOURCE_SPEC["revel_file"].get("user_hint")

    def test_spliceai_snv_has_user_hint(self):
        assert RESOURCE_SPEC["spliceai_snv_file"].get("user_hint")

    def test_alphamissense_has_user_hint(self):
        assert RESOURCE_SPEC["alphamissense_file"].get("user_hint")

    def test_clinvar_has_user_hint(self):
        assert RESOURCE_SPEC["clinvar_vcf"].get("user_hint")


# ---------------------------------------------------------------------------
# bwa_index_prefix — prefix resource type
# ---------------------------------------------------------------------------

def _make_real_bwa_index(prefix: str, bwt_size: int = _BWA_BWT_MIN_BYTES + 1024) -> None:
    """Write all five BWA index files at the given prefix with a realistic .bwt size."""
    for suffix in _BWA_INDEX_SUFFIXES:
        p = Path(prefix + suffix)
        if suffix == ".bwt":
            p.write_bytes(b"x" * bwt_size)
        else:
            p.write_bytes(b"header_data" * 10)   # small but non-empty


class TestBwaIndexPrefixResource:
    def test_in_resource_spec(self):
        assert "bwa_index_prefix" in RESOURCE_SPEC

    def test_is_not_required(self):
        assert RESOURCE_SPEC["bwa_index_prefix"].get("required") is False

    def test_not_provided_returns_not_provided(self):
        result = check_resource("bwa_index_prefix", None)
        assert result.status == "not_provided"
        assert result.required is False

    # --- All five files present and .bwt large enough ---

    def test_valid_full_index_returns_found(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix)
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "found"
        assert result.path == prefix

    def test_valid_full_index_ok_message(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix)
        result = check_resource("bwa_index_prefix", prefix)
        assert any("OK" in m for m in result.messages)

    # --- Missing files ---

    def test_no_files_returns_missing(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "missing"

    def test_only_amb_present_returns_missing(self, tmp_path):
        """A single .amb file must not satisfy the check — all 5 are required."""
        prefix = str(tmp_path / "genome.fa")
        Path(prefix + ".amb").write_bytes(b"x" * 100)
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "missing"

    def test_four_of_five_files_returns_missing(self, tmp_path):
        """Missing any one of the five files must return missing."""
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix)
        Path(prefix + ".sa").unlink()   # remove one
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "missing"

    def test_missing_error_lists_missing_filenames(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        result = check_resource("bwa_index_prefix", prefix)
        assert any(".amb" in m for m in result.messages)

    def test_missing_error_contains_bwa_index_command(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        result = check_resource("bwa_index_prefix", prefix)
        assert any("bwa index" in m for m in result.messages)

    # --- Placeholder / too-small .bwt ---

    def test_placeholder_bwt_returns_bad_format(self, tmp_path):
        """All five files exist but .bwt is only a few bytes — should be bad_format."""
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            Path(prefix + suffix).write_bytes(b"placeholder")
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "bad_format"

    def test_bwt_exactly_at_minimum_returns_bad_format(self, tmp_path):
        """A .bwt exactly at the minimum (not above it) is still rejected."""
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            p = Path(prefix + suffix)
            if suffix == ".bwt":
                p.write_bytes(b"x" * _BWA_BWT_MIN_BYTES)   # exactly at threshold
            else:
                p.write_bytes(b"x" * 100)
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "bad_format"

    def test_bwt_one_byte_above_minimum_returns_found(self, tmp_path):
        """One byte above the threshold should pass."""
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix, bwt_size=_BWA_BWT_MIN_BYTES + 1)
        result = check_resource("bwa_index_prefix", prefix)
        assert result.status == "found"

    def test_placeholder_error_message_mentions_size(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            Path(prefix + suffix).write_bytes(b"x")
        result = check_resource("bwa_index_prefix", prefix)
        assert any("bytes" in m for m in result.messages)

    def test_placeholder_error_message_mentions_bwt(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            Path(prefix + suffix).write_bytes(b"x")
        result = check_resource("bwa_index_prefix", prefix)
        assert any(".bwt" in m for m in result.messages)

    # --- Preservation through check_all_resources (regression) ---

    def test_check_all_resources_includes_bwa_index_prefix(self):
        results = check_all_resources({})
        assert "bwa_index_prefix" in results

    def test_bwa_index_prefix_preserved_when_valid(self, tmp_path):
        """
        Regression: bwa_index_prefix must survive check_all_resources() so
        write_paths() writes it back to paths.yaml.
        """
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix)
        results = check_all_resources({"bwa_index_prefix": prefix})
        assert "bwa_index_prefix" in results
        assert results["bwa_index_prefix"].status == "found"
        assert results["bwa_index_prefix"].path == prefix

    def test_bwa_index_prefix_preserved_when_missing(self, tmp_path):
        """Even a missing prefix must be recorded (status: missing) in results."""
        prefix = str(tmp_path / "nonexistent")
        results = check_all_resources({"bwa_index_prefix": prefix})
        assert "bwa_index_prefix" in results
        assert results["bwa_index_prefix"].status == "missing"
        assert results["bwa_index_prefix"].path == prefix

    def test_bwa_index_prefix_preserved_when_bad_format(self, tmp_path):
        """A placeholder index must be recorded as bad_format, not dropped."""
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            Path(prefix + suffix).write_bytes(b"x")
        results = check_all_resources({"bwa_index_prefix": prefix})
        assert "bwa_index_prefix" in results
        assert results["bwa_index_prefix"].status == "bad_format"

    # --- _check_prefix_resource directly ---

    def test_check_prefix_resource_directly_found(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        _make_real_bwa_index(prefix)
        result = _check_prefix_resource("bwa_index_prefix", prefix, {})
        assert result.status == "found"

    def test_check_prefix_resource_directly_missing(self, tmp_path):
        prefix = str(tmp_path / "no_index")
        result = _check_prefix_resource("bwa_index_prefix", prefix, {})
        assert result.status == "missing"

    def test_check_prefix_resource_directly_bad_format(self, tmp_path):
        prefix = str(tmp_path / "genome.fa")
        for suffix in _BWA_INDEX_SUFFIXES:
            Path(prefix + suffix).write_bytes(b"tiny")
        result = _check_prefix_resource("bwa_index_prefix", prefix, {})
        assert result.status == "bad_format"
