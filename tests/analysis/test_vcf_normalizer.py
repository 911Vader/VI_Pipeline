"""
Tests for analysis/vcf_normalizer.py

Strategy:
  - _get_bcftools() and _get_reference_fasta() are pure config lookups → tested directly.
  - normalize_vcf() calls subprocess → tested with unittest.mock.patch.

No real bcftools binary or FASTA is used.

Run: pytest tests/analysis/test_vcf_normalizer.py -v
"""

import subprocess
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from analysis.vcf_normalizer import (
    normalize_vcf,
    _get_bcftools,
    _get_reference_fasta,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_paths(bcftools_path="/usr/bin/bcftools", fasta_path="/data/GRCh38.fa") -> dict:
    def entry(path):
        if path:
            return {"path": path, "status": "found"}
        return {"path": None, "status": "not_provided"}

    return {
        "tools": {
            "bcftools": entry(bcftools_path),
        },
        "resources": {
            "reference_fasta": entry(fasta_path),
        },
    }


MINIMAL_PATHS = _make_paths()


# ---------------------------------------------------------------------------
# _get_bcftools
# ---------------------------------------------------------------------------

class TestGetBcftools:
    def test_returns_path_from_paths_config(self):
        result = _get_bcftools(MINIMAL_PATHS)
        assert result == "/usr/bin/bcftools"

    def test_falls_back_to_path_when_not_in_config(self):
        cfg = _make_paths(bcftools_path=None)
        with patch("shutil.which", return_value="/usr/local/bin/bcftools"):
            result = _get_bcftools(cfg)
        assert result == "/usr/local/bin/bcftools"

    def test_raises_when_not_in_config_and_not_on_path(self):
        cfg = _make_paths(bcftools_path=None)
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="bcftools"):
                _get_bcftools(cfg)

    def test_error_message_contains_install_hint(self):
        cfg = _make_paths(bcftools_path=None)
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="conda"):
                _get_bcftools(cfg)

    def test_error_message_mentions_setup(self):
        cfg = _make_paths(bcftools_path=None)
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="python cli.py --setup"):
                _get_bcftools(cfg)

    def test_prefers_config_over_path_lookup(self):
        """When paths_config has bcftools, shutil.which must not be called."""
        with patch("shutil.which") as mock_which:
            result = _get_bcftools(MINIMAL_PATHS)
        mock_which.assert_not_called()
        assert result == "/usr/bin/bcftools"


# ---------------------------------------------------------------------------
# _get_reference_fasta
# ---------------------------------------------------------------------------

class TestGetReferenceFasta:
    def test_returns_fasta_from_paths_config(self):
        result = _get_reference_fasta(MINIMAL_PATHS)
        assert result == "/data/GRCh38.fa"

    def test_raises_when_fasta_missing(self):
        cfg = _make_paths(fasta_path=None)
        with pytest.raises(RuntimeError, match="reference FASTA"):
            _get_reference_fasta(cfg)

    def test_raises_when_status_not_found(self):
        cfg = {"resources": {"reference_fasta": {"path": "/data/GRCh38.fa", "status": "missing"}}}
        with pytest.raises(RuntimeError, match="reference FASTA"):
            _get_reference_fasta(cfg)

    def test_error_message_mentions_setup(self):
        cfg = _make_paths(fasta_path=None)
        with pytest.raises(RuntimeError, match="python cli.py --setup"):
            _get_reference_fasta(cfg)


# ---------------------------------------------------------------------------
# normalize_vcf — mocked subprocess
# ---------------------------------------------------------------------------

class TestNormalizeVcf:
    def test_success_returns_output_path(self, tmp_path):
        out = tmp_path / "sample.norm.vcf"
        out.write_text("##fileformat=VCFv4.2\n")  # simulate bcftools output

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

        assert result == str(out)

    def test_success_calls_bcftools_norm(self, tmp_path):
        out = tmp_path / "sample.norm.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)
            cmd = mock_run.call_args[0][0]

        assert "bcftools" in cmd[0]
        assert "norm" in cmd

    def test_command_contains_split_multiallelic_flag(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)
            cmd = mock_run.call_args[0][0]

        assert "-m" in cmd
        assert "-any" in cmd

    def test_command_contains_fasta_flag(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)
            cmd = mock_run.call_args[0][0]

        assert "-f" in cmd
        assert "/data/GRCh38.fa" in cmd

    def test_command_contains_output_flag(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)
            cmd = mock_run.call_args[0][0]

        assert "-o" in cmd
        assert str(out) in cmd

    def test_command_contains_input_vcf(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)
            cmd = mock_run.call_args[0][0]

        assert "/input/sample.vcf" in cmd

    def test_creates_output_directory(self, tmp_path):
        out = tmp_path / "subdir" / "nested" / "out.vcf"
        assert not out.parent.exists()

        def _create_output(*args, **kwargs):
            out.write_text("##fileformat=VCFv4.2\n")
            return MagicMock(returncode=0, stdout="", stderr="")

        with patch("subprocess.run", side_effect=_create_output):
            normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

        assert out.parent.exists()

    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1, stdout="", stderr="ERROR: could not open reference\n"
            )
            with pytest.raises(RuntimeError, match="bcftools norm failed"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_nonzero_exit_includes_stderr_in_error(self, tmp_path):
        out = tmp_path / "out.vcf"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1, stdout="", stderr="ERROR: reference mismatch at chr1\n"
            )
            with pytest.raises(RuntimeError, match="reference mismatch"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_nonzero_exit_suggests_setup(self, tmp_path):
        out = tmp_path / "out.vcf"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error")
            with pytest.raises(RuntimeError, match="python cli.py --setup"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_timeout_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd_stub = ["/usr/bin/bcftools"]

        with patch("subprocess.run",
                   side_effect=subprocess.TimeoutExpired(cmd_stub, 7200)):
            with pytest.raises(RuntimeError, match="timed out"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_timeout_error_mentions_split_by_chromosome(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd_stub = ["/usr/bin/bcftools"]

        with patch("subprocess.run",
                   side_effect=subprocess.TimeoutExpired(cmd_stub, 7200)):
            with pytest.raises(RuntimeError, match="bcftools view"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_file_not_found_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"

        with patch("subprocess.run", side_effect=FileNotFoundError("not found")):
            with pytest.raises(RuntimeError, match="not found at runtime"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_output_not_produced_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        # bcftools exits 0 but does not create the output file

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            with pytest.raises(RuntimeError, match="output file was not produced"):
                normalize_vcf("/input/sample.vcf", str(out), MINIMAL_PATHS)

    def test_missing_bcftools_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        cfg = _make_paths(bcftools_path=None)

        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="bcftools"):
                normalize_vcf("/input/sample.vcf", str(out), cfg)

    def test_missing_fasta_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        cfg = _make_paths(fasta_path=None)

        with pytest.raises(RuntimeError, match="reference FASTA"):
            normalize_vcf("/input/sample.vcf", str(out), cfg)
