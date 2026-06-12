"""
Tests for analysis/upstream/fastq_handler.py

All subprocess calls are mocked — no real alignment or variant calling occurs.
Tests cover:
  - Tool validation (bwa, samtools missing → RuntimeError)
  - Resource validation (fasta missing, BWA index missing)
  - Docker helper (_get_docker)
  - BWA alignment failure (non-zero exit, timeout, FileNotFoundError)
  - samtools index failure
  - DeepVariant failure (non-zero exit, timeout, FileNotFoundError, output missing)
  - Output file not produced
  - Successful paired-end run
  - Successful single-end run

Run: pytest tests/analysis/upstream/test_fastq_handler.py -v
"""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest

from analysis.upstream.fastq_handler import (
    run_fastq_to_vcf,
    _require_tool,
    _get_docker,
    _get_reference_fasta,
    _get_bwa_index_prefix,
    _run_bwa_align,
    _run_samtools_index,
    _run_deepvariant,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_paths(
    docker_path="/usr/local/bin/docker",
    fasta_path="/data/GRCh38.fa",
    bwa_index_prefix=None,
) -> dict:
    def entry(path, status=None):
        if path:
            return {"path": path, "status": status or "found"}
        return {"path": None, "status": "not_provided"}

    resources = {
        "reference_fasta": entry(fasta_path),
    }
    if bwa_index_prefix:
        resources["bwa_index_prefix"] = entry(bwa_index_prefix)
    else:
        resources["bwa_index_prefix"] = {"path": None, "status": "not_provided"}

    return {
        "tools": {
            "docker": entry(docker_path),
        },
        "resources": resources,
    }


MINIMAL_PATHS = _make_paths()
MINIMAL_SETTINGS = {"vep": {"fork": 4}}

_OK_PROC = MagicMock(returncode=0, stdout=b"", stderr=b"")


# ---------------------------------------------------------------------------
# _require_tool
# ---------------------------------------------------------------------------

class TestRequireTool:
    def test_found_tool_returns_path(self):
        with patch("shutil.which", return_value="/usr/bin/bwa"):
            result = _require_tool("bwa")
        assert result == "/usr/bin/bwa"

    def test_missing_tool_raises_runtime_error(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="not found on PATH"):
                _require_tool("bwa")

    def test_error_mentions_tool_name(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="bwa"):
                _require_tool("bwa")

    def test_error_suggests_conda_install(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="conda install"):
                _require_tool("samtools")


# ---------------------------------------------------------------------------
# _get_docker
# ---------------------------------------------------------------------------

class TestGetDocker:
    def test_returns_path_from_paths_config(self):
        result = _get_docker(MINIMAL_PATHS)
        assert result == "/usr/local/bin/docker"

    def test_falls_back_to_which_when_not_in_config(self):
        with patch("shutil.which", return_value="/usr/bin/docker"):
            result = _get_docker({})
        assert result == "/usr/bin/docker"


# ---------------------------------------------------------------------------
# _get_reference_fasta
# ---------------------------------------------------------------------------

class TestGetReferenceFasta:
    def test_returns_path_from_paths_config(self):
        result = _get_reference_fasta(MINIMAL_PATHS)
        assert result == "/data/GRCh38.fa"

    def test_missing_fasta_raises_runtime_error(self):
        cfg = _make_paths(fasta_path=None)
        with pytest.raises(RuntimeError, match="reference FASTA"):
            _get_reference_fasta(cfg)

    def test_error_suggests_setup(self):
        cfg = _make_paths(fasta_path=None)
        with pytest.raises(RuntimeError, match="python cli.py --setup"):
            _get_reference_fasta(cfg)


# ---------------------------------------------------------------------------
# _get_bwa_index_prefix
# ---------------------------------------------------------------------------

class TestGetBwaIndexPrefix:
    def test_explicit_prefix_in_paths_config(self):
        cfg = _make_paths(bwa_index_prefix="/data/bwa_index/genome.fa")
        result = _get_bwa_index_prefix(cfg, "/data/GRCh38.fa")
        assert result == "/data/bwa_index/genome.fa"

    def test_falls_back_to_fasta_when_index_files_present(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"genome.fa{ext}").write_bytes(b"index")
        result = _get_bwa_index_prefix(MINIMAL_PATHS, str(fasta))
        assert result == str(fasta)

    def test_missing_index_raises_runtime_error(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        with pytest.raises(RuntimeError, match="BWA index not found"):
            _get_bwa_index_prefix(MINIMAL_PATHS, str(fasta))

    def test_missing_index_error_shows_build_command(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        with pytest.raises(RuntimeError, match="bwa index"):
            _get_bwa_index_prefix(MINIMAL_PATHS, str(fasta))

    def test_partial_index_raises_runtime_error(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        for ext in [".amb", ".ann"]:
            (tmp_path / f"genome.fa{ext}").write_bytes(b"x")
        with pytest.raises(RuntimeError, match="BWA index not found"):
            _get_bwa_index_prefix(MINIMAL_PATHS, str(fasta))

    def test_status_not_found_in_config_falls_back_to_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"genome.fa{ext}").write_bytes(b"x")
        cfg = _make_paths(fasta_path=str(fasta))
        cfg["resources"]["bwa_index_prefix"] = {
            "path": "/data/index/genome.fa",
            "status": "missing",
        }
        result = _get_bwa_index_prefix(cfg, str(fasta))
        assert result == str(fasta)

    def test_empty_path_in_config_falls_back_to_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"genome.fa{ext}").write_bytes(b"x")
        cfg = _make_paths(fasta_path=str(fasta))
        cfg["resources"]["bwa_index_prefix"] = {"path": None, "status": "found"}
        result = _get_bwa_index_prefix(cfg, str(fasta))
        assert result == str(fasta)

    def test_found_prefix_is_used_not_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"genome.fa{ext}").write_bytes(b"x")
        explicit_prefix = "/data/separate_index/genome.fa"
        cfg = _make_paths(fasta_path=str(fasta), bwa_index_prefix=explicit_prefix)
        result = _get_bwa_index_prefix(cfg, str(fasta))
        assert result == explicit_prefix


# ---------------------------------------------------------------------------
# _run_bwa_align
# ---------------------------------------------------------------------------

class TestRunBwaAlign:
    def _call(self, tmp_path, fastq_r2=None, bwa_rc=0, sort_rc=0):
        sorted_bam = str(tmp_path / "s.sorted.bam")
        if bwa_rc == 0 and sort_rc == 0:
            Path(sorted_bam).write_bytes(b"BAM")

        bwa_proc = MagicMock()
        bwa_proc.stdout = MagicMock()
        bwa_proc.returncode = bwa_rc
        bwa_proc.stderr.read.return_value = b"bwa stderr" if bwa_rc != 0 else b""
        bwa_proc.wait.return_value = None

        sort_result = MagicMock()
        sort_result.returncode = sort_rc
        sort_result.stderr = "sort stderr" if sort_rc != 0 else ""

        with patch("subprocess.Popen", return_value=bwa_proc), \
             patch("subprocess.run", return_value=sort_result):
            return _run_bwa_align(
                "/usr/bin/bwa", "/usr/bin/samtools", "/data/genome.fa",
                "R1.fastq.gz", fastq_r2,
                "SAMPLE", sorted_bam, 4,
            )

    def test_successful_paired_end_returns_bam_path(self, tmp_path):
        result = self._call(tmp_path, fastq_r2="R2.fastq.gz")
        assert result.endswith(".sorted.bam")

    def test_successful_single_end_returns_bam_path(self, tmp_path):
        result = self._call(tmp_path, fastq_r2=None)
        assert result.endswith(".sorted.bam")

    def test_bwa_nonzero_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="BWA alignment failed"):
            self._call(tmp_path, bwa_rc=1)

    def test_sort_nonzero_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="samtools sort failed"):
            self._call(tmp_path, sort_rc=1)

    def test_bam_not_produced_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "missing.sorted.bam")
        bwa_proc = MagicMock()
        bwa_proc.stdout = MagicMock()
        bwa_proc.returncode = 0
        bwa_proc.stderr.read.return_value = b""
        bwa_proc.wait.return_value = None
        sort_result = MagicMock(returncode=0, stderr="")
        with patch("subprocess.Popen", return_value=bwa_proc), \
             patch("subprocess.run", return_value=sort_result):
            with pytest.raises(RuntimeError, match="sorted BAM was not produced"):
                _run_bwa_align(
                    "/usr/bin/bwa", "/usr/bin/samtools", "/data/genome.fa",
                    "R1.fastq.gz", None, "SAMPLE", sorted_bam, 4,
                )

    def test_timeout_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "s.sorted.bam")
        bwa_proc = MagicMock()
        bwa_proc.stdout = MagicMock()
        bwa_proc.returncode = 0
        bwa_proc.stderr.read.return_value = b""
        bwa_proc.wait.return_value = None
        with patch("subprocess.Popen", return_value=bwa_proc), \
             patch("subprocess.run", side_effect=subprocess.TimeoutExpired("bwa", 100)):
            with pytest.raises(RuntimeError, match="timed out"):
                _run_bwa_align(
                    "/usr/bin/bwa", "/usr/bin/samtools", "/data/genome.fa",
                    "R1.fastq.gz", None, "SAMPLE", sorted_bam, 4,
                )

    def test_file_not_found_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "s.sorted.bam")
        with patch("subprocess.Popen", side_effect=FileNotFoundError("bwa not found")):
            with pytest.raises(RuntimeError, match="not found at runtime"):
                _run_bwa_align(
                    "/nonexistent/bwa", "/usr/bin/samtools", "/data/genome.fa",
                    "R1.fastq.gz", None, "SAMPLE", sorted_bam, 4,
                )


# ---------------------------------------------------------------------------
# _run_samtools_index
# ---------------------------------------------------------------------------

class TestRunSamtoolsIndex:
    def test_success_does_not_raise(self):
        result = MagicMock(returncode=0, stderr="")
        with patch("subprocess.run", return_value=result):
            _run_samtools_index("/usr/bin/samtools", "sample.sorted.bam")

    def test_nonzero_exit_raises_runtime_error(self):
        result = MagicMock(returncode=1, stderr="error indexing")
        with patch("subprocess.run", return_value=result):
            with pytest.raises(RuntimeError, match="samtools index failed"):
                _run_samtools_index("/usr/bin/samtools", "sample.sorted.bam")


# ---------------------------------------------------------------------------
# _run_deepvariant
# ---------------------------------------------------------------------------

class TestRunDeepVariant:
    def _call(self, tmp_path, rc=0, write_vcf=True,
              raise_timeout=False, raise_fnf=False):
        out_dir = str(tmp_path)
        fasta_path = "/data/genome.fa"
        bam_path = str(tmp_path / "SAMPLE.sorted.bam")
        out_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")

        if write_vcf and rc == 0 and not raise_timeout and not raise_fnf:
            Path(out_vcf).write_bytes(b"BGZF")

        mock_result = MagicMock(returncode=rc, stderr="dv error" if rc != 0 else "")

        if raise_fnf:
            side_effect = FileNotFoundError("docker not found")
        elif raise_timeout:
            side_effect = subprocess.TimeoutExpired("docker", 100)
        else:
            side_effect = None

        with patch("subprocess.run",
                   side_effect=side_effect,
                   return_value=mock_result if side_effect is None else None):
            return _run_deepvariant(
                "/usr/bin/docker", fasta_path, bam_path,
                out_dir, "SAMPLE", 4, "WGS", "1.10.0", 43200,
            )

    def test_success_returns_vcf_gz_path(self, tmp_path):
        result = self._call(tmp_path)
        assert result.endswith(".dv.vcf.gz")

    def test_success_path_is_within_out_dir(self, tmp_path):
        result = self._call(tmp_path)
        assert str(tmp_path) in result

    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="DeepVariant exited with code"):
            self._call(tmp_path, rc=1)

    def test_timeout_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="timed out"):
            self._call(tmp_path, raise_timeout=True)

    def test_docker_not_found_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="Docker executable not found"):
            self._call(tmp_path, raise_fnf=True)

    def test_output_not_produced_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="output VCF was not produced"):
            self._call(tmp_path, write_vcf=False)

    def test_timeout_message_mentions_settings_yaml(self, tmp_path):
        with pytest.raises(RuntimeError, match="settings.yaml"):
            self._call(tmp_path, raise_timeout=True)


# ---------------------------------------------------------------------------
# run_fastq_to_vcf — full pipeline (all sub-steps mocked)
# ---------------------------------------------------------------------------

class TestRunFastqToVcf:
    def _run_success(self, tmp_path, fastq_r2=None):
        raw_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")
        sorted_bam = str(tmp_path / "SAMPLE.sorted.bam")
        Path(sorted_bam).write_bytes(b"BAM")
        Path(raw_vcf).write_bytes(b"BGZF")

        with patch("shutil.which", return_value="/usr/bin/bwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_bwa_index_prefix", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._run_bwa_align", return_value=sorted_bam), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_deepvariant", return_value=raw_vcf):
            return run_fastq_to_vcf(
                fastq_r1="R1.fastq.gz",
                fastq_r2=fastq_r2,
                sample_id="SAMPLE",
                output_dir=str(tmp_path),
                paths_config=MINIMAL_PATHS,
                settings=MINIMAL_SETTINGS,
            )

    def test_paired_end_returns_vcf_gz_path(self, tmp_path):
        result = self._run_success(tmp_path, fastq_r2="R2.fastq.gz")
        assert result.endswith(".vcf.gz")

    def test_single_end_returns_vcf_gz_path(self, tmp_path):
        result = self._run_success(tmp_path, fastq_r2=None)
        assert result.endswith(".vcf.gz")

    def test_bwa_missing_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._require_tool",
                   side_effect=RuntimeError("'bwa' is required for FASTQ processing but was not found on PATH")):
            with pytest.raises(RuntimeError, match="bwa"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_samtools_missing_raises_runtime_error(self, tmp_path):
        def require_side_effect(name):
            if name == "samtools":
                raise RuntimeError(f"'{name}' is required for FASTQ processing but was not found on PATH")
            return f"/usr/bin/{name}"
        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=require_side_effect):
            with pytest.raises(RuntimeError, match="samtools"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_missing_fasta_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta",
                   side_effect=RuntimeError("reference FASTA is required")):
            with pytest.raises(RuntimeError, match="reference FASTA"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_missing_bwa_index_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_bwa_index_prefix",
                   side_effect=RuntimeError("BWA index not found")):
            with pytest.raises(RuntimeError, match="BWA index"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_alignment_failure_propagates(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_bwa_index_prefix", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._run_bwa_align",
                   side_effect=RuntimeError("BWA alignment failed (exit code 1)")):
            with pytest.raises(RuntimeError, match="BWA alignment failed"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_variant_calling_failure_propagates(self, tmp_path):
        sorted_bam = str(tmp_path / "SAMPLE.sorted.bam")
        Path(sorted_bam).write_bytes(b"BAM")
        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_bwa_index_prefix", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._run_bwa_align", return_value=sorted_bam), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_deepvariant",
                   side_effect=RuntimeError("DeepVariant exited with code 1")):
            with pytest.raises(RuntimeError, match="DeepVariant"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_output_dir_created(self, tmp_path):
        new_dir = str(tmp_path / "deep" / "nested")
        raw_vcf_path = str(Path(new_dir) / "SAMPLE.dv.vcf.gz")
        sorted_bam_path = str(Path(new_dir) / "SAMPLE.sorted.bam")

        with patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_bwa_index_prefix", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._run_bwa_align", return_value=sorted_bam_path), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_deepvariant", return_value=raw_vcf_path):
            run_fastq_to_vcf(
                "R1.fastq.gz", None, "SAMPLE", new_dir,
                MINIMAL_PATHS, MINIMAL_SETTINGS,
            )
        assert Path(new_dir).exists()
