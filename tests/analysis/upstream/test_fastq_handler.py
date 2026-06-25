"""
Tests for analysis/upstream/fastq_handler.py

All subprocess calls are mocked -- no real alignment or variant calling occurs.
Tests cover:
  - Tool validation (minibwa, samtools missing -> RuntimeError)
  - Resource validation (fasta missing, minibwa index missing)
  - Docker helper (_get_docker)
  - minibwa alignment failure (non-zero exit, timeout, FileNotFoundError)
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
    _get_minibwa_exe,
    _get_minibwa_index_prefix,
    _run_minibwa_align,
    _run_samtools_index,
    _run_deepvariant,
    _detect_gpu,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_paths(
    docker_path="/usr/local/bin/docker",
    fasta_path="/data/GRCh38.fa",
    minibwa_exe_path="/usr/bin/minibwa",
    minibwa_index_prefix=None,
) -> dict:
    def entry(path, status=None):
        if path:
            return {"path": path, "status": status or "found"}
        return {"path": None, "status": "not_provided"}

    resources = {
        "reference_fasta": entry(fasta_path),
        "minibwa_exe": entry(minibwa_exe_path),
    }
    if minibwa_index_prefix:
        resources["minibwa_index_prefix"] = entry(minibwa_index_prefix)
    else:
        resources["minibwa_index_prefix"] = {"path": None, "status": "not_provided"}

    return {
        "tools": {
            "docker": entry(docker_path),
        },
        "resources": resources,
    }


MINIMAL_PATHS = _make_paths()
MINIMAL_SETTINGS = {"vep": {"fork": 4}}


# ---------------------------------------------------------------------------
# _require_tool
# ---------------------------------------------------------------------------

class TestRequireTool:
    def test_found_tool_returns_path(self):
        with patch("shutil.which", return_value="/usr/bin/samtools"):
            result = _require_tool("samtools")
        assert result == "/usr/bin/samtools"

    def test_missing_tool_raises_runtime_error(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="not found on PATH"):
                _require_tool("samtools")

    def test_error_mentions_tool_name(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="samtools"):
                _require_tool("samtools")

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
# _get_minibwa_exe
# ---------------------------------------------------------------------------

class TestGetMinibwaExe:
    def test_returns_path_from_paths_config(self):
        cfg = _make_paths(minibwa_exe_path="/tools/minibwa/minibwa")
        with patch("pathlib.Path.is_file", return_value=True):
            result = _get_minibwa_exe(cfg)
        assert result == "/tools/minibwa/minibwa"

    def test_falls_back_to_which(self):
        cfg = _make_paths(minibwa_exe_path=None)
        with patch("shutil.which", return_value="/usr/local/bin/minibwa"):
            result = _get_minibwa_exe(cfg)
        assert result == "/usr/local/bin/minibwa"

    def test_missing_raises_runtime_error(self):
        cfg = _make_paths(minibwa_exe_path=None)
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="minibwa"):
                _get_minibwa_exe(cfg)

    def test_error_suggests_build_from_source(self):
        cfg = _make_paths(minibwa_exe_path=None)
        with patch("shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="github.com/lh3/minibwa"):
                _get_minibwa_exe(cfg)


# ---------------------------------------------------------------------------
# _get_minibwa_index_prefix
# ---------------------------------------------------------------------------

class TestGetMinibwaIndexPrefix:
    def test_explicit_prefix_in_paths_config(self):
        cfg = _make_paths(minibwa_index_prefix="/data/minibwa_index/genome")
        result = _get_minibwa_index_prefix(cfg, "/data/GRCh38.fa")
        assert result == "/data/minibwa_index/genome"

    def test_falls_back_to_fasta_when_mbw_present(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        (tmp_path / "genome.fa.mbw").write_bytes(b"index")
        result = _get_minibwa_index_prefix(MINIMAL_PATHS, str(fasta))
        assert result == str(fasta)

    def test_missing_index_raises_runtime_error(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        with pytest.raises(RuntimeError, match="minibwa index not found"):
            _get_minibwa_index_prefix(MINIMAL_PATHS, str(fasta))

    def test_missing_index_error_shows_build_command(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        with pytest.raises(RuntimeError, match="minibwa index -l"):
            _get_minibwa_index_prefix(MINIMAL_PATHS, str(fasta))

    def test_status_not_found_in_config_falls_back_to_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        (tmp_path / "genome.fa.mbw").write_bytes(b"x")
        cfg = _make_paths(fasta_path=str(fasta))
        cfg["resources"]["minibwa_index_prefix"] = {
            "path": "/data/index/genome",
            "status": "missing",
        }
        result = _get_minibwa_index_prefix(cfg, str(fasta))
        assert result == str(fasta)

    def test_empty_path_in_config_falls_back_to_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        (tmp_path / "genome.fa.mbw").write_bytes(b"x")
        cfg = _make_paths(fasta_path=str(fasta))
        cfg["resources"]["minibwa_index_prefix"] = {"path": None, "status": "found"}
        result = _get_minibwa_index_prefix(cfg, str(fasta))
        assert result == str(fasta)

    def test_found_prefix_is_used_not_fasta(self, tmp_path):
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        (tmp_path / "genome.fa.mbw").write_bytes(b"x")
        explicit_prefix = "/data/separate_index/genome"
        cfg = _make_paths(fasta_path=str(fasta), minibwa_index_prefix=explicit_prefix)
        result = _get_minibwa_index_prefix(cfg, str(fasta))
        assert result == explicit_prefix


# ---------------------------------------------------------------------------
# _run_minibwa_align
# ---------------------------------------------------------------------------

class TestRunMinibwaAlign:
    def _call(self, tmp_path, fastq_r2=None, align_rc=0, sort_rc=0):
        sorted_bam = str(tmp_path / "s.sorted.bam")

        align_proc = MagicMock()
        align_proc.stdout = MagicMock()
        align_proc.returncode = align_rc
        align_proc.stderr.read.return_value = b"minibwa stderr" if align_rc != 0 else b""
        align_proc.wait.return_value = None

        sort_result = MagicMock()
        sort_result.returncode = sort_rc
        sort_result.stderr = "sort stderr" if sort_rc != 0 else ""

        if align_rc == 0 and sort_rc == 0:
            Path(sorted_bam).write_bytes(b"BAM")

        with patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("subprocess.Popen", return_value=align_proc), \
             patch("subprocess.run", return_value=sort_result):
            return _run_minibwa_align(
                "/usr/bin/minibwa", "/usr/bin/samtools", "/data/genome",
                "R1.fastq.gz", fastq_r2,
                "SAMPLE", sorted_bam, 4,
            )

    def test_successful_paired_end_returns_bam_path(self, tmp_path):
        result = self._call(tmp_path, fastq_r2="R2.fastq.gz")
        assert result.endswith(".sorted.bam")

    def test_successful_single_end_returns_bam_path(self, tmp_path):
        result = self._call(tmp_path, fastq_r2=None)
        assert result.endswith(".sorted.bam")

    def test_aligner_nonzero_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="minibwa alignment failed"):
            self._call(tmp_path, align_rc=1)

    def test_sort_nonzero_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="samtools sort failed"):
            self._call(tmp_path, sort_rc=1)

    def test_bam_not_produced_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "missing.sorted.bam")
        align_proc = MagicMock()
        align_proc.stdout = MagicMock()
        align_proc.returncode = 0
        align_proc.stderr.read.return_value = b""
        align_proc.wait.return_value = None
        sort_result = MagicMock(returncode=0, stderr="")
        with patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("subprocess.Popen", return_value=align_proc), \
             patch("subprocess.run", return_value=sort_result):
            with pytest.raises(RuntimeError, match="sorted BAM was not produced"):
                _run_minibwa_align(
                    "/usr/bin/minibwa", "/usr/bin/samtools", "/data/genome",
                    "R1.fastq.gz", None, "SAMPLE", sorted_bam, 4,
                )

    def test_timeout_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "s.sorted.bam")
        align_proc = MagicMock()
        align_proc.stdout = MagicMock()
        align_proc.returncode = 0
        align_proc.stderr.read.return_value = b""
        align_proc.wait.return_value = None
        with patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("subprocess.Popen", return_value=align_proc), \
             patch("subprocess.run", side_effect=subprocess.TimeoutExpired("sort", 100)):
            with pytest.raises(RuntimeError, match="timed out"):
                _run_minibwa_align(
                    "/usr/bin/minibwa", "/usr/bin/samtools", "/data/genome",
                    "R1.fastq.gz", None, "SAMPLE", sorted_bam, 4,
                )

    def test_file_not_found_raises_runtime_error(self, tmp_path):
        sorted_bam = str(tmp_path / "s.sorted.bam")
        with patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("subprocess.Popen", side_effect=FileNotFoundError("minibwa not found")):
            with pytest.raises(RuntimeError, match="not found at runtime"):
                _run_minibwa_align(
                    "/nonexistent/minibwa", "/usr/bin/samtools", "/data/genome",
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

    def test_gpu_image_tag_and_flag_in_command(self, tmp_path):
        out_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")
        Path(out_vcf).write_bytes(b"BGZF")
        captured = []
        with patch("subprocess.run",
                   side_effect=lambda cmd, **kw: captured.append(cmd) or MagicMock(returncode=0, stderr="")):
            _run_deepvariant(
                "/usr/bin/docker", "/data/genome.fa", str(tmp_path / "SAMPLE.sorted.bam"),
                str(tmp_path), "SAMPLE", 4, "WGS", "1.10.0", 0, use_gpu=True,
            )
        cmd = captured[0]
        assert "--gpus" in cmd and "all" in cmd
        assert any("1.10.0-gpu" in t for t in cmd)

    def test_cpu_mode_no_gpu_flag(self, tmp_path):
        out_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")
        Path(out_vcf).write_bytes(b"BGZF")
        captured = []
        with patch("subprocess.run",
                   side_effect=lambda cmd, **kw: captured.append(cmd) or MagicMock(returncode=0, stderr="")):
            _run_deepvariant(
                "/usr/bin/docker", "/data/genome.fa", str(tmp_path / "SAMPLE.sorted.bam"),
                str(tmp_path), "SAMPLE", 4, "WGS", "1.10.0", 0, use_gpu=False,
            )
        cmd = captured[0]
        assert "--gpus" not in cmd
        assert not any("1.10.0-gpu" in t for t in cmd)

    def test_disable_small_model_flag_in_command(self, tmp_path):
        out_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")
        Path(out_vcf).write_bytes(b"BGZF")
        captured = []
        with patch("subprocess.run",
                   side_effect=lambda cmd, **kw: captured.append(cmd) or MagicMock(returncode=0, stderr="")):
            _run_deepvariant(
                "/usr/bin/docker", "/data/genome.fa", str(tmp_path / "SAMPLE.sorted.bam"),
                str(tmp_path), "SAMPLE", 4, "WGS", "1.10.0", 0, disable_small_model=True,
            )
        assert "--disable_small_model" in captured[0]

    def test_disable_small_model_false_omits_flag(self, tmp_path):
        out_vcf = str(tmp_path / "SAMPLE.dv.vcf.gz")
        Path(out_vcf).write_bytes(b"BGZF")
        captured = []
        with patch("subprocess.run",
                   side_effect=lambda cmd, **kw: captured.append(cmd) or MagicMock(returncode=0, stderr="")):
            _run_deepvariant(
                "/usr/bin/docker", "/data/genome.fa", str(tmp_path / "SAMPLE.sorted.bam"),
                str(tmp_path), "SAMPLE", 4, "WGS", "1.10.0", 0,
            )
        assert "--disable_small_model" not in captured[0]


# ---------------------------------------------------------------------------
# _detect_gpu
# ---------------------------------------------------------------------------

class TestDetectGpu:
    def test_returns_true_when_nvidia_smi_succeeds(self):
        with patch("shutil.which", return_value="/usr/bin/nvidia-smi"), \
             patch("subprocess.run", return_value=MagicMock(returncode=0, stdout="NVIDIA RTX A6000\n")):
            assert _detect_gpu() is True

    def test_returns_false_when_nvidia_smi_not_found(self):
        with patch("shutil.which", return_value=None):
            assert _detect_gpu() is False

    def test_returns_false_when_nvidia_smi_fails(self):
        with patch("shutil.which", return_value="/usr/bin/nvidia-smi"), \
             patch("subprocess.run", return_value=MagicMock(returncode=1, stdout="")):
            assert _detect_gpu() is False

    def test_returns_false_when_nvidia_smi_raises(self):
        with patch("shutil.which", return_value="/usr/bin/nvidia-smi"), \
             patch("subprocess.run", side_effect=Exception("timeout")):
            assert _detect_gpu() is False


# ---------------------------------------------------------------------------
# run_fastq_to_vcf -- full pipeline (all sub-steps mocked)
# ---------------------------------------------------------------------------

class TestRunFastqToVcf:
    def _run_success(self, tmp_path, fastq_r2=None):
        out = tmp_path
        sorted_bam = str(out / "SAMPLE.sorted.bam")
        dedup_bam  = str(out / "SAMPLE.dedup.bam")
        raw_vcf    = str(out / "SAMPLE.dv.vcf.gz")

        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_minibwa_index_prefix", return_value="/data/genome"), \
             patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("analysis.upstream.fastq_handler._run_minibwa_align", return_value=sorted_bam), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_markdup", return_value=dedup_bam), \
             patch("analysis.upstream.fastq_handler._run_flagstat"), \
             patch("analysis.upstream.fastq_handler._detect_gpu", return_value=False), \
             patch("analysis.upstream.fastq_handler._run_deepvariant", return_value=raw_vcf):
            return run_fastq_to_vcf(
                fastq_r1="R1.fastq.gz",
                fastq_r2=fastq_r2,
                sample_id="SAMPLE",
                output_dir=str(out),
                paths_config=MINIMAL_PATHS,
                settings=MINIMAL_SETTINGS,
            )

    def test_paired_end_returns_vcf_gz_path(self, tmp_path):
        result = self._run_success(tmp_path, fastq_r2="R2.fastq.gz")
        assert result.endswith(".vcf.gz")

    def test_single_end_returns_vcf_gz_path(self, tmp_path):
        result = self._run_success(tmp_path, fastq_r2=None)
        assert result.endswith(".vcf.gz")

    def test_minibwa_missing_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._get_minibwa_exe",
                   side_effect=RuntimeError("minibwa executable not found")):
            with pytest.raises(RuntimeError, match="minibwa"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_samtools_missing_raises_runtime_error(self, tmp_path):
        def require_side_effect(name):
            if name == "samtools":
                raise RuntimeError(
                    f"'{name}' is required for FASTQ processing but was not found on PATH"
                )
            return f"/usr/bin/{name}"

        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=require_side_effect):
            with pytest.raises(RuntimeError, match="samtools"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_missing_fasta_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta",
                   side_effect=RuntimeError("reference FASTA is required")):
            with pytest.raises(RuntimeError, match="reference FASTA"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_missing_minibwa_index_raises_runtime_error(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_minibwa_index_prefix",
                   side_effect=RuntimeError("minibwa index not found")):
            with pytest.raises(RuntimeError, match="minibwa index"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_alignment_failure_propagates(self, tmp_path):
        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_minibwa_index_prefix", return_value="/data/genome"), \
             patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("analysis.upstream.fastq_handler._run_minibwa_align",
                   side_effect=RuntimeError("minibwa alignment failed (exit code 1)")):
            with pytest.raises(RuntimeError, match="minibwa alignment failed"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_variant_calling_failure_propagates(self, tmp_path):
        sorted_bam = str(tmp_path / "SAMPLE.sorted.bam")
        dedup_bam  = str(tmp_path / "SAMPLE.dedup.bam")
        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_minibwa_index_prefix", return_value="/data/genome"), \
             patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("analysis.upstream.fastq_handler._run_minibwa_align", return_value=sorted_bam), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_markdup", return_value=dedup_bam), \
             patch("analysis.upstream.fastq_handler._run_flagstat"), \
             patch("analysis.upstream.fastq_handler._detect_gpu", return_value=False), \
             patch("analysis.upstream.fastq_handler._run_deepvariant",
                   side_effect=RuntimeError("DeepVariant exited with code 1")):
            with pytest.raises(RuntimeError, match="DeepVariant"):
                run_fastq_to_vcf(
                    "R1.fastq.gz", None, "SAMPLE", str(tmp_path),
                    MINIMAL_PATHS, MINIMAL_SETTINGS,
                )

    def test_output_dir_created(self, tmp_path):
        new_dir = str(tmp_path / "deep" / "nested")
        raw_vcf_path    = str(Path(new_dir) / "SAMPLE.dv.vcf.gz")
        sorted_bam_path = str(Path(new_dir) / "SAMPLE.sorted.bam")
        dedup_bam_path  = str(Path(new_dir) / "SAMPLE.dedup.bam")

        with patch("analysis.upstream.fastq_handler._get_minibwa_exe", return_value="/usr/bin/minibwa"), \
             patch("analysis.upstream.fastq_handler._require_tool", side_effect=lambda n: f"/usr/bin/{n}"), \
             patch("analysis.upstream.fastq_handler._get_docker", return_value="/usr/local/bin/docker"), \
             patch("analysis.upstream.fastq_handler._get_reference_fasta", return_value="/data/GRCh38.fa"), \
             patch("analysis.upstream.fastq_handler._get_minibwa_index_prefix", return_value="/data/genome"), \
             patch("analysis.upstream.fastq_handler.shutil.which", return_value=None), \
             patch("analysis.upstream.fastq_handler._run_minibwa_align", return_value=sorted_bam_path), \
             patch("analysis.upstream.fastq_handler._run_samtools_index"), \
             patch("analysis.upstream.fastq_handler._run_markdup", return_value=dedup_bam_path), \
             patch("analysis.upstream.fastq_handler._run_flagstat"), \
             patch("analysis.upstream.fastq_handler._detect_gpu", return_value=False), \
             patch("analysis.upstream.fastq_handler._run_deepvariant", return_value=raw_vcf_path):
            run_fastq_to_vcf(
                "R1.fastq.gz", None, "SAMPLE", new_dir,
                MINIMAL_PATHS, MINIMAL_SETTINGS,
            )
        assert Path(new_dir).exists()
