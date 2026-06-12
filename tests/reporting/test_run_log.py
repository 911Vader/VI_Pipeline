"""Tests for reporting/run_log.py"""
import json
from pathlib import Path

import pytest

from reporting.run_log import RunLog, init_run_log, PIPELINE_VERSION


def _basic_metadata():
    return {"sample_id": "S001", "case_id": "C001"}


# ---------------------------------------------------------------------------
# RunLog accumulation
# ---------------------------------------------------------------------------

class TestRunLogAccumulation:
    def test_add_step_recorded(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_step("parse_vcf", "completed", detail="100 variants")
        d = log.to_dict()
        assert len(d["steps_completed"]) == 1
        assert d["steps_completed"][0]["step"] == "parse_vcf"
        assert d["steps_completed"][0]["status"] == "completed"
        assert d["steps_completed"][0]["detail"] == "100 variants"

    def test_add_step_no_detail(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_step("filter", "completed")
        d = log.to_dict()
        assert "detail" not in d["steps_completed"][0]

    def test_add_warning(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_warning("SpliceAI scores not available")
        d = log.to_dict()
        assert "SpliceAI scores not available" in d["warnings"]

    def test_add_error(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_error("VCF parse failed")
        d = log.to_dict()
        assert "VCF parse failed" in d["errors"]

    def test_add_output_file(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_output_file("xlsx", "results/S001.xlsx")
        d = log.to_dict()
        assert d["output_files"]["xlsx"] == "results/S001.xlsx"

    def test_set_summary(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.set_summary({"retained": 5, "excluded": 2})
        d = log.to_dict()
        assert d["summary"]["retained"] == 5

    def test_multiple_steps_ordered(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        log.add_step("step_a")
        log.add_step("step_b")
        d = log.to_dict()
        steps = [s["step"] for s in d["steps_completed"]]
        assert steps == ["step_a", "step_b"]


# ---------------------------------------------------------------------------
# RunLog structure / required fields
# ---------------------------------------------------------------------------

class TestRunLogStructure:
    def test_to_dict_required_keys(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        d = log.to_dict()
        required = {
            "run_id", "pipeline_version", "started_at", "completed_at",
            "metadata", "steps_completed", "output_files", "summary",
            "warnings", "errors",
        }
        assert required.issubset(d.keys())

    def test_pipeline_version_set(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        d = log.to_dict()
        assert d["pipeline_version"] == PIPELINE_VERSION

    def test_run_id_is_8_char_hex(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        d = log.to_dict()
        run_id = d["run_id"]
        assert len(run_id) == 8
        assert all(c in "0123456789abcdef" for c in run_id)

    def test_metadata_preserved(self, tmp_path):
        meta = {"sample_id": "S1", "case_id": "C99"}
        log = RunLog(str(tmp_path / "run.log.json"), meta)
        d = log.to_dict()
        assert d["metadata"]["case_id"] == "C99"

    def test_completed_at_none_before_finalize(self, tmp_path):
        log = RunLog(str(tmp_path / "run.log.json"), _basic_metadata())
        d = log.to_dict()
        assert d["completed_at"] is None


# ---------------------------------------------------------------------------
# finalize() — writes file
# ---------------------------------------------------------------------------

class TestRunLogFinalize:
    def test_writes_json_file(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = RunLog(path, _basic_metadata())
        result = log.finalize()
        assert Path(result).exists()

    def test_returns_path(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = RunLog(path, _basic_metadata())
        result = log.finalize()
        assert result == path

    def test_json_parseable(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = RunLog(path, _basic_metadata())
        log.add_step("parse_vcf", "completed")
        log.finalize()
        with open(path) as f:
            data = json.load(f)
        assert data["steps_completed"][0]["step"] == "parse_vcf"

    def test_completed_at_set_after_finalize(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = RunLog(path, _basic_metadata())
        log.finalize()
        with open(path) as f:
            data = json.load(f)
        assert data["completed_at"] is not None

    def test_creates_parent_dirs(self, tmp_path):
        path = str(tmp_path / "deep" / "nested" / "run.log.json")
        log = RunLog(path, _basic_metadata())
        log.finalize()
        assert Path(path).exists()


# ---------------------------------------------------------------------------
# Context manager
# ---------------------------------------------------------------------------

class TestRunLogContextManager:
    def test_context_manager_writes_on_exit(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        with RunLog(path, _basic_metadata()) as log:
            log.add_step("parse_vcf")
        assert Path(path).exists()

    def test_exception_captured_as_error(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        try:
            with RunLog(path, _basic_metadata()) as log:
                raise RuntimeError("something broke")
        except RuntimeError:
            pass
        with open(path) as f:
            data = json.load(f)
        assert any("RuntimeError" in e for e in data["errors"])

    def test_exception_not_suppressed(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        with pytest.raises(ValueError):
            with RunLog(path, _basic_metadata()) as log:
                raise ValueError("deliberate")


# ---------------------------------------------------------------------------
# init_run_log factory
# ---------------------------------------------------------------------------

class TestInitRunLog:
    def test_returns_run_log_instance(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = init_run_log(path, _basic_metadata())
        assert isinstance(log, RunLog)

    def test_with_paths_config(self, tmp_path):
        path = str(tmp_path / "run.log.json")
        log = init_run_log(path, _basic_metadata(), paths_config={"ref": "/ref.fa"})
        assert isinstance(log, RunLog)
