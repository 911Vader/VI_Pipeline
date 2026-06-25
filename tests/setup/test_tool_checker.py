"""
Tests for setup/tool_checker.py

These tests use mocking to avoid requiring real tools to be installed.
Run: pytest tests/setup/test_tool_checker.py -v
"""

import pytest
from unittest.mock import patch, MagicMock
from setup.tool_checker import (
    check_tool,
    check_all_tools,
    _find_executable,
    _parse_version,
    _version_ok,
    ToolCheckResult,
)


class TestFindExecutable:
    def test_finds_existing_tool(self):
        with patch("shutil.which", return_value="/usr/bin/python3"):
            result = _find_executable(["python3", "python"])
        assert result == "/usr/bin/python3"

    def test_returns_none_when_not_found(self):
        with patch("shutil.which", return_value=None):
            result = _find_executable(["nonexistent_tool_xyz"])
        assert result is None

    def test_tries_all_candidates(self):
        def which_side(name):
            return "/usr/bin/bcftools" if name == "bcftools" else None
        with patch("shutil.which", side_effect=which_side):
            result = _find_executable(["notfound", "bcftools"])
        assert result == "/usr/bin/bcftools"


class TestParseVersion:
    def test_parse_python_version(self):
        raw = "Python 3.11.2"
        result = _parse_version(raw, "python")
        assert result == (3, 11)

    def test_parse_bcftools_version(self):
        raw = "bcftools 1.18\nUsing htslib 1.18"
        result = _parse_version(raw, "bcftools")
        assert result == (1, 18)

    def test_parse_docker_version(self):
        raw = "Docker version 24.0.5, build ced0996"
        result = _parse_version(raw, "docker")
        assert result == (24, 0)

    def test_parse_vep_version(self):
        raw = "ensembl-vep : 112.0"
        result = _parse_version(raw, "vep")
        assert result == (112, 0)

    def test_returns_none_for_unparseable(self):
        result = _parse_version("garbage output with no version", "python")
        assert result is None


class TestVersionOk:
    def test_version_meets_minimum(self):
        assert _version_ok((3, 11), (3, 10)) is True

    def test_version_equals_minimum(self):
        assert _version_ok((3, 10), (3, 10)) is True

    def test_version_below_minimum(self):
        assert _version_ok((3, 9), (3, 10)) is False

    def test_major_version_higher(self):
        assert _version_ok((4, 0), (3, 10)) is True


class TestCheckTool:
    def test_found_and_valid(self):
        with (
            patch("shutil.which", return_value="/usr/bin/python3"),
            patch("subprocess.run") as mock_run,
        ):
            mock_run.return_value = MagicMock(
                stdout="Python 3.11.2\n",
                stderr="",
                returncode=0,
            )
            result = check_tool("python")
        assert result.status == "found"
        assert result.path == "/usr/bin/python3"
        assert result.version_ok is True

    def test_not_found_on_path(self):
        with patch("shutil.which", return_value=None):
            result = check_tool("vep")
        assert result.status == "missing"
        assert result.path is None
        assert "not found" in result.message.lower()

    def test_version_too_old(self):
        with (
            patch("shutil.which", return_value="/usr/bin/bcftools"),
            patch("subprocess.run") as mock_run,
        ):
            mock_run.return_value = MagicMock(
                stdout="bcftools 1.10\n",
                stderr="",
                returncode=0,
            )
            result = check_tool("bcftools")
        assert result.status == "version_fail"
        assert result.version_ok is False
        assert "1.17" in result.message or "upgrade" in result.message.lower()

    def test_hint_path_used(self):
        with (
            patch("pathlib.Path.is_file", return_value=True),
            patch("subprocess.run") as mock_run,
        ):
            mock_run.return_value = MagicMock(
                stdout="Python 3.12.0\n",
                stderr="",
                returncode=0,
            )
            result = check_tool("python", hint_path="/opt/custom/python3")
        assert result.path == "/opt/custom/python3"

    def test_check_all_tools_returns_all_keys(self):
        with (
            patch("shutil.which", return_value=None),
        ):
            results = check_all_tools()
        expected_keys = {"docker", "conda", "python", "nextflow", "vep",
                         "bcftools", "bgzip", "tabix", "minibwa", "samtools"}
        assert set(results.keys()) == expected_keys
