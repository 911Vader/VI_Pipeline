"""
Tests for analysis/annotator.py

Strategy:
  - build_vep_command() is pure logic (no subprocess) → tested directly.
  - run_vep() calls subprocess.run → tested with unittest.mock.patch.
  - annotate() orchestrates both → integration test with both mocked.

No real VEP is invoked. All filesystem checks that subprocess would trigger
are avoided by mocking at the subprocess boundary.

Run: pytest tests/analysis/test_annotator.py -v
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from analysis.annotator import (
    build_vep_command,
    run_vep,
    annotate,
    _get_required_path,
    _get_optional_path,
    _build_plugin_args,
    _build_custom_args,
    _resolve_timeout,
    _output_is_complete,
    _count_vcf_records,
)

# ---------------------------------------------------------------------------
# Fixtures — minimal valid paths_config and settings
# ---------------------------------------------------------------------------

def _make_paths(
    vep_path="/opt/vep/vep",
    cache_dir="/data/vep_cache",
    fasta="/data/GRCh38.fa",
    plugins_dir=None,
    loftee_path=None,
    clinvar=None,
    cadd_snv=None,
    cadd_indel=None,
    revel=None,
    spliceai_snv=None,
    spliceai_indel=None,
    alphamissense=None,
    loftee_human_ancestor_fa=None,
    loftee_conservation_file=None,
    loftee_gerp_bigwig=None,
) -> dict:
    """Construct a minimal paths_config dict matching the paths.yaml structure."""

    def entry(path):
        if path:
            return {"path": path, "status": "found"}
        return {"path": None, "status": "not_provided"}

    return {
        "tools": {
            "vep": entry(vep_path),
        },
        "resources": {
            "vep_cache_dir":            entry(cache_dir),
            "reference_fasta":          entry(fasta),
            "vep_plugins_dir":          entry(plugins_dir),
            "loftee_path":              entry(loftee_path),
            "clinvar_vcf":              entry(clinvar),
            "cadd_snv_file":            entry(cadd_snv),
            "cadd_indel_file":          entry(cadd_indel),
            "revel_file":               entry(revel),
            "spliceai_snv_file":        entry(spliceai_snv),
            "spliceai_indel_file":      entry(spliceai_indel),
            "alphamissense_file":       entry(alphamissense),
            "loftee_human_ancestor_fa": entry(loftee_human_ancestor_fa),
            "loftee_conservation_file": entry(loftee_conservation_file),
            "loftee_gerp_bigwig":       entry(loftee_gerp_bigwig),
        },
    }


MINIMAL_PATHS = _make_paths()

MINIMAL_SETTINGS = {
    "vep": {
        "fork": 4,
        "distance": 5000,
        "extra_flags": [],
    }
}


# ---------------------------------------------------------------------------
# _get_required_path
# ---------------------------------------------------------------------------

class TestGetRequiredPath:
    def test_found_path_returned(self):
        cfg = _make_paths()
        result = _get_required_path(cfg, "tools", "vep", "VEP")
        assert result == "/opt/vep/vep"

    def test_missing_path_raises_runtime_error(self):
        cfg = _make_paths(vep_path=None)
        with pytest.raises(RuntimeError, match="required for annotation"):
            _get_required_path(cfg, "tools", "vep", "VEP executable")

    def test_wrong_status_raises_runtime_error(self):
        cfg = {"tools": {"vep": {"path": "/usr/bin/vep", "status": "missing"}}}
        with pytest.raises(RuntimeError, match="required for annotation"):
            _get_required_path(cfg, "tools", "vep", "VEP executable")

    def test_error_message_contains_section_and_key(self):
        cfg = _make_paths(cache_dir=None)
        with pytest.raises(RuntimeError, match="resources.vep_cache_dir"):
            _get_required_path(cfg, "resources", "vep_cache_dir", "VEP cache")

    def test_error_message_contains_setup_instruction(self):
        cfg = _make_paths(vep_path=None)
        with pytest.raises(RuntimeError, match="python cli.py --setup"):
            _get_required_path(cfg, "tools", "vep", "VEP executable")


# ---------------------------------------------------------------------------
# _get_optional_path
# ---------------------------------------------------------------------------

class TestGetOptionalPath:
    def test_found_optional_path_returned(self):
        cfg = _make_paths(clinvar="/data/clinvar.vcf.gz")
        result = _get_optional_path(cfg, "resources", "clinvar_vcf")
        assert result == "/data/clinvar.vcf.gz"

    def test_absent_optional_returns_none(self):
        result = _get_optional_path(MINIMAL_PATHS, "resources", "clinvar_vcf")
        assert result is None

    def test_wrong_status_returns_none(self):
        cfg = {"resources": {"clinvar_vcf": {"path": "/data/clinvar.vcf.gz", "status": "missing"}}}
        result = _get_optional_path(cfg, "resources", "clinvar_vcf")
        assert result is None

    def test_non_dict_entry_returns_none(self):
        cfg = {"resources": {"clinvar_vcf": None}}
        result = _get_optional_path(cfg, "resources", "clinvar_vcf")
        assert result is None

    def test_missing_status_with_path_logs_warning(self, caplog):
        """Configured path with status='missing' must return None AND log a warning."""
        import logging
        cfg = {"resources": {"clinvar_vcf": {"path": "/data/clinvar.vcf.gz", "status": "missing"}}}
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            result = _get_optional_path(cfg, "resources", "clinvar_vcf")
        assert result is None
        assert any("missing" in r.message.lower() for r in caplog.records)

    def test_bad_format_status_with_path_logs_warning(self, caplog):
        """Configured path with status='bad_format' must return None AND log a warning."""
        import logging
        cfg = {"resources": {"revel_file": {"path": "/data/revel.tsv.gz", "status": "bad_format"}}}
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            result = _get_optional_path(cfg, "resources", "revel_file")
        assert result is None
        assert any("bad_format" in r.message for r in caplog.records)

    def test_not_provided_status_no_warning(self, caplog):
        """status='not_provided' with no path must NOT log a warning."""
        import logging
        cfg = {"resources": {"revel_file": {"path": None, "status": "not_provided"}}}
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            _get_optional_path(cfg, "resources", "revel_file")
        assert len(caplog.records) == 0

    def test_warning_includes_path_and_status(self, caplog):
        """The warning message must contain both the path and the status string."""
        import logging
        cfg = {"resources": {"revel_file": {"path": "/my/revel.tsv.gz", "status": "missing"}}}
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            _get_optional_path(cfg, "resources", "revel_file")
        combined = " ".join(r.message for r in caplog.records)
        assert "/my/revel.tsv.gz" in combined
        assert "missing" in combined


# ---------------------------------------------------------------------------
# _build_plugin_args
# ---------------------------------------------------------------------------

class TestBuildPluginArgs:
    def test_no_plugins_returns_empty_list(self):
        args = _build_plugin_args(MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert args == []

    def test_cadd_added_when_both_files_present(self):
        cfg = _make_paths(cadd_snv="/data/cadd_snv.tsv.gz", cadd_indel="/data/cadd_indel.tsv.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert "--plugin" in args
        cadd_idx = args.index("--plugin")
        assert "CADD" in args[cadd_idx + 1]
        assert "/data/cadd_snv.tsv.gz" in args[cadd_idx + 1]
        assert "/data/cadd_indel.tsv.gz" in args[cadd_idx + 1]

    def test_cadd_skipped_when_only_snv_present(self):
        cfg = _make_paths(cadd_snv="/data/cadd_snv.tsv.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert not any("CADD" in a for a in args)

    def test_cadd_skipped_when_only_indel_present(self):
        cfg = _make_paths(cadd_indel="/data/cadd_indel.tsv.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert not any("CADD" in a for a in args)

    def test_revel_added_when_present(self):
        cfg = _make_paths(revel="/data/revel.tsv.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert "--plugin" in args
        revel_idx = [i for i, a in enumerate(args) if a == "--plugin" and "REVEL" in args[i + 1]]
        assert revel_idx, "REVEL plugin not found in args"

    def test_revel_absent_when_not_configured(self):
        args = _build_plugin_args(MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert not any("REVEL" in a for a in args)

    def test_spliceai_added_when_both_files_present(self):
        cfg = _make_paths(
            spliceai_snv="/data/spliceai_snv.vcf.gz",
            spliceai_indel="/data/spliceai_indel.vcf.gz",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("SpliceAI" in a for a in args)

    def test_spliceai_enabled_snv_only_when_indel_absent(self):
        # SNV-only mode: SpliceAI is enabled with just the SNV file.
        # Indel coverage is unavailable but SNV annotations still fire.
        cfg = _make_paths(spliceai_snv="/data/spliceai_snv.vcf.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("SpliceAI" in a for a in args)
        spliceai_arg = next(a for a in args if "SpliceAI" in a)
        assert "snv=" in spliceai_arg
        assert "indel=" not in spliceai_arg

    def test_alphamissense_added_when_present(self):
        cfg = _make_paths(alphamissense="/data/alphamissense.tsv.gz")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("AlphaMissense" in a for a in args)

    def test_loftee_added_when_loftee_path_configured(self):
        cfg = _make_paths(loftee_path="/opt/loftee_grch38")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("LoF" in a for a in args)

    def test_loftee_includes_loftee_path(self):
        cfg = _make_paths(loftee_path="/opt/loftee_grch38")
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        lof_arg = next((a for a in args if "LoF" in a), None)
        assert lof_arg is not None
        assert "loftee_path:/opt/loftee_grch38" in lof_arg

    def test_loftee_includes_human_ancestor_fa_when_present(self):
        cfg = _make_paths(
            loftee_path="/opt/loftee_grch38",
            loftee_human_ancestor_fa="/data/human_ancestor_GRCh38.fa.gz",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        lof_arg = next(a for a in args if "LoF" in a)
        assert "human_ancestor_fa:/data/human_ancestor_GRCh38.fa.gz" in lof_arg

    def test_loftee_includes_conservation_file_when_present(self):
        cfg = _make_paths(
            loftee_path="/opt/loftee_grch38",
            loftee_conservation_file="/data/loftee.sql",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        lof_arg = next(a for a in args if "LoF" in a)
        assert "conservation_file:/data/loftee.sql" in lof_arg

    def test_loftee_includes_gerp_bigwig_when_present(self):
        cfg = _make_paths(
            loftee_path="/opt/loftee_grch38",
            loftee_gerp_bigwig="/data/gerp.bw",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        lof_arg = next(a for a in args if "LoF" in a)
        assert "gerp_bigwig:/data/gerp.bw" in lof_arg

    def test_loftee_absent_when_loftee_path_not_configured(self):
        args = _build_plugin_args(MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert not any("LoF" in a for a in args)

    def test_loftee_warns_when_human_ancestor_fa_missing(self, caplog):
        import logging
        cfg = _make_paths(loftee_path="/opt/loftee_grch38")
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("human_ancestor_fa" in r.message for r in caplog.records)

    def test_loftee_warns_when_conservation_file_missing(self, caplog):
        import logging
        cfg = _make_paths(loftee_path="/opt/loftee_grch38")
        with caplog.at_level(logging.WARNING, logger="analysis.annotator"):
            _build_plugin_args(cfg, MINIMAL_SETTINGS)
        assert any("conservation_file" in r.message for r in caplog.records)

    def test_loftee_all_data_files_in_plugin_string(self):
        cfg = _make_paths(
            loftee_path="/opt/loftee_grch38",
            loftee_human_ancestor_fa="/data/human_ancestor.fa.gz",
            loftee_conservation_file="/data/loftee.sql",
            loftee_gerp_bigwig="/data/gerp.bw",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        lof_arg = next(a for a in args if "LoF" in a)
        assert "loftee_path:/opt/loftee_grch38" in lof_arg
        assert "human_ancestor_fa:/data/human_ancestor.fa.gz" in lof_arg
        assert "conservation_file:/data/loftee.sql" in lof_arg
        assert "gerp_bigwig:/data/gerp.bw" in lof_arg

    def test_multiple_plugins_all_present(self):
        cfg = _make_paths(
            revel="/data/revel.tsv.gz",
            cadd_snv="/data/cadd_snv.tsv.gz",
            cadd_indel="/data/cadd_indel.tsv.gz",
            alphamissense="/data/am.tsv.gz",
        )
        args = _build_plugin_args(cfg, MINIMAL_SETTINGS)
        plugin_vals = [args[i + 1] for i, a in enumerate(args) if a == "--plugin"]
        assert any("CADD" in v for v in plugin_vals)
        assert any("REVEL" in v for v in plugin_vals)
        assert any("AlphaMissense" in v for v in plugin_vals)


# ---------------------------------------------------------------------------
# _build_custom_args
# ---------------------------------------------------------------------------

class TestBuildCustomArgs:
    def test_clinvar_added_when_present(self):
        cfg = _make_paths(clinvar="/data/clinvar.vcf.gz")
        args = _build_custom_args(cfg)
        assert "--custom" in args
        custom_idx = args.index("--custom")
        val = args[custom_idx + 1]
        assert "/data/clinvar.vcf.gz" in val
        assert "ClinVar" in val
        assert "CLNSIG" in val
        assert "CLNREVSTAT" in val
        assert "CLNDN" in val

    def test_clinvar_absent_when_not_configured(self):
        args = _build_custom_args(MINIMAL_PATHS)
        assert "--custom" not in args

    def test_custom_format_is_vcf_exact(self):
        cfg = _make_paths(clinvar="/data/clinvar.vcf.gz")
        args = _build_custom_args(cfg)
        custom_idx = args.index("--custom")
        val = args[custom_idx + 1]
        # Format: file,name,type,match_type,coordsys,fields...
        parts = val.split(",")
        assert parts[1] == "ClinVar"
        assert parts[2] == "vcf"
        assert parts[3] == "exact"


# ---------------------------------------------------------------------------
# build_vep_command — core command structure
# ---------------------------------------------------------------------------

class TestBuildVepCommand:
    def test_command_starts_with_vep_executable(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert cmd[0] == "/opt/vep/vep"

    def test_input_file_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--input_file" in cmd
        assert cmd[cmd.index("--input_file") + 1] == "in.vcf"

    def test_output_file_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--output_file" in cmd
        assert cmd[cmd.index("--output_file") + 1] == "out.vcf"

    def test_vcf_output_format_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--vcf" in cmd

    def test_offline_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--offline" in cmd

    def test_cache_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--cache" in cmd

    def test_assembly_is_grch38(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--assembly" in cmd
        assert cmd[cmd.index("--assembly") + 1] == "GRCh38"

    def test_dir_cache_is_set(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--dir_cache" in cmd
        assert cmd[cmd.index("--dir_cache") + 1] == "/data/vep_cache"

    def test_fasta_is_set(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--fasta" in cmd
        assert cmd[cmd.index("--fasta") + 1] == "/data/GRCh38.fa"

    def test_fork_count_from_settings(self):
        settings = {"vep": {"fork": 8, "distance": 5000, "extra_flags": []}}
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, settings)
        assert "--fork" in cmd
        assert cmd[cmd.index("--fork") + 1] == "8"

    def test_fork_default_when_not_in_settings(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, {})
        assert "--fork" in cmd
        assert cmd[cmd.index("--fork") + 1] == "4"

    def test_hgvs_flags_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--hgvs" in cmd
        assert "--hgvsg" in cmd

    def test_symbol_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--symbol" in cmd

    def test_canonical_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--canonical" in cmd

    def test_af_flags_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--af" in cmd
        assert "--af_gnomad" in cmd
        assert "--max_af" in cmd

    def test_no_stats_flag_present(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--no_stats" in cmd

    def test_plugins_dir_added_when_configured(self):
        cfg = _make_paths(plugins_dir="/opt/vep/plugins")
        cmd = build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)
        assert "--dir_plugins" in cmd
        assert cmd[cmd.index("--dir_plugins") + 1] == "/opt/vep/plugins"

    def test_plugins_dir_absent_when_not_configured(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--dir_plugins" not in cmd

    def test_extra_flags_appended(self):
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["extra_flags"] = ["--pick", "--everything"]
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, settings)
        assert "--pick" in cmd
        assert "--everything" in cmd

    def test_missing_vep_exe_raises_runtime_error(self):
        cfg = _make_paths(vep_path=None)
        with pytest.raises(RuntimeError, match="VEP executable"):
            build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)

    def test_missing_cache_dir_raises_runtime_error(self):
        cfg = _make_paths(cache_dir=None)
        with pytest.raises(RuntimeError, match="VEP cache directory"):
            build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)

    def test_missing_fasta_raises_runtime_error(self):
        cfg = _make_paths(fasta=None)
        with pytest.raises(RuntimeError, match="GRCh38 reference FASTA"):
            build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)

    def test_error_message_contains_setup_instructions(self):
        cfg = _make_paths(vep_path=None)
        with pytest.raises(RuntimeError, match="python cli.py --setup"):
            build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)

    def test_optional_resources_absent_does_not_raise(self):
        """All optional resources absent → command still builds successfully."""
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert cmd[0] == "/opt/vep/vep"

    def test_clinvar_custom_in_command_when_configured(self):
        cfg = _make_paths(clinvar="/data/clinvar.vcf.gz")
        cmd = build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)
        assert "--custom" in cmd
        custom_val = cmd[cmd.index("--custom") + 1]
        assert "clinvar.vcf.gz" in custom_val

    def test_no_custom_in_command_when_clinvar_absent(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert "--custom" not in cmd

    def test_command_is_list_of_strings(self):
        cmd = build_vep_command("in.vcf", "out.vcf", MINIMAL_PATHS, MINIMAL_SETTINGS)
        assert isinstance(cmd, list)
        assert all(isinstance(t, str) for t in cmd)

    def test_all_plugins_present_in_command(self):
        cfg = _make_paths(
            plugins_dir="/opt/plugins",
            loftee_path="/opt/loftee_grch38",
            clinvar="/data/clinvar.vcf.gz",
            cadd_snv="/data/cadd_snv.tsv.gz",
            cadd_indel="/data/cadd_indel.tsv.gz",
            revel="/data/revel.tsv.gz",
            spliceai_snv="/data/spliceai_snv.vcf.gz",
            spliceai_indel="/data/spliceai_indel.vcf.gz",
            alphamissense="/data/am.tsv.gz",
            loftee_human_ancestor_fa="/data/human_ancestor.fa.gz",
            loftee_conservation_file="/data/loftee.sql",
            loftee_gerp_bigwig="/data/gerp.bw",
        )
        cmd = build_vep_command("in.vcf", "out.vcf", cfg, MINIMAL_SETTINGS)
        cmd_str = " ".join(cmd)
        assert "CADD" in cmd_str
        assert "REVEL" in cmd_str
        assert "SpliceAI" in cmd_str
        assert "AlphaMissense" in cmd_str
        assert "ClinVar" in cmd_str
        assert "LoF" in cmd_str


# ---------------------------------------------------------------------------
# run_vep — mocked subprocess
# ---------------------------------------------------------------------------

class TestRunVep:
    def test_success_returns_output_path(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")  # simulate VEP output
        cmd = ["/opt/vep/vep", "--input_file", "in.vcf", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = run_vep(cmd, str(out))

        assert result == str(out)

    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="ERROR: Could not open cache\n",
            )
            with pytest.raises(RuntimeError, match="exited with error code 1"):
                run_vep(cmd, str(out))

    def test_nonzero_exit_error_contains_vep_stderr(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="ERROR: Could not open cache at /bad/path\n",
            )
            with pytest.raises(RuntimeError, match="Could not open cache"):
                run_vep(cmd, str(out))

    def test_nonzero_exit_suggests_setup(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error")
            with pytest.raises(RuntimeError, match="python cli.py --setup"):
                run_vep(cmd, str(out))

    def test_file_not_found_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd = ["/nonexistent/vep", "--output_file", str(out)]

        with patch("subprocess.run", side_effect=FileNotFoundError("not found")):
            with pytest.raises(RuntimeError, match="not found at runtime"):
                run_vep(cmd, str(out))

    def test_timeout_raises_runtime_error(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run", side_effect=sp.TimeoutExpired(cmd, 3600)):
            with pytest.raises(RuntimeError, match="timed out"):
                run_vep(cmd, str(out))

    def test_output_file_not_created_raises_runtime_error(self, tmp_path):
        out = tmp_path / "out.vcf"
        # VEP exits 0 but does not write output_vcf (file does not exist)
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            with pytest.raises(RuntimeError, match="output file was not produced"):
                run_vep(cmd, str(out))


# ---------------------------------------------------------------------------
# annotate() — integration (both build and run mocked)
# ---------------------------------------------------------------------------

class TestAnnotate:
    def test_annotate_returns_output_path(self, tmp_path):
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = annotate("in.vcf", str(out), MINIMAL_PATHS, MINIMAL_SETTINGS)

        assert result == str(out)

    def test_annotate_creates_output_directory(self, tmp_path):
        out = tmp_path / "subdir" / "nested" / "out.vcf"
        out.parent.mkdir(parents=True)
        out.write_text("##fileformat=VCFv4.2\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            annotate("in.vcf", str(out), MINIMAL_PATHS, MINIMAL_SETTINGS)

        assert out.parent.exists()

    def test_annotate_raises_on_missing_vep(self, tmp_path):
        out = tmp_path / "out.vcf"
        cfg = _make_paths(vep_path=None)
        with pytest.raises(RuntimeError, match="VEP executable"):
            annotate("in.vcf", str(out), cfg, MINIMAL_SETTINGS)


# ---------------------------------------------------------------------------
# _resolve_timeout
# ---------------------------------------------------------------------------

class TestResolveTimeout:
    def test_default_returns_6_hours_in_seconds(self):
        assert _resolve_timeout({}) == 6 * 3600

    def test_explicit_hours_converted_to_seconds(self):
        settings = {"vep": {"timeout_hours": 10}}
        assert _resolve_timeout(settings) == 10 * 3600

    def test_zero_returns_none(self):
        settings = {"vep": {"timeout_hours": 0}}
        assert _resolve_timeout(settings) is None

    def test_negative_returns_none(self):
        settings = {"vep": {"timeout_hours": -1}}
        assert _resolve_timeout(settings) is None

    def test_float_hours_accepted(self):
        settings = {"vep": {"timeout_hours": 1.5}}
        assert _resolve_timeout(settings) == int(1.5 * 3600)

    def test_invalid_string_falls_back_to_6h(self):
        settings = {"vep": {"timeout_hours": "notanumber"}}
        assert _resolve_timeout(settings) == 6 * 3600

    def test_none_value_falls_back_to_6h(self):
        settings = {"vep": {"timeout_hours": None}}
        assert _resolve_timeout(settings) == 6 * 3600

    def test_missing_vep_key_uses_default(self):
        settings = {"other": {}}
        assert _resolve_timeout(settings) == 6 * 3600


# ---------------------------------------------------------------------------
# _output_is_complete
# ---------------------------------------------------------------------------

class TestOutputIsComplete:
    def test_returns_false_when_file_does_not_exist(self, tmp_path):
        assert _output_is_complete(str(tmp_path / "missing.vcf")) is False

    def test_returns_false_when_file_empty(self, tmp_path):
        f = tmp_path / "empty.vcf"
        f.write_text("")
        assert _output_is_complete(str(f)) is False

    def test_returns_false_when_no_csq_header(self, tmp_path):
        f = tmp_path / "no_csq.vcf"
        f.write_text("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        assert _output_is_complete(str(f)) is False

    def test_returns_false_when_csq_header_but_no_records(self, tmp_path):
        """CSQ header alone is not enough — VEP may crash before writing any records."""
        f = tmp_path / "headers_only.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotation\">\n"
        )
        assert _output_is_complete(str(f)) is False

    def test_returns_true_when_csq_header_and_data_record_present(self, tmp_path):
        f = tmp_path / "annotated.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotation\">\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=foo\n"
        )
        assert _output_is_complete(str(f)) is True

    def test_returns_true_for_gzipped_csq_file_with_records(self, tmp_path):
        import gzip as gz
        f = tmp_path / "annotated.vcf.gz"
        content = (
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotation\">\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=foo\n"
        )
        with gz.open(str(f), "wt") as fh:
            fh.write(content)
        assert _output_is_complete(str(f)) is True

    def test_returns_false_for_gzipped_csq_header_only(self, tmp_path):
        """Gzipped file with only headers (no records) must return False."""
        import gzip as gz
        f = tmp_path / "headers_only.vcf.gz"
        with gz.open(str(f), "wt") as fh:
            fh.write(
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotation\">\n"
            )
        assert _output_is_complete(str(f)) is False

    def test_returns_false_for_gzipped_without_csq(self, tmp_path):
        import gzip as gz
        f = tmp_path / "no_csq.vcf.gz"
        with gz.open(str(f), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n")
        assert _output_is_complete(str(f)) is False

    def test_stops_at_data_line_even_if_csq_appears_later(self, tmp_path):
        """CSQ after data line should not be found (stops scanning at first non-## line)."""
        f = tmp_path / "csq_late.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"late\">\n"
        )
        assert _output_is_complete(str(f)) is False

    def test_returns_false_on_unreadable_file(self, tmp_path):
        """Should not raise even if the file cannot be read."""
        result = _output_is_complete("/proc/self/mem")
        assert result is False


# ---------------------------------------------------------------------------
# _count_vcf_records
# ---------------------------------------------------------------------------

class TestCountVcfRecords:
    def test_counts_data_lines(self, tmp_path):
        f = tmp_path / "test.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\n"
            "chr1\t200\t.\tG\tC\t.\tPASS\t.\n"
        )
        assert _count_vcf_records(str(f)) == 2

    def test_empty_file_returns_zero(self, tmp_path):
        f = tmp_path / "empty.vcf"
        f.write_text("")
        assert _count_vcf_records(str(f)) == 0

    def test_missing_file_returns_none(self, tmp_path):
        assert _count_vcf_records(str(tmp_path / "missing.vcf")) is None


# ---------------------------------------------------------------------------
# run_vep — timeout configuration
# ---------------------------------------------------------------------------

class TestRunVepTimeout:
    def test_timeout_from_settings_passed_to_subprocess(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]
        settings = {"vep": {"timeout_hours": 3}}

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            out.write_text("##fileformat=VCFv4.2\n")
            run_vep(cmd, str(out), settings)
            _, kwargs = mock_run.call_args
            assert kwargs.get("timeout") == 3 * 3600

    def test_timeout_disabled_when_zero(self, tmp_path):
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]
        settings = {"vep": {"timeout_hours": 0}}

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            out.write_text("##fileformat=VCFv4.2\n")
            run_vep(cmd, str(out), settings)
            _, kwargs = mock_run.call_args
            assert kwargs.get("timeout") is None

    def test_timeout_error_message_shows_current_timeout(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]
        settings = {"vep": {"timeout_hours": 4, "fork": 4}}

        with patch("subprocess.run", side_effect=sp.TimeoutExpired(cmd, 4 * 3600)):
            with pytest.raises(RuntimeError, match="4.0 hour"):
                run_vep(cmd, str(out), settings)

    def test_timeout_error_message_shows_increase_instruction(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]
        settings = {"vep": {"timeout_hours": 2, "fork": 4}}

        with patch("subprocess.run", side_effect=sp.TimeoutExpired(cmd, 2 * 3600)):
            with pytest.raises(RuntimeError, match="timeout_hours"):
                run_vep(cmd, str(out), settings)

    def test_timeout_error_message_shows_bcftools_split_command(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]
        settings = {"vep": {"timeout_hours": 1, "fork": 4}}

        with patch("subprocess.run", side_effect=sp.TimeoutExpired(cmd, 3600)):
            with pytest.raises(RuntimeError, match="bcftools"):
                run_vep(cmd, str(out), settings)

    def test_default_timeout_applied_when_settings_empty(self, tmp_path):
        import subprocess as sp
        out = tmp_path / "out.vcf"
        cmd = ["/opt/vep/vep", "--output_file", str(out)]

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            out.write_text("##fileformat=VCFv4.2\n")
            run_vep(cmd, str(out), {})
            _, kwargs = mock_run.call_args
            assert kwargs.get("timeout") == 6 * 3600


# ---------------------------------------------------------------------------
# annotate() — skip-if-complete behaviour
# ---------------------------------------------------------------------------

class TestAnnotateSkipIfComplete:
    def _make_complete_vcf(self, path):
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotation\">\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=foo\n"
        )

    def test_skips_vep_when_output_complete_and_skip_enabled(self, tmp_path):
        out = tmp_path / "out.vcf"
        self._make_complete_vcf(out)
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["skip_if_annotated_exists"] = True

        with patch("subprocess.run") as mock_run:
            result = annotate("in.vcf", str(out), MINIMAL_PATHS, settings)
            mock_run.assert_not_called()

        assert result == str(out)

    def test_runs_vep_when_output_complete_but_skip_disabled(self, tmp_path):
        out = tmp_path / "out.vcf"
        self._make_complete_vcf(out)
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["skip_if_annotated_exists"] = False

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            annotate("in.vcf", str(out), MINIMAL_PATHS, settings)
            mock_run.assert_called_once()

    def test_runs_vep_when_output_missing(self, tmp_path):
        out = tmp_path / "out.vcf"
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["skip_if_annotated_exists"] = True

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            out.write_text("##fileformat=VCFv4.2\n")  # VEP creates file
            annotate("in.vcf", str(out), MINIMAL_PATHS, settings)
            mock_run.assert_called_once()

    def test_runs_vep_when_output_partial_no_csq(self, tmp_path):
        """Partial output (no CSQ header) must not be mistaken for a complete run."""
        out = tmp_path / "out.vcf"
        out.write_text("##fileformat=VCFv4.2\n")  # exists but no CSQ line
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["skip_if_annotated_exists"] = True

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            annotate("in.vcf", str(out), MINIMAL_PATHS, settings)
            mock_run.assert_called_once()

    def test_skip_returns_existing_path(self, tmp_path):
        out = tmp_path / "out.vcf"
        self._make_complete_vcf(out)
        settings = dict(MINIMAL_SETTINGS)
        settings["vep"] = dict(MINIMAL_SETTINGS["vep"])
        settings["vep"]["skip_if_annotated_exists"] = True

        with patch("subprocess.run"):
            result = annotate("in.vcf", str(out), MINIMAL_PATHS, settings)

        assert result == str(out)

    def test_default_skip_setting_is_true(self, tmp_path):
        """When skip_if_annotated_exists is not set, default is True (skip)."""
        out = tmp_path / "out.vcf"
        self._make_complete_vcf(out)
        # MINIMAL_SETTINGS has no skip_if_annotated_exists key
        settings = {
            "vep": {
                "fork": 4,
                "distance": 5000,
                "extra_flags": [],
            }
        }

        with patch("subprocess.run") as mock_run:
            result = annotate("in.vcf", str(out), MINIMAL_PATHS, settings)
            mock_run.assert_not_called()

        assert result == str(out)
