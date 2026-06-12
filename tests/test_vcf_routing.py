"""
tests/test_vcf_routing.py — Tests for the automatic annotated/unannotated VCF
routing logic introduced in v1.1.

Tests the decision point in cli.py that chooses between:
  - Annotated path:   CSQ header present → skip VEP → parse directly
  - Unannotated path: No CSQ header      → run VEP  → parse annotated output

The tests mock at module boundaries:
  - analysis.vcf_validator.validate_vcf  (filesystem check)
  - analysis.vcf_validator.has_csq_annotation  (header peek)
  - analysis.annotator.annotate          (subprocess boundary)
  - analysis.parser.parse_vep_output     (reads real VCF data)
  - downstream steps (filter, classify, rank, reporting)

No real VCF files, VEP executables, or filesystem I/O are needed.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch, call
import pytest

from analysis.vcf_validator import has_csq_annotation


# ---------------------------------------------------------------------------
# has_csq_annotation — already tested in test_vcf_validator.py
# These tests focus specifically on the routing logic triggered by its return value.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Helpers — minimal Variant-like objects for mocking parse output
# ---------------------------------------------------------------------------

def _make_mock_variant(excluded=False, acmg_classification="VUS"):
    v = MagicMock()
    v.excluded = excluded
    v.deprioritized = False
    v.deprioritization_reason = None
    v.acmg_classification = acmg_classification
    v.acmg_evidence_codes = []
    v.gene_symbol = "GENE1"
    v.in_gene_list = False
    v.gene_relevance = ""
    v.reviewer_flags = []
    # Numeric fields must be None (not MagicMock) so reviewer_flags comparisons work
    v.gnomad_af = None
    v.cadd_phred = None
    v.revel_score = None
    v.spliceai_max_score = None
    v.loftee_lof = None
    v.loftee_filter = None
    v.clinvar_significance = None
    v.clinvar_review_status = None
    v.depth = None
    v.genotype_quality = None
    v.vaf = None
    v.filter_status = "PASS"
    return v


# ---------------------------------------------------------------------------
# Shared mock config that satisfies assert_setup_complete()
# ---------------------------------------------------------------------------

_MOCK_PATHS_CONFIG = {
    "setup_completed": True,
    "tools": {
        "vep": {"path": "/opt/vep/vep", "status": "found"},
    },
    "resources": {
        "vep_cache_dir":   {"path": "/data/vep_cache", "status": "found"},
        "reference_fasta": {"path": "/data/GRCh38.fa",  "status": "found"},
        "vep_plugins_dir": {"path": None, "status": "not_provided"},
        "clinvar_vcf":     {"path": None, "status": "not_provided"},
        "cadd_snv_file":   {"path": None, "status": "not_provided"},
        "cadd_indel_file": {"path": None, "status": "not_provided"},
        "revel_file":      {"path": None, "status": "not_provided"},
        "spliceai_snv_file":   {"path": None, "status": "not_provided"},
        "spliceai_indel_file": {"path": None, "status": "not_provided"},
        "alphamissense_file": {"path": None, "status": "not_provided"},
    },
}

_MOCK_SETTINGS = {
    "vep": {"fork": 4, "distance": 5000, "extra_flags": []},
    "reporting": {"max_candidates_in_report": 20},
    "filtering": {"require_pass": True, "max_af": 0.01},
    "acmg": {"ba1_af_threshold": 0.05, "insilico_thresholds": {}},
    "reviewer_flags": {},
}


def _run_cli_main(tmp_path, input_vcf, extra_args=None):
    """
    Import and call cli.main() with a controlled argv.
    Returns the integer exit code.
    """
    import importlib
    import cli
    importlib.reload(cli)  # fresh state each test

    args = [
        "vi_pipeline",
        "--input", input_vcf,
        "--output", str(tmp_path / "results"),
        "--sample-id", "STEST",
    ]
    if extra_args:
        args += extra_args

    with patch.object(sys, "argv", args):
        return cli.main()


# ---------------------------------------------------------------------------
# Annotated VCF path — CSQ present, VEP should NOT be called
# ---------------------------------------------------------------------------

class TestAnnotatedVcfPath:

    def _full_patch(self, tmp_path, has_csq=True):
        """Context manager stack that patches the full pipeline for one CLI run."""
        variants = [_make_mock_variant()]
        return (
            patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG),
            patch("yaml.safe_load", return_value=_MOCK_SETTINGS),
            patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"),
            patch("analysis.vcf_validator.has_csq_annotation", return_value=has_csq),
            patch("analysis.annotator.annotate", return_value="/output/STEST.vep.vcf"),
            patch("analysis.parser.parse_vep_output", return_value=variants),
            patch("analysis.filter_engine.apply_filters", return_value=variants),
            patch("analysis.acmg.classifier.classify_all"),
            patch("analysis.ranker.rank_variants", return_value=variants),
            patch("reporting.json_export.write_json", return_value=str(tmp_path / "out.json")),
            patch("reporting.spreadsheet.write_spreadsheet", return_value=str(tmp_path / "out.xlsx")),
            patch("reporting.pdf_report.write_pdf", return_value=str(tmp_path / "out.pdf")),
        )

    def test_annotated_vcf_exits_zero(self, tmp_path):
        patches = self._full_patch(tmp_path, has_csq=True)
        with patches[0], patches[1], patches[2], patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11]:
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")
        assert rc == 0

    def test_annotated_vcf_skips_annotation(self, tmp_path):
        patches = self._full_patch(tmp_path, has_csq=True)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_annotate, patches[5], patches[6], \
             patches[7], patches[8], patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        mock_annotate.assert_not_called()

    def test_annotated_vcf_parses_input_directly(self, tmp_path):
        patches = self._full_patch(tmp_path, has_csq=True)
        with patches[0], patches[1], patches[2], patches[3], patches[4], \
             patches[5] as mock_parse, patches[6], patches[7], patches[8], \
             patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        # Should parse the original input, not an intermediate annotated file
        mock_parse.assert_called_once_with("/input/sample.vcf")


# ---------------------------------------------------------------------------
# Unannotated VCF path — no CSQ, VEP MUST be called
# ---------------------------------------------------------------------------

class TestUnannotatedVcfPath:

    def _full_patch(self, tmp_path, annotated_out=None):
        variants = [_make_mock_variant()]
        annotated_out = annotated_out or str(tmp_path / "results" / "STEST.vep.vcf")
        return (
            patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG),
            patch("yaml.safe_load", return_value=_MOCK_SETTINGS),
            patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"),
            patch("analysis.vcf_validator.has_csq_annotation", return_value=False),
            patch("analysis.annotator.annotate", return_value=annotated_out),
            patch("analysis.parser.parse_vep_output", return_value=variants),
            patch("analysis.filter_engine.apply_filters", return_value=variants),
            patch("analysis.acmg.classifier.classify_all"),
            patch("analysis.ranker.rank_variants", return_value=variants),
            patch("reporting.json_export.write_json", return_value=str(tmp_path / "out.json")),
            patch("reporting.spreadsheet.write_spreadsheet", return_value=str(tmp_path / "out.xlsx")),
            patch("reporting.pdf_report.write_pdf", return_value=str(tmp_path / "out.pdf")),
        )

    def test_unannotated_vcf_exits_zero(self, tmp_path):
        patches = self._full_patch(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11]:
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")
        assert rc == 0

    def test_unannotated_vcf_calls_annotate(self, tmp_path):
        patches = self._full_patch(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_annotate, patches[5], patches[6], \
             patches[7], patches[8], patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        mock_annotate.assert_called_once()

    def test_unannotated_annotate_called_with_correct_input(self, tmp_path):
        patches = self._full_patch(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_annotate, patches[5], patches[6], \
             patches[7], patches[8], patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        # First positional arg to annotate() must be the original input VCF
        call_args = mock_annotate.call_args
        assert call_args[0][0] == "/input/sample.vcf"

    def test_unannotated_output_vcf_uses_sample_id(self, tmp_path):
        """Annotated VCF output path must include sample_id and .vep.vcf suffix."""
        patches = self._full_patch(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_annotate, patches[5], patches[6], \
             patches[7], patches[8], patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        output_arg = mock_annotate.call_args[0][1]
        assert "STEST" in output_arg
        assert output_arg.endswith(".vep.vcf")

    def test_unannotated_parse_uses_annotated_output(self, tmp_path):
        """Parser must receive the VEP output path, not the original input."""
        annotated_out = str(tmp_path / "results" / "STEST.vep.vcf")
        patches = self._full_patch(tmp_path, annotated_out=annotated_out)
        with patches[0], patches[1], patches[2], patches[3], patches[4], \
             patches[5] as mock_parse, patches[6], patches[7], patches[8], \
             patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        mock_parse.assert_called_once_with(annotated_out)

    def test_unannotated_output_does_not_overwrite_input(self, tmp_path):
        """The annotated VCF output path must never equal the input path."""
        patches = self._full_patch(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_annotate, patches[5], patches[6], \
             patches[7], patches[8], patches[9], patches[10], patches[11]:
            _run_cli_main(tmp_path, "/input/sample.vcf")
        output_arg = mock_annotate.call_args[0][1]
        assert output_arg != "/input/sample.vcf"


# ---------------------------------------------------------------------------
# Annotation failure path
# ---------------------------------------------------------------------------

class TestAnnotationFailurePath:

    def test_annotation_failure_returns_nonzero(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=False), \
             patch("analysis.annotator.annotate",
                   side_effect=RuntimeError("VEP exited with error code 1")), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")
        assert rc != 0

    def test_annotation_failure_does_not_call_parser(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=False), \
             patch("analysis.annotator.annotate",
                   side_effect=RuntimeError("VEP failed")), \
             patch("analysis.parser.parse_vep_output") as mock_parse, \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            _run_cli_main(tmp_path, "/input/sample.vcf")
        mock_parse.assert_not_called()

    def test_validation_failure_returns_nonzero(self, tmp_path):
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf",
                   side_effect=ValueError("Input file not found")):
            rc = _run_cli_main(tmp_path, "/nonexistent/sample.vcf")
        assert rc != 0


# ---------------------------------------------------------------------------
# Missing VEP / missing resources — annotate() raises, CLI handles it
# ---------------------------------------------------------------------------

class TestMissingVepResources:

    def test_missing_vep_exe_returns_nonzero(self, tmp_path):
        variants = [_make_mock_variant()]
        no_vep_config = {
            **_MOCK_PATHS_CONFIG,
            "tools": {"vep": {"path": None, "status": "missing"}},
        }
        with patch("setup.config_writer.assert_setup_complete", return_value=no_vep_config), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=False), \
             patch("analysis.annotator.annotate",
                   side_effect=RuntimeError("VEP executable is required but not available")), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")
        assert rc != 0

    def test_missing_fasta_returns_nonzero(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=False), \
             patch("analysis.annotator.annotate",
                   side_effect=RuntimeError("GRCh38 reference FASTA is required")), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")
        assert rc != 0


# ---------------------------------------------------------------------------
# Preservation of existing annotated VCF behavior
# ---------------------------------------------------------------------------

class TestAnnotatedPathPreservation:
    """
    Regression tests: verify that the annotated-VCF path (the v1.0 behavior)
    continues to work exactly as before after the v1.1 routing was added.
    """

    def test_annotated_path_calls_parse_filter_classify_rank_report(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/ann.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=True), \
             patch("analysis.annotator.annotate") as mock_annotate, \
             patch("analysis.parser.parse_vep_output", return_value=variants) as mock_parse, \
             patch("analysis.filter_engine.apply_filters", return_value=variants) as mock_filter, \
             patch("analysis.acmg.classifier.classify_all") as mock_classify, \
             patch("analysis.ranker.rank_variants", return_value=variants) as mock_rank, \
             patch("reporting.json_export.write_json", return_value="out.json") as mock_json, \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/ann.vcf")

        assert rc == 0
        mock_annotate.assert_not_called()
        mock_parse.assert_called_once_with("/input/ann.vcf")
        mock_filter.assert_called_once()
        mock_classify.assert_called_once()
        mock_rank.assert_called_once()
        mock_json.assert_called_once()

    def test_annotated_path_run_log_records_detection_step(self, tmp_path):
        """The run log must record the vcf_detection step even for annotated VCFs."""
        variants = [_make_mock_variant()]
        captured_steps = []

        def fake_add_step(name, status="completed", detail=None):
            captured_steps.append(name)

        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf", return_value="/input/ann.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=True), \
             patch("analysis.annotator.annotate"), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"), \
             patch("reporting.run_log.RunLog.add_step", side_effect=fake_add_step), \
             patch("reporting.run_log.RunLog.finalize", return_value="run.log.json"):
            _run_cli_main(tmp_path, "/input/ann.vcf")

        assert "vcf_detection" in captured_steps


# ---------------------------------------------------------------------------
# FASTQ routing through CLI — v1.2 additions
# ---------------------------------------------------------------------------

class TestFastqCliRouting:
    """
    Tests for the FASTQ input branch added in v1.2.
    The FASTQ path:
      detect_input_type → FASTQ → run_fastq_to_vcf → annotate → parse → rest
    """

    def _full_patch_fastq(self, tmp_path, raw_vcf=None, annotated_vcf=None):
        variants = [_make_mock_variant()]
        raw_vcf = raw_vcf or str(tmp_path / "results" / "STEST.raw.vcf")
        annotated_vcf = annotated_vcf or str(tmp_path / "results" / "STEST.vep.vcf")
        return (
            patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG),
            patch("yaml.safe_load", return_value=_MOCK_SETTINGS),
            patch("analysis.upstream.fastq_handler.run_fastq_to_vcf", return_value=raw_vcf),
            patch("analysis.annotator.annotate", return_value=annotated_vcf),
            patch("analysis.parser.parse_vep_output", return_value=variants),
            patch("analysis.filter_engine.apply_filters", return_value=variants),
            patch("analysis.acmg.classifier.classify_all"),
            patch("analysis.ranker.rank_variants", return_value=variants),
            patch("reporting.json_export.write_json", return_value=str(tmp_path / "out.json")),
            patch("reporting.spreadsheet.write_spreadsheet", return_value=str(tmp_path / "out.xlsx")),
            patch("reporting.pdf_report.write_pdf", return_value=str(tmp_path / "out.pdf")),
        )

    def test_fastq_input_exits_zero(self, tmp_path):
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            rc = _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        assert rc == 0

    def test_fastq_calls_run_fastq_to_vcf(self, tmp_path):
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2] as mock_fq, patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        mock_fq.assert_called_once()

    def test_fastq_passes_r1_correctly(self, tmp_path):
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2] as mock_fq, patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        assert mock_fq.call_args.kwargs["fastq_r1"] == "/input/sample.fastq.gz"

    def test_fastq_single_end_r2_is_none(self, tmp_path):
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2] as mock_fq, patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        assert mock_fq.call_args.kwargs["fastq_r2"] is None

    def test_fastq_always_calls_annotate(self, tmp_path):
        """FASTQ output is always unannotated — annotate() must always be called."""
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2], patches[3] as mock_ann, patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        mock_ann.assert_called_once()

    def test_fastq_parse_receives_annotated_vcf(self, tmp_path):
        annotated = str(tmp_path / "results" / "STEST.vep.vcf")
        patches = self._full_patch_fastq(tmp_path, annotated_vcf=annotated)
        with patches[0], patches[1], patches[2], patches[3], patches[4] as mock_parse, \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        mock_parse.assert_called_once_with(annotated)

    def test_fastq_handler_failure_returns_nonzero(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.upstream.fastq_handler.run_fastq_to_vcf",
                   side_effect=RuntimeError("BWA alignment failed")), \
             patch("analysis.annotator.annotate", return_value="out.vep.vcf"), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        assert rc != 0

    def test_fastq_annotation_failure_returns_nonzero(self, tmp_path):
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.upstream.fastq_handler.run_fastq_to_vcf",
                   return_value="/tmp/raw.vcf"), \
             patch("analysis.annotator.annotate",
                   side_effect=RuntimeError("VEP exited with error code 1")), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json", return_value="out.json"), \
             patch("reporting.spreadsheet.write_spreadsheet", return_value="out.xlsx"), \
             patch("reporting.pdf_report.write_pdf", return_value="out.pdf"):
            rc = _run_cli_main(tmp_path, "/input/sample.fastq.gz")
        assert rc != 0

    def test_fq_gz_extension_also_routes_to_fastq(self, tmp_path):
        patches = self._full_patch_fastq(tmp_path)
        with patches[0], patches[1], patches[2] as mock_fq, patches[3], patches[4], \
             patches[5], patches[6], patches[7], patches[8], patches[9], patches[10]:
            _run_cli_main(tmp_path, "/input/sample.fq.gz")
        mock_fq.assert_called_once()


class TestBamCramAndUnknownStillBlocked:
    """Regression: BAM/CRAM and unknown types still fail after v1.2 changes."""

    def test_bam_returns_nonzero(self, tmp_path):
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS):
            rc = _run_cli_main(tmp_path, "/input/sample.bam")
        assert rc != 0

    def test_unknown_extension_returns_nonzero(self, tmp_path):
        with patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS):
            rc = _run_cli_main(tmp_path, "/input/sample.xyz")
        assert rc != 0


# ---------------------------------------------------------------------------
# Normalization routing — v1.3 additions
# ---------------------------------------------------------------------------

class TestNormalizationRouting:
    """
    Tests for the multi-allelic normalization step added in v1.3.

    The raw VCF path (unannotated, no CSQ) now checks has_multiallelic_records():
      - True  → normalize_vcf() called → normalized path passed to annotate()
      - False → normalize_vcf() NOT called → original path passed to annotate()

    The annotated VCF path (CSQ present) is never normalized.
    """

    def _base_patches(self, tmp_path):
        """Return the standard patch set for a raw unannotated VCF run."""
        variants = [_make_mock_variant()]
        annotated_out = str(tmp_path / "results" / "STEST.vep.vcf")
        return (
            patch("setup.config_writer.assert_setup_complete", return_value=_MOCK_PATHS_CONFIG),
            patch("yaml.safe_load", return_value=_MOCK_SETTINGS),
            patch("analysis.vcf_validator.validate_vcf", return_value="/input/sample.vcf"),
            patch("analysis.vcf_validator.has_csq_annotation", return_value=False),
            patch("analysis.vcf_validator.has_multiallelic_records"),
            patch("analysis.vcf_normalizer.normalize_vcf"),
            patch("analysis.annotator.annotate", return_value=annotated_out),
            patch("analysis.parser.parse_vep_output", return_value=variants),
            patch("analysis.filter_engine.apply_filters", return_value=variants),
            patch("analysis.acmg.classifier.classify_all"),
            patch("analysis.ranker.rank_variants", return_value=variants),
            patch("reporting.json_export.write_json",
                  return_value=str(tmp_path / "out.json")),
            patch("reporting.spreadsheet.write_spreadsheet",
                  return_value=str(tmp_path / "out.xlsx")),
            patch("reporting.pdf_report.write_pdf",
                  return_value=str(tmp_path / "out.pdf")),
        )

    def test_multiallelic_detected_calls_normalize(self, tmp_path):
        """When has_multiallelic_records returns True, normalize_vcf must be called."""
        norm_out = str(tmp_path / "results" / "STEST.norm.vcf")
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.return_value = norm_out
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")

        assert rc == 0
        mock_norm.assert_called_once()

    def test_no_multiallelic_skips_normalize(self, tmp_path):
        """When has_multiallelic_records returns False, normalize_vcf must NOT be called."""
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = False
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")

        assert rc == 0
        mock_norm.assert_not_called()

    def test_normalized_path_passed_to_annotate(self, tmp_path):
        """annotate() must receive the normalized VCF path, not the original input."""
        norm_out = str(tmp_path / "results" / "STEST.norm.vcf")
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6] as mock_annotate, \
             patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.return_value = norm_out
            _run_cli_main(tmp_path, "/input/sample.vcf")

        annotate_input = mock_annotate.call_args[0][0]
        assert annotate_input == norm_out

    def test_original_path_passed_to_annotate_when_no_normalization(self, tmp_path):
        """When normalization is skipped, annotate() must still receive the original path."""
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5], \
             patches[6] as mock_annotate, \
             patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = False
            _run_cli_main(tmp_path, "/input/sample.vcf")

        annotate_input = mock_annotate.call_args[0][0]
        assert annotate_input == "/input/sample.vcf"

    def test_normalized_output_filename_contains_sample_id(self, tmp_path):
        """The .norm.vcf output path must embed the sample ID."""
        norm_out = str(tmp_path / "results" / "STEST.norm.vcf")
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.return_value = norm_out
            _run_cli_main(tmp_path, "/input/sample.vcf")

        norm_call_args = mock_norm.call_args[0]
        norm_output_path = norm_call_args[1]
        assert "STEST" in norm_output_path
        assert norm_output_path.endswith(".norm.vcf")

    def test_normalized_output_does_not_overwrite_input(self, tmp_path):
        """Normalized VCF path must never equal the original input path."""
        norm_out = str(tmp_path / "results" / "STEST.norm.vcf")
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.return_value = norm_out
            _run_cli_main(tmp_path, "/input/sample.vcf")

        norm_output_path = mock_norm.call_args[0][1]
        assert norm_output_path != "/input/sample.vcf"

    def test_normalization_failure_returns_nonzero(self, tmp_path):
        """If normalize_vcf raises, the CLI must return non-zero."""
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6], patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.side_effect = RuntimeError("bcftools norm failed (exit code 1)")
            rc = _run_cli_main(tmp_path, "/input/sample.vcf")

        assert rc != 0

    def test_normalization_failure_does_not_call_annotate(self, tmp_path):
        """If normalize_vcf raises, annotate() must not be called."""
        patches = self._base_patches(tmp_path)
        with patches[0], patches[1], patches[2], patches[3], \
             patches[4] as mock_multi, \
             patches[5] as mock_norm, \
             patches[6] as mock_annotate, \
             patches[7], patches[8], patches[9], \
             patches[10], patches[11], patches[12], patches[13]:
            mock_multi.return_value = True
            mock_norm.side_effect = RuntimeError("bcftools norm failed")
            _run_cli_main(tmp_path, "/input/sample.vcf")

        mock_annotate.assert_not_called()

    def test_annotated_vcf_path_never_calls_normalize(self, tmp_path):
        """CSQ-annotated VCF must bypass normalization entirely."""
        variants = [_make_mock_variant()]
        with patch("setup.config_writer.assert_setup_complete",
                   return_value=_MOCK_PATHS_CONFIG), \
             patch("yaml.safe_load", return_value=_MOCK_SETTINGS), \
             patch("analysis.vcf_validator.validate_vcf",
                   return_value="/input/annotated.vcf"), \
             patch("analysis.vcf_validator.has_csq_annotation", return_value=True), \
             patch("analysis.vcf_validator.has_multiallelic_records") as mock_multi, \
             patch("analysis.vcf_normalizer.normalize_vcf") as mock_norm, \
             patch("analysis.annotator.annotate"), \
             patch("analysis.parser.parse_vep_output", return_value=variants), \
             patch("analysis.filter_engine.apply_filters", return_value=variants), \
             patch("analysis.acmg.classifier.classify_all"), \
             patch("analysis.ranker.rank_variants", return_value=variants), \
             patch("reporting.json_export.write_json",
                   return_value=str(tmp_path / "out.json")), \
             patch("reporting.spreadsheet.write_spreadsheet",
                   return_value=str(tmp_path / "out.xlsx")), \
             patch("reporting.pdf_report.write_pdf",
                   return_value=str(tmp_path / "out.pdf")):
            _run_cli_main(tmp_path, "/input/annotated.vcf")

        mock_multi.assert_not_called()
        mock_norm.assert_not_called()
