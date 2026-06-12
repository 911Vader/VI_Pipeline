"""
validator.py — Main orchestration for the validation run.

Workflow:
  1. Load truth set TSV → List[TruthEntry]
  2. Parse annotated VCF  → List[Variant]
  3. Apply filter_engine  → sets excluded / deprioritized flags
  4. Apply classifier     → sets acmg_evidence_codes + acmg_classification
  5. Apply reviewer_flags → sets reviewer_flags + classification_confidence
  6. Build truth index    → Dict[(chrom, pos, ref, alt) → TruthEntry]
  7. ConcordanceEngine    → List[VariantResult] + ValidationMetrics
  8. ToolConcordanceEngine→ ToolConcordanceReport
  9. Build discordance log
 10. Write all output files

The validator classifies ALL variants from the VCF (not just truth set entries)
so that tool concordance statistics are computed over the full dataset.

Usage (API):
    from validation.validator import run_validation
    run_validation(
        truth_set_path="test_data/truth_set.tsv",
        vcf_path="runs/pt1/pt1.norm.vep.vcf",
        output_dir="validation_out/",
        settings=settings_dict,   # from config/settings.yaml
    )
"""

import datetime
import logging
from pathlib import Path
from typing import List, Optional

import yaml

from analysis.parser import parse_vep_output
from analysis.filter_engine import apply_filters
from analysis.acmg.classifier import classify_all
from analysis.reviewer_flags import annotate_reviewer_flags

from validation.truth_set import load_truth_set, build_truth_index
from validation.concordance import ConcordanceEngine
from validation.tool_concordance import compute_tool_concordance
from validation.discordance_log import build_discordance_entries
from validation.reporter import (
    write_pdf_summary,
    write_discordance_xlsx,
    write_tool_concordance_xlsx,
    write_raw_results,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_validation(
    truth_set_path: str,
    vcf_path: str,
    output_dir: str,
    settings: Optional[dict] = None,
    gene_list: Optional[List[str]] = None,
    sample_id: str = "VALIDATION_RUN",
) -> dict:
    """
    Execute a full validation run and write all output artefacts.

    Parameters
    ----------
    truth_set_path: Path to the curated truth set TSV.
    vcf_path:       Path to a VEP-annotated VCF (may be pre-existing).
    output_dir:     Directory for all output files.
    settings:       Parsed settings.yaml dict. If None, loaded from
                    config/settings.yaml.
    gene_list:      Optional gene list for filter_engine.
    sample_id:      Run identifier for metadata.

    Returns
    -------
    dict with keys: pdf, discordance_xlsx, tool_concordance_xlsx,
                    raw_tsv, raw_json, verdict, concordance_rate,
                    n_critical_failures.
    """
    run_start = datetime.datetime.utcnow()
    out_dir   = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    settings = _load_settings(settings)

    logger.info("=== Validation run started: %s ===", run_start.isoformat())
    logger.info("Truth set: %s", truth_set_path)
    logger.info("VCF:       %s", vcf_path)
    logger.info("Output:    %s", out_dir)

    # 1. Load truth set
    truth_entries = load_truth_set(truth_set_path)
    truth_index   = build_truth_index(truth_entries)
    logger.info("Truth set: %d entries", len(truth_entries))

    # 2. Parse VCF
    logger.info("Parsing VCF …")
    variants = parse_vep_output(vcf_path)
    logger.info("VCF parsed: %d raw variants", len(variants))

    # 3. Filter
    logger.info("Applying filters …")
    apply_filters(variants, settings, gene_list=gene_list)

    # 4. Classify
    logger.info("Classifying …")
    classify_all(variants, settings)

    # 5. Reviewer flags
    logger.info("Annotating reviewer flags …")
    annotate_reviewer_flags(variants, settings, gene_list=gene_list)

    # 6. Concordance
    logger.info("Computing concordance …")
    engine  = ConcordanceEngine(truth_index, variants)
    results, metrics = engine.run()

    # 7. Tool concordance (over all VCF variants — larger sample)
    logger.info("Computing tool concordance …")
    tool_report = compute_tool_concordance(variants, settings)

    # 8. Discordance log
    discord_entries = build_discordance_entries(results)
    logger.info(
        "%d discordant variants logged (%d critical failures)",
        len(discord_entries), metrics.n_critical_failure,
    )

    # 9. Write outputs
    metadata = {
        "sample_id":      sample_id,
        "truth_set_path": truth_set_path,
        "vcf_path":       vcf_path,
        "run_date":       run_start.strftime("%Y-%m-%d %H:%M UTC"),
    }

    pdf_path  = str(out_dir / "validation_summary.pdf")
    xl1_path  = str(out_dir / "discordance_log.xlsx")
    xl2_path  = str(out_dir / "tool_concordance.xlsx")

    write_pdf_summary(pdf_path, metrics, tool_report, metadata)
    write_discordance_xlsx(xl1_path, discord_entries)
    write_tool_concordance_xlsx(xl2_path, tool_report)
    tsv_path, json_path = write_raw_results(str(out_dir), results)

    verdict, _ = metrics.pass_fail_verdict()
    rate = metrics.concordance_rate

    logger.info(
        "=== Validation complete: %s | Concordance %.1f%% | %d critical failures ===",
        verdict, (rate or 0) * 100, metrics.n_critical_failure,
    )

    return {
        "pdf":                    pdf_path,
        "discordance_xlsx":       xl1_path,
        "tool_concordance_xlsx":  xl2_path,
        "raw_tsv":                tsv_path,
        "raw_json":               json_path,
        "verdict":                verdict,
        "concordance_rate":       rate,
        "n_critical_failures":    metrics.n_critical_failure,
        "metrics":                metrics,
        "tool_concordance":       tool_report,
        "results":                results,
    }


# ---------------------------------------------------------------------------
# Settings loader
# ---------------------------------------------------------------------------

def _load_settings(settings: Optional[dict]) -> dict:
    if settings is not None:
        return settings
    settings_path = Path(__file__).parent.parent / "config" / "settings.yaml"
    if settings_path.exists():
        with settings_path.open() as fh:
            return yaml.safe_load(fh) or {}
    logger.warning("config/settings.yaml not found — using empty settings dict")
    return {}
