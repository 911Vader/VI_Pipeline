"""
validate.py — CLI for the VI_Pipeline internal validation framework.

Usage:
    python validate.py --truth-set test_data/truth_set.tsv \\
                       --vcf runs/pt1/pt1.norm.vep.vcf \\
                       --output-dir validation_out/

    # With optional gene list and custom sample label:
    python validate.py --truth-set test_data/truth_set.tsv \\
                       --vcf runs/pt1/pt1.norm.vep.vcf \\
                       --output-dir validation_out/ \\
                       --gene-list resources/gene_list.txt \\
                       --sample-id PT1_VALIDATION

Output files written to <output-dir>/:
    validation_summary.pdf    — PDF with metrics and pass/fail verdict
    discordance_log.xlsx      — full discordance log
    tool_concordance.xlsx     — pairwise tool concordance matrix
    raw_results.tsv           — every truth set variant, expected vs actual
    raw_results.json          — same data, JSON format
"""

import argparse
import logging
import sys
from pathlib import Path


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)-8s] %(name)s — %(message)s",
        datefmt="%H:%M:%S",
    )


def _load_gene_list(path: str):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Gene list not found: {path}")
    genes = []
    with p.open() as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.append(line)
    return genes


def main() -> int:
    parser = argparse.ArgumentParser(
        description="VI_Pipeline internal validation — compare pipeline output to curated truth set",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--truth-set", required=True, metavar="TSV",
        help="Curated truth set TSV (see test_data/truth_set.tsv for format)",
    )
    parser.add_argument(
        "--vcf", required=True, metavar="VCF",
        help="VEP-annotated VCF to validate (may be pre-existing from a pipeline run)",
    )
    parser.add_argument(
        "--output-dir", required=True, metavar="DIR",
        help="Directory for all output files (created if absent)",
    )
    parser.add_argument(
        "--gene-list", metavar="FILE",
        help="Optional gene list (one symbol per line) for filter_engine context",
    )
    parser.add_argument(
        "--sample-id", default="VALIDATION_RUN", metavar="ID",
        help="Run identifier shown in the PDF report (default: VALIDATION_RUN)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()
    _configure_logging(args.verbose)

    log = logging.getLogger("validate")

    # Validate inputs exist before importing heavy modules
    for label, path in [("truth set", args.truth_set), ("VCF", args.vcf)]:
        if not Path(path).exists():
            log.error("%s not found: %s", label, path)
            return 1

    gene_list = None
    if args.gene_list:
        try:
            gene_list = _load_gene_list(args.gene_list)
            log.info("Gene list loaded: %d genes from %s", len(gene_list), args.gene_list)
        except FileNotFoundError as exc:
            log.error("%s", exc)
            return 1

    from validation.validator import run_validation

    try:
        result = run_validation(
            truth_set_path = args.truth_set,
            vcf_path       = args.vcf,
            output_dir     = args.output_dir,
            gene_list      = gene_list,
            sample_id      = args.sample_id,
        )
    except Exception as exc:
        log.error("Validation run failed: %s", exc, exc_info=True)
        return 1

    # Print summary to stdout
    rate = result["concordance_rate"]
    rate_str = f"{rate*100:.1f}%" if rate is not None else "N/A"
    m = result["metrics"]

    print()
    print("=" * 60)
    print(f"  VALIDATION COMPLETE — Verdict: {result['verdict']}")
    print("=" * 60)
    print(f"  Truth set variants : {m.n_truth_total}")
    print(f"  Matched in VCF     : {m.n_matched}  (not found: {m.n_not_found})")
    print(f"  Concordant         : {m.n_concordant}  ({rate_str})")
    print(f"  Borderline         : {m.n_borderline}")
    print(f"  Direction errors   : {m.n_direction_error}")
    print(f"  CRITICAL FAILURES  : {m.n_critical_failure}")
    print(f"  Sensitivity (P/LP) : {_pct(m.sensitivity)}")
    print(f"  Specificity (B/LB) : {_pct(m.specificity)}")
    print(f"  Pipeline VUS rate  : {_pct(m.vus_rate)}")
    print()
    print("  Output files:")
    for key in ("pdf", "discordance_xlsx", "tool_concordance_xlsx", "raw_tsv", "raw_json"):
        print(f"    {key:<25} {result[key]}")
    print("=" * 60)
    print()

    # Non-zero exit on FAIL
    return 0 if result["verdict"] != "FAIL" else 2


def _pct(v) -> str:
    return f"{v*100:.1f}%" if v is not None else "N/A"


if __name__ == "__main__":
    sys.exit(main())
