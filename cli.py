"""
cli.py — VI_Pipeline command-line entry point.

Usage examples:

  Setup (run once before first analysis):
    python cli.py --setup
    python cli.py --setup --force

  Analysis (VCF input, version 1):
    python cli.py --input sample.vcf.gz \
                  --gene-list config/gene_list.txt \
                  --output results/ \
                  --sample-id SAMPLE_001 \
                  --case-id CASE_001

  Using a metadata file instead of inline flags:
    python cli.py --input sample.vcf.gz \
                  --gene-list config/gene_list.txt \
                  --output results/ \
                  --metadata metadata.yaml
"""

import argparse
import logging
import sys
from pathlib import Path

import yaml


def _stage(step: int, total: int, label: str) -> None:
    """Emit a structured stage marker parsed by the GUI progress bar."""
    print(f"##STAGE## {step}/{total} {label}", flush=True)


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=level,
        stream=sys.stdout,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vi_pipeline",
        description=(
            "VI_Pipeline — Forensic genomic variant prioritisation pipeline\n"
            "GRCh38 | Single sample | SNVs and small indels | ACMG/AMP classification"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Setup mode ---
    setup_group = parser.add_argument_group("Setup")
    setup_group.add_argument(
        "--setup",
        action="store_true",
        help="Run tool and resource checks, write config/paths.yaml. "
             "Run this once before first analysis.",
    )
    setup_group.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of setup even if paths.yaml already shows success.",
    )

    # --- Analysis mode ---
    analysis_group = parser.add_argument_group("Analysis")
    analysis_group.add_argument(
        "--input", "-i",
        metavar="FILE",
        help="Input file: VCF, VCF.GZ, FASTQ, FASTQ.GZ, or BAM.",
    )
    analysis_group.add_argument(
        "--input2",
        metavar="FILE",
        help="R2 FASTQ file for paired-end input. Only used when --input is a FASTQ.",
    )
    analysis_group.add_argument(
        "--gene-list", "-g",
        metavar="FILE",
        default="config/gene_list.txt",
        help="Gene list file (one HGNC symbol per line). Default: config/gene_list.txt",
    )
    analysis_group.add_argument(
        "--output", "-o",
        metavar="DIR",
        default="results",
        help="Output directory. Created if it does not exist. Default: results/",
    )

    # --- Metadata ---
    meta_group = parser.add_argument_group("Metadata")
    meta_group.add_argument("--sample-id", metavar="ID", help="Sample identifier.")
    meta_group.add_argument("--case-id", metavar="ID", help="Case identifier (optional).")
    meta_group.add_argument(
        "--metadata",
        metavar="FILE",
        help="YAML file with sample_id, case_id, phenotype_terms. "
             "Overrides --sample-id / --case-id if both supplied.",
    )

    # --- General ---
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging.")
    parser.add_argument("--version", action="version", version="VI_Pipeline 1.0.0")

    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()
    _configure_logging(args.verbose)
    log = logging.getLogger("cli")

    # --- Setup mode ---
    if args.setup:
        from setup.setup_runner import run_setup
        success = run_setup(force=args.force)
        return 0 if success else 1

    # --- Analysis mode ---
    if not args.input:
        parser.print_help()
        print("\nError: --input is required for analysis. Run --setup first if you haven't.")
        return 1

    # Guard: require setup to have completed
    try:
        from setup.config_writer import assert_setup_complete
        paths_config = assert_setup_complete()
    except RuntimeError as exc:
        print(f"\nSetup required:\n  {exc}")
        return 1

    # --- Load settings ---
    settings_path = Path("config/settings.yaml")
    settings: dict = {}
    if settings_path.exists():
        with open(settings_path) as fh:
            settings = yaml.safe_load(fh) or {}

    # --- Resolve metadata ---
    metadata: dict = {
        "sample_id": args.sample_id or "UNKNOWN",
        "case_id":   args.case_id or "",
        "input_vcf": args.input,
        "gene_list_file": args.gene_list,
    }
    if args.metadata:
        meta_path = Path(args.metadata)
        if not meta_path.exists():
            print(f"\nError: metadata file not found: {args.metadata}")
            return 1
        with open(meta_path) as fh:
            file_meta = yaml.safe_load(fh) or {}
        metadata.update(file_meta)

    sample_id = metadata.get("sample_id", "UNKNOWN")
    out_dir   = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Run log ---
    from reporting.run_log import RunLog
    log_path = out_dir / f"{sample_id}.run.log.json"
    run_log  = RunLog(str(log_path), metadata, paths_config)

    log.info("Analysis mode — input: %s  sample: %s", args.input, sample_id)

    try:
        from analysis.input_router import detect_input_type, InputType
        from analysis.parser import parse_vep_output

        input_type = detect_input_type(args.input)
        log.info("Input type detected: %s", input_type.value)
        run_log.add_step("detect_input", "completed", detail=f"Input type: {input_type.value}")

        # ----------------------------------------------------------------
        # Stage tracker — emits ##STAGE## markers consumed by the GUI.
        # _TAIL = 5 shared stages after vcf_to_parse is determined:
        #   parse · filter · classify · rank · reports
        # ----------------------------------------------------------------
        _TAIL = 5

        class _Prog:
            def __init__(self, total: int) -> None:
                self.n = 0
                self.total = total
            def next(self, label: str) -> None:
                self.n += 1
                _stage(self.n, self.total, label)
            def expand(self) -> None:
                self.total += 1

        # ----------------------------------------------------------------
        # VCF / VCF.GZ path — validate, auto-annotate if CSQ absent
        # ----------------------------------------------------------------
        if input_type in (InputType.VCF, InputType.VCF_GZ):
            from analysis.vcf_validator import (
                validate_vcf, has_csq_annotation, has_multiallelic_records,
            )
            from analysis.vcf_normalizer import normalize_vcf
            from analysis.annotator import annotate

            validate_vcf(args.input, paths_config)

            if has_csq_annotation(args.input):
                log.info("Detected annotated VCF, skipping VEP")
                run_log.add_step("vcf_detection", "completed",
                                 detail="Pre-annotated VCF detected, VEP skipped")
                prog = _Prog(_TAIL)
                vcf_to_parse = args.input
            else:
                log.info("Detected unannotated VCF, running VEP annotation first")
                run_log.add_step("vcf_detection", "completed",
                                 detail="Unannotated VCF detected, running VEP")

                # Peek for multi-allelic records now (fast) so total stage count is known
                # before emitting the first stage marker.
                _needs_norm = has_multiallelic_records(args.input)
                prog = _Prog((2 if _needs_norm else 1) + _TAIL)

                vcf_for_annotation = args.input
                if _needs_norm:
                    log.info("Multi-allelic records detected — normalizing with bcftools norm")
                    prog.next("Normalising VCF")
                    norm_vcf = str(out_dir / f"{sample_id}.norm.vcf")
                    vcf_for_annotation = normalize_vcf(args.input, norm_vcf, paths_config, threads=int(settings.get("vep", {}).get("fork", 4)))
                    run_log.add_output_file("normalized_vcf", vcf_for_annotation)
                    run_log.add_step("normalize_vcf", "completed",
                                     detail=f"Normalized VCF: {vcf_for_annotation}")
                else:
                    log.info("No multi-allelic records detected — skipping normalization")

                # Use a distinct output name when input was normalized so that a
                # stale pre-normalization vep.vcf is never reused by skip_if_annotated_exists.
                vep_suffix = "norm.vep.vcf" if _needs_norm else "vep.vcf"
                annotated_vcf = out_dir / f"{sample_id}.{vep_suffix}"
                prog.next("Running VEP annotation")
                vcf_to_parse = annotate(
                    vcf_for_annotation, str(annotated_vcf), paths_config, settings
                )
                run_log.add_output_file("annotated_vcf", vcf_to_parse)
                run_log.add_step("vep_annotation", "completed",
                                 detail=f"Annotated VCF: {vcf_to_parse}")

        # ----------------------------------------------------------------
        # FASTQ path — align → call variants → normalise? → annotate → parse
        # ----------------------------------------------------------------
        elif input_type == InputType.FASTQ:
            from analysis.vcf_validator import has_multiallelic_records
            from analysis.vcf_normalizer import normalize_vcf
            from analysis.annotator import annotate
            from analysis.upstream.fastq_handler import run_fastq_to_vcf

            log.info("FASTQ input — running alignment and variant calling")
            run_log.add_step("fastq_detected", "completed",
                             detail=f"R1={args.input}  R2={args.input2 or 'none (single-end)'}")

            # 4 FASTQ steps (sync → align → markdup → DV) + VEP + _TAIL = 10 stages.
            # Will expand to 11 if normalisation is needed after DeepVariant.
            prog = _Prog(4 + _TAIL)

            prog.next("Syncing reads  (seqkit)")
            raw_vcf = run_fastq_to_vcf(
                fastq_r1=args.input,
                fastq_r2=args.input2,
                sample_id=sample_id,
                output_dir=str(out_dir),
                paths_config=paths_config,
                settings=settings,
                stage_callbacks=[
                    lambda: prog.next("Aligning reads  (BWA)"),
                    lambda: prog.next("Marking duplicates  (samtools)"),
                    lambda: prog.next("Calling variants  (DeepVariant)"),
                ],
            )
            run_log.add_output_file("raw_vcf", raw_vcf)
            run_log.add_step("fastq_to_vcf", "completed", detail=f"Raw VCF: {raw_vcf}")

            vcf_for_annotation = raw_vcf
            if has_multiallelic_records(raw_vcf):
                log.info(
                    "Multi-allelic records detected in DeepVariant output — "
                    "normalizing with bcftools norm"
                )
                prog.expand()   # total becomes 9
                prog.next("Normalising VCF")
                norm_vcf = str(out_dir / f"{sample_id}.norm.vcf")
                vcf_for_annotation = normalize_vcf(raw_vcf, norm_vcf, paths_config, threads=int(settings.get("vep", {}).get("fork", 4)))
                run_log.add_output_file("normalized_vcf", vcf_for_annotation)
                run_log.add_step("normalize_vcf", "completed",
                                 detail=f"Normalized VCF: {vcf_for_annotation}")
                vep_suffix = "norm.vep.vcf"
            else:
                log.info("No multi-allelic records in DeepVariant output — skipping normalization")
                vep_suffix = "vep.vcf"

            log.info("Annotating VCF produced from FASTQ input")
            annotated_vcf = out_dir / f"{sample_id}.{vep_suffix}"
            prog.next("Running VEP annotation")
            vcf_to_parse = annotate(vcf_for_annotation, str(annotated_vcf), paths_config, settings)
            run_log.add_output_file("annotated_vcf", vcf_to_parse)
            run_log.add_step("vep_annotation", "completed",
                             detail=f"Annotated VCF: {vcf_to_parse}")

        # ----------------------------------------------------------------
        # BAM/CRAM path — index if needed → call variants → normalise? → annotate
        # ----------------------------------------------------------------
        elif input_type == InputType.BAM_CRAM:
            from analysis.vcf_validator import has_multiallelic_records
            from analysis.vcf_normalizer import normalize_vcf
            from analysis.annotator import annotate
            from analysis.upstream.bam_handler import run_bam_to_vcf

            log.info("BAM/CRAM input — running variant calling")
            run_log.add_step("bam_detected", "completed",
                             detail=f"BAM={args.input}")

            # 2 stages minimum (DV + VEP), expands to 3 if normalisation needed.
            prog = _Prog(2 + _TAIL)

            prog.next("Calling variants  (DeepVariant)")
            raw_vcf = run_bam_to_vcf(
                bam_path=args.input,
                sample_id=sample_id,
                output_dir=str(out_dir),
                paths_config=paths_config,
                settings=settings,
            )
            run_log.add_output_file("raw_vcf", raw_vcf)
            run_log.add_step("bam_to_vcf", "completed", detail=f"Raw VCF: {raw_vcf}")

            vcf_for_annotation = raw_vcf
            if has_multiallelic_records(raw_vcf):
                log.info(
                    "Multi-allelic records detected in DeepVariant output — "
                    "normalizing with bcftools norm"
                )
                prog.expand()   # total becomes 3 + _TAIL
                prog.next("Normalising VCF")
                norm_vcf = str(out_dir / f"{sample_id}.norm.vcf")
                vcf_for_annotation = normalize_vcf(raw_vcf, norm_vcf, paths_config, threads=int(settings.get("vep", {}).get("fork", 4)))
                run_log.add_output_file("normalized_vcf", vcf_for_annotation)
                run_log.add_step("normalize_vcf", "completed",
                                 detail=f"Normalized VCF: {vcf_for_annotation}")
                vep_suffix = "norm.vep.vcf"
            else:
                log.info("No multi-allelic records in DeepVariant output — skipping normalization")
                vep_suffix = "vep.vcf"

            log.info("Annotating VCF produced from BAM input")
            annotated_vcf = out_dir / f"{sample_id}.{vep_suffix}"
            prog.next("Running VEP annotation")
            vcf_to_parse = annotate(vcf_for_annotation, str(annotated_vcf), paths_config, settings)
            run_log.add_output_file("annotated_vcf", vcf_to_parse)
            run_log.add_step("vep_annotation", "completed",
                             detail=f"Annotated VCF: {vcf_to_parse}")

        # ----------------------------------------------------------------
        # Unknown extension
        # ----------------------------------------------------------------
        else:
            raise ValueError(
                f"Unrecognised input file type: {args.input}\n"
                f"\n"
                f"  Supported formats:\n"
                f"    .vcf / .vcf.gz               — annotated or unannotated VCF\n"
                f"    .fastq / .fastq.gz / .fq / .fq.gz  — paired or single-end reads\n"
                f"    .bam / .cram                 — sorted BAM or CRAM"
            )

        # ----------------------------------------------------------------
        # Common tail — parse → filter → classify → rank → reports
        # ----------------------------------------------------------------

        # 1. Parse
        prog.next("Parsing variants")
        log.info("Parsing VCF: %s", vcf_to_parse)
        variants = parse_vep_output(vcf_to_parse)
        run_log.add_step("parse_vcf", "completed", detail=f"{len(variants)} variants parsed")
        log.info("Parsed %d variants", len(variants))

        # 1b. ClinGen gene validity lookup (fast — merged into parse stage)
        clingen_entry = paths_config.get("resources", {}).get("clingen_gene_validity_csv", {})
        clingen_csv = clingen_entry.get("path") if clingen_entry.get("status") == "found" else None
        if clingen_csv:
            from analysis.clingen import load_clingen, annotate_clingen
            gene_dict = load_clingen(clingen_csv)
            n_clingen = annotate_clingen(variants, gene_dict)
            run_log.add_step("clingen_lookup", "completed",
                             detail=f"{n_clingen} variants matched a curated gene")
        else:
            log.info("ClinGen gene validity not configured — skipping")
            run_log.add_step("clingen_lookup", "skipped", detail="not configured")

        # 2. Load gene list
        gene_list = []
        gene_list_path = Path(args.gene_list)
        if gene_list_path.exists():
            with open(gene_list_path) as fh:
                gene_list = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
            run_log.add_step("load_gene_list", "completed", detail=f"{len(gene_list)} genes")
        else:
            run_log.add_warning(f"Gene list not found: {args.gene_list}")
            log.warning("Gene list not found: %s", args.gene_list)

        # 3. Filter
        prog.next("Filtering variants")
        from analysis.filter_engine import apply_filters
        variants = apply_filters(variants, settings, gene_list=gene_list)
        n_retained = sum(1 for v in variants if not v.excluded)
        n_excluded = sum(1 for v in variants if v.excluded)
        run_log.add_step("filter", "completed",
                         detail=f"{n_retained} retained, {n_excluded} excluded")
        log.info("Filter: %d retained, %d excluded", n_retained, n_excluded)

        # 4. ACMG classification + reviewer flags
        prog.next("Classifying variants  (ACMG)")
        from analysis.acmg.classifier import classify_all
        classify_all(variants, settings)
        run_log.add_step("acmg_classify", "completed")

        from analysis.reviewer_flags import annotate_reviewer_flags
        annotate_reviewer_flags(variants, settings, gene_list=gene_list)
        run_log.add_step("reviewer_flags", "completed")

        # 5. Rank
        prog.next("Ranking candidates")
        from analysis.ranker import rank_variants
        variants = rank_variants(variants, settings)
        run_log.add_step("rank", "completed")

        # 6. Reports — each format is isolated so one failure never blocks the rest.
        prog.next("Writing reports")

        json_path = out_dir / f"{sample_id}.json"
        tsv_path  = out_dir / f"{sample_id}.variants.tsv"
        xlsx_path = out_dir / f"{sample_id}.xlsx"
        pdf_path  = out_dir / f"{sample_id}.pdf"

        written_json = None
        try:
            from reporting.json_export import write_json
            written_json = write_json(variants, str(json_path), metadata, settings)
            run_log.add_output_file("json", written_json)
            run_log.add_step("report_json", "completed")
        except Exception as exc:
            log.warning("JSON export failed (other reports unaffected): %s", exc)
            run_log.add_step("report_json", "failed", detail=str(exc))

        written_tsv = None
        try:
            from reporting.tsv_export import write_tsv
            written_tsv = write_tsv(variants, str(tsv_path), metadata, settings)
            run_log.add_output_file("tsv", written_tsv)
            run_log.add_step("report_tsv", "completed")
        except Exception as exc:
            log.warning("TSV export failed (other reports unaffected): %s", exc)
            run_log.add_step("report_tsv", "failed", detail=str(exc))

        written_xlsx = None
        try:
            from reporting.spreadsheet import write_spreadsheet
            written_xlsx = write_spreadsheet(variants, str(xlsx_path), metadata, settings)
            run_log.add_output_file("xlsx", written_xlsx)
            run_log.add_step("report_xlsx", "completed")
        except Exception as exc:
            log.warning("XLSX export failed (other reports unaffected): %s", exc)
            run_log.add_step("report_xlsx", "failed", detail=str(exc))

        written_pdf = None
        try:
            from reporting.pdf_report import write_pdf
            written_pdf = write_pdf(variants, str(pdf_path), metadata, settings)
            run_log.add_output_file("pdf", written_pdf)
            run_log.add_step("report_pdf", "completed")
        except Exception as exc:
            log.warning("PDF report failed (analysis data still in JSON/TSV): %s", exc)
            run_log.add_step("report_pdf", "failed", detail=str(exc))

        run_log.set_summary({
            "total_variants": len(variants),
            "retained":       n_retained,
            "excluded":       n_excluded,
        })

        # Register run log path before finalize() so it appears in the JSON.
        run_log.add_output_file("run_log", str(log_path))
        log_written = run_log.finalize()

        print(f"\nVI_Pipeline complete.")
        print(f"  TSV:     {written_tsv  or '(failed — see log)'}")
        print(f"  JSON:    {written_json or '(failed — see log)'}")
        print(f"  XLSX:    {written_xlsx or '(failed — see log)'}")
        print(f"  PDF:     {written_pdf  or '(failed — see log)'}")
        print(f"  Run log: {log_written}")
        return 0

    except Exception as exc:
        run_log.add_error(f"{type(exc).__name__}: {exc}")
        run_log.finalize()
        log.error("Pipeline failed: %s", exc, exc_info=args.verbose)
        print(f"\nError: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
