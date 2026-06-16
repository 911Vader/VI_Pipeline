"""
bam_handler.py — Call variants from a pre-aligned BAM/CRAM using DeepVariant.

The BAM must be coordinate-sorted against GRCh38. A .bai index is created
automatically if one is not already present alongside the BAM.

Usage in cli.py:
    from analysis.upstream.bam_handler import run_bam_to_vcf
    raw_vcf = run_bam_to_vcf(bam_path, sample_id, output_dir, paths_config, settings)
"""

import logging
import subprocess
from pathlib import Path

from analysis.upstream.fastq_handler import (
    _get_docker,
    _get_reference_fasta,
    _require_tool,
    _run_deepvariant,
    _run_flagstat,
    _run_samtools_index,
)

logger = logging.getLogger(__name__)

_DV_SHARDS_DEFAULT  = 4
_DV_VERSION_DEFAULT = "1.10.0"
_DV_MODEL_DEFAULT   = "WGS"
_DV_TIMEOUT_DEFAULT = 43200   # 12 h


def _check_bam_sort_order(samtools_exe: str, bam_path: str) -> None:
    """
    Raise RuntimeError if the BAM is not coordinate-sorted.

    Reads the @HD header line and inspects the SO (sort order) tag.
    If the tag is absent or the header is missing, logs a warning rather
    than failing — some valid BAMs omit @HD entirely.
    """
    result = subprocess.run(
        [samtools_exe, "view", "-H", bam_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Could not read BAM header: {bam_path}\n"
            f"  {result.stderr.strip()}"
        )
    for line in result.stdout.splitlines():
        if not line.startswith("@HD"):
            continue
        fields = {
            kv.split(":", 1)[0]: kv.split(":", 1)[1]
            for kv in line.split("\t")[1:]
            if ":" in kv
        }
        so = fields.get("SO", "").lower()
        if so == "coordinate":
            return
        if so:
            raise RuntimeError(
                f"BAM is not coordinate-sorted (SO:{so}): {bam_path}\n"
                f"  DeepVariant requires a coordinate-sorted BAM.\n"
                f"  Sort with:  samtools sort -@ 8 -o sorted.bam {bam_path}"
            )
        logger.warning(
            "BAM @HD header has no SO tag — cannot verify sort order: %s", bam_path,
        )
        return
    logger.warning(
        "BAM has no @HD header line — cannot verify sort order: %s", bam_path,
    )


def run_bam_to_vcf(
    bam_path: str,
    sample_id: str,
    output_dir: str,
    paths_config: dict,
    settings: dict,
) -> str:
    """
    Call variants from a sorted BAM/CRAM with DeepVariant.

    Parameters
    ----------
    bam_path:     Path to the coordinate-sorted BAM (or CRAM) file.
    sample_id:    Used to name output files.
    output_dir:   Directory where outputs are written.
    paths_config: Parsed config/paths.yaml.
    settings:     Parsed config/settings.yaml.

    Returns
    -------
    Host path to the DeepVariant output VCF (.dv.vcf.gz).
    """
    samtools_exe = _require_tool("samtools")
    docker_exe   = _get_docker(paths_config)
    fasta_path   = _get_reference_fasta(paths_config)

    dv_cfg        = settings.get("deepvariant", {})
    dv_threads    = int(dv_cfg.get("num_shards",    _DV_SHARDS_DEFAULT))
    version       = str(dv_cfg.get("version",        _DV_VERSION_DEFAULT))
    model_type    = str(dv_cfg.get("model_type",     _DV_MODEL_DEFAULT))
    timeout_h     = float(dv_cfg.get("timeout_hours", _DV_TIMEOUT_DEFAULT / 3600))
    timeout       = int(timeout_h * 3600) if timeout_h > 0 else None
    visual_report = bool(dv_cfg.get("visual_report", False))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    logger.info(
        "BAM -> VCF  sample=%s  caller=DeepVariant-%s  shards=%d",
        sample_id, version, dv_threads,
    )

    # Validate sort order before DeepVariant — a name-sorted BAM produces a
    # cryptic DV error with no hint about the real cause.
    _check_bam_sort_order(samtools_exe, bam_path)

    # Create BAM index if neither <name>.bam.bai nor <name>.bai exists.
    bam = Path(bam_path)
    has_index = (bam.parent / (bam.name + ".bai")).exists() or \
                bam.with_suffix(".bai").exists()
    if not has_index:
        logger.info("BAM index not found — creating with samtools index: %s", bam_path)
        _run_samtools_index(samtools_exe, bam_path)

    # BAM QC before the expensive DeepVariant step.
    _run_flagstat(samtools_exe, bam_path, str(out), sample_id)

    raw_vcf = _run_deepvariant(
        docker_exe=docker_exe,
        fasta_path=fasta_path,
        bam_path=bam_path,
        out_dir=str(out),
        sample_id=sample_id,
        num_shards=dv_threads,
        model_type=model_type,
        version=version,
        timeout=timeout,
        visual_report=visual_report,
    )

    logger.info("BAM variant calling complete: %s", raw_vcf)
    return raw_vcf
