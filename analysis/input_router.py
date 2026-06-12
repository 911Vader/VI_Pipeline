"""
input_router.py — Detect input file type and route to the appropriate handler.

Supported input types:
  - VCF (.vcf)                              → validated and passed to annotation
  - VCF.GZ / VCF.BGZ (.vcf.gz, .vcf.bgz)  → validated + tabix index checked
  - BCF (.bcf)                              → treated as VCF_GZ
  - FASTQ (.fastq, .fq, .fastq.gz, .fq.gz) → handled by fastq_handler.py in cli.py
  - BAM / CRAM (.bam, .cram)               → handled by bam_handler.py in cli.py

Note on route_input():
  route_input() is a VCF-only routing helper used by legacy tests.
  In v1.2 and later, cli.py calls detect_input_type() directly and
  dispatches to fastq_handler.py for FASTQ input and bam_handler.py for BAM.
  route_input() still raises NotImplementedError for FASTQ/BAM to preserve
  that contract.

All other extensions raise ValueError.
"""

import logging
from enum import Enum
from pathlib import Path
from typing import Optional

from analysis.vcf_validator import validate_vcf

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Input type catalogue
# ---------------------------------------------------------------------------

class InputType(Enum):
    VCF      = "vcf"
    VCF_GZ   = "vcf_gz"
    FASTQ    = "fastq"
    BAM_CRAM = "bam_cram"
    UNKNOWN  = "unknown"


# Extension → InputType mapping.
# Checked longest-first so ".vcf.gz" is matched before ".gz".
_EXTENSION_MAP = [
    (".vcf.gz",    InputType.VCF_GZ),
    (".vcf.bgz",   InputType.VCF_GZ),
    (".vcf",       InputType.VCF),
    (".bcf",       InputType.VCF_GZ),   # binary VCF; treated as compressed VCF
    (".fastq.gz",  InputType.FASTQ),
    (".fq.gz",     InputType.FASTQ),
    (".fastq",     InputType.FASTQ),
    (".fq",        InputType.FASTQ),
    (".bam",       InputType.BAM_CRAM),
    (".cram",      InputType.BAM_CRAM),
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def detect_input_type(input_path: str) -> InputType:
    """
    Determine the input file type from the file extension.

    Does not access the filesystem — extension-only detection.
    Case-insensitive.

    Parameters
    ----------
    input_path: Path string (absolute or relative).

    Returns
    -------
    InputType enum value.
    """
    lower = Path(input_path).name.lower()
    for suffix, input_type in _EXTENSION_MAP:
        if lower.endswith(suffix):
            return input_type
    return InputType.UNKNOWN


def route_input(input_path: str, paths_config: Optional[dict] = None) -> str:
    """
    Detect input type, validate it, and return a VCF path ready for annotation.

    Parameters
    ----------
    input_path:   Path to the input file supplied by the user.
    paths_config: Validated paths from config/paths.yaml. Passed through to
                  vcf_validator for future normalisation steps.

    Returns
    -------
    Absolute path string to a validated VCF (or VCF.GZ) file.

    Raises
    ------
    NotImplementedError: For FASTQ and BAM/CRAM inputs (version 1 limitation).
    ValueError:          For unrecognised file types, or if VCF validation fails.
    """
    input_type = detect_input_type(input_path)
    logger.info("Detected input type: %s for %s", input_type.value, input_path)

    if input_type in (InputType.VCF, InputType.VCF_GZ):
        validated_path = validate_vcf(input_path, paths_config)
        logger.info("Input routed to VCF analysis path: %s", validated_path)
        return validated_path

    if input_type == InputType.FASTQ:
        raise NotImplementedError(
            f"FASTQ input is not supported in version 1 of this pipeline.\n"
            f"  Provided: {input_path}\n"
            f"\n"
            f"  Version 1 accepts VCF or VCF.GZ files only.\n"
            f"  If you have a FASTQ file, you will need to align and call variants\n"
            f"  first to produce a VCF, then re-run the pipeline with that VCF.\n"
            f"\n"
            f"  FASTQ support is planned for a future version."
        )

    if input_type == InputType.BAM_CRAM:
        raise NotImplementedError(
            f"BAM/CRAM input is not supported in version 1 of this pipeline.\n"
            f"  Provided: {input_path}\n"
            f"\n"
            f"  Version 1 accepts VCF or VCF.GZ files only.\n"
            f"  If you have a BAM or CRAM file, you will need to call variants\n"
            f"  first to produce a VCF, then re-run the pipeline with that VCF.\n"
            f"\n"
            f"  BAM/CRAM support is planned for a future version."
        )

    # InputType.UNKNOWN
    suffix = Path(input_path).suffix or "(no extension)"
    raise ValueError(
        f"Unrecognised input file type: {input_path}\n"
        f"  Extension detected: {suffix}\n"
        f"\n"
        f"  Supported in version 1:\n"
        f"    .vcf        — plain VCF\n"
        f"    .vcf.gz     — bgzipped VCF (requires .tbi tabix index)\n"
        f"\n"
        f"  Not yet supported (planned for future versions):\n"
        f"    .fastq / .fq / .fastq.gz / .fq.gz\n"
        f"    .bam / .cram"
    )
