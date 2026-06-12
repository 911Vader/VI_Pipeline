"""
vcf_validator.py — Validate an input VCF or VCF.GZ before annotation.

Version 1 checks (in order):
  1. File exists
  2. Path is a regular file, not a directory
  3. File is readable
  4. File is not empty
  5. If VCF.GZ: tabix index (.tbi) is present alongside the file
  6. First line of the VCF starts with ##fileformat=VCF

What version 1 does NOT do:
  - Left-normalisation (requires bcftools — deferred to a later phase)
  - Contig name validation (chr-prefix check — deferred)
  - Sample column presence check — deferred

paths_config is accepted now so the signature is stable for future phases
when bcftools normalisation will be added. It is not used in version 1.
"""

import gzip
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def validate_vcf(vcf_path: str, paths_config: Optional[dict] = None) -> str:
    """
    Validate an input VCF or VCF.GZ file.

    Parameters
    ----------
    vcf_path:    Absolute or relative path to the input VCF.
    paths_config: Validated tool/resource paths from config/paths.yaml.
                  Accepted for future normalisation steps; unused in version 1.

    Returns
    -------
    The same vcf_path string, unchanged, if all checks pass.

    Raises
    ------
    ValueError: With a beginner-friendly message describing exactly what is wrong
                and how to fix it.
    """
    _check_exists(vcf_path)
    _check_is_file(vcf_path)
    _check_readable(vcf_path)
    _check_not_empty(vcf_path)

    is_gz = _is_compressed(vcf_path)

    if is_gz:
        _check_tabix_index(vcf_path)

    _check_vcf_header(vcf_path, is_gz)

    logger.info("VCF validation passed: %s", vcf_path)
    return vcf_path


# ---------------------------------------------------------------------------
# Internal checks — each raises ValueError with a specific message on failure
# ---------------------------------------------------------------------------

def _is_compressed(vcf_path: str) -> bool:
    """Return True if the path ends with .gz or .bgz."""
    lower = vcf_path.lower()
    return lower.endswith(".vcf.gz") or lower.endswith(".vcf.bgz")


def _check_exists(vcf_path: str) -> None:
    if not Path(vcf_path).exists():
        raise ValueError(
            f"Input file not found: {vcf_path}\n"
            f"  Check that the path is correct and the file exists."
        )


def _check_is_file(vcf_path: str) -> None:
    if not Path(vcf_path).is_file():
        raise ValueError(
            f"Input path is not a file: {vcf_path}\n"
            f"  Expected a VCF (.vcf) or compressed VCF (.vcf.gz) file.\n"
            f"  A directory was found at this path instead."
        )


def _check_readable(vcf_path: str) -> None:
    if not os.access(vcf_path, os.R_OK):
        raise ValueError(
            f"Input file cannot be read: {vcf_path}\n"
            f"  Check file permissions.\n"
            f"  Try: ls -l {vcf_path}"
        )


def _check_not_empty(vcf_path: str) -> None:
    if Path(vcf_path).stat().st_size == 0:
        raise ValueError(
            f"Input file is empty (0 bytes): {vcf_path}\n"
            f"  The file exists but contains no data."
        )


def _check_tabix_index(vcf_gz_path: str) -> None:
    """
    Confirm a .tbi tabix index exists alongside the bgzipped VCF.

    VEP and bcftools both require a tabix index when reading bgzipped input.
    If the index is missing and tabix is on PATH, it is created automatically.
    """
    tbi = Path(vcf_gz_path + ".tbi")
    if tbi.exists():
        logger.debug("Tabix index found: %s", tbi)
        return

    # Attempt to auto-create the index — common case when the user downloads
    # a VCF without its accompanying .tbi file.
    tabix_exe = shutil.which("tabix")
    if tabix_exe:
        logger.info(
            "Tabix index not found — creating automatically: %s.tbi "
            "(this may take a moment for large files)",
            vcf_gz_path,
        )
        try:
            result = subprocess.run(
                [tabix_exe, "-p", "vcf", vcf_gz_path],
                capture_output=True,
                text=True,
                timeout=600,
            )
            if result.returncode == 0 and tbi.exists():
                logger.info("Tabix index created: %s", tbi)
                return
            logger.warning(
                "tabix exited %d while creating index: %s",
                result.returncode, result.stderr.strip(),
            )
        except subprocess.TimeoutExpired:
            logger.warning("tabix timed out while creating index for %s", vcf_gz_path)
        except Exception as exc:
            logger.warning("tabix auto-creation failed: %s", exc)

    raise ValueError(
        f"Tabix index not found for compressed VCF.\n"
        f"  VCF file : {vcf_gz_path}\n"
        f"  Expected : {tbi}\n"
        f"\n"
        f"  To create the index, run:\n"
        f"    tabix -p vcf {vcf_gz_path}\n"
        f"\n"
        f"  If the file was not compressed with bgzip, you will need to recompress it:\n"
        f"    bgzip {vcf_gz_path.rstrip('.gz')}\n"
        f"    tabix -p vcf {vcf_gz_path}"
    )


def _check_vcf_header(vcf_path: str, is_gz: bool) -> None:
    """
    Read the first line and confirm it starts with '##fileformat=VCF'.

    This is the mandatory first line of any valid VCF file (VCF spec §1.1).
    """
    try:
        if is_gz:
            with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
                first_line = fh.readline().rstrip("\n")
        else:
            with open(vcf_path, "r", encoding="utf-8", errors="replace") as fh:
                first_line = fh.readline().rstrip("\n")
    except (OSError, gzip.BadGzipFile) as exc:
        raise ValueError(
            f"Could not read input file: {vcf_path}\n"
            f"  Error: {exc}\n"
            f"\n"
            f"  If this is a bgzipped VCF, confirm it was compressed with bgzip "
            f"(not standard gzip):\n"
            f"    file {vcf_path}"
        ) from exc

    if not first_line.startswith("##fileformat=VCF"):
        # Provide a helpful preview of what was found
        preview = first_line[:80] if first_line else "(empty line)"
        raise ValueError(
            f"Input file does not appear to be a valid VCF: {vcf_path}\n"
            f"  The first line of a valid VCF must start with '##fileformat=VCF'\n"
            f"  Found instead: {preview!r}\n"
            f"\n"
            f"  Common causes:\n"
            f"    - File is not a VCF (wrong extension or wrong file)\n"
            f"    - File is bgzipped but was opened without .gz extension\n"
            f"    - File header was accidentally removed or truncated"
        )

    logger.debug("VCF header OK: %s", first_line.strip())


# ---------------------------------------------------------------------------
# CSQ detection
# ---------------------------------------------------------------------------

def has_csq_annotation(vcf_path: str) -> bool:
    """
    Return True if the VCF header contains a VEP CSQ INFO line.

    Reads only the header lines (stops at the first non-## line) so it is
    fast even on large files.  Safe to call on plain or bgzipped VCF files.

    This is used by the CLI to decide whether VEP annotation is needed before
    parsing.  A VCF is considered annotated when its header contains a line
    matching:  ##INFO=<ID=CSQ,...>
    """
    is_gz = _is_compressed(vcf_path)
    try:
        if is_gz:
            fh = gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace")
        else:
            fh = open(vcf_path, "r", encoding="utf-8", errors="replace")
        with fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                if line.startswith("##INFO=<ID=CSQ"):
                    logger.debug("CSQ annotation header found: %s", vcf_path)
                    return True
    except (OSError, gzip.BadGzipFile) as exc:
        logger.warning("Could not read VCF header for CSQ detection: %s — %s", vcf_path, exc)
    return False


# ---------------------------------------------------------------------------
# Multi-allelic detection
# ---------------------------------------------------------------------------

def has_multiallelic_records(vcf_path: str, max_records: int = 10_000) -> bool:
    """
    Return True if any of the first max_records data lines has a multi-allelic
    ALT field (i.e. the ALT column contains a comma, e.g. 'T,G').

    Reads only up to max_records data lines for performance on large files.
    Safe to call on plain or bgzipped VCF files.
    Returns False on any read error — never raises.

    This is used by the CLI to decide whether bcftools norm is needed before
    VEP annotation.  Only called for raw (unannotated) VCFs.
    """
    is_gz = _is_compressed(vcf_path)
    try:
        open_fn = gzip.open if is_gz else open
        with open_fn(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            checked = 0
            for line in fh:
                if line.startswith("#"):
                    continue
                # VCF columns: CHROM POS ID REF ALT ...
                parts = line.split("\t", 5)
                if len(parts) >= 5 and "," in parts[4]:
                    logger.debug(
                        "Multi-allelic record found in %s at %s:%s ALT=%s",
                        vcf_path, parts[0], parts[1], parts[4].rstrip(),
                    )
                    return True
                checked += 1
                if checked >= max_records:
                    break
    except (OSError, gzip.BadGzipFile) as exc:
        logger.warning(
            "Could not read VCF for multi-allelic detection: %s — %s", vcf_path, exc
        )
    return False
