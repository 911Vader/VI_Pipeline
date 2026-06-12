"""
vcf_normalizer.py — Split multi-allelic VCF records and left-align indels.

Uses bcftools norm.  Only invoked for raw (unannotated) VCFs where
has_multiallelic_records() returns True.

Design decisions:
  - bcftools is used because it is already a required tool for variant calling
    and is validated during setup.  No extra dependency is introduced.
  - The GRCh38 reference FASTA is taken from paths_config so we use exactly
    the same reference that VEP will use for annotation.
  - '-m -any' splits ALL multi-allelic records into separate biallelic lines,
    which is what VEP and the downstream parser expect.
  - Left-alignment of indels ('-f {fasta}') is included so that the
    normalised VCF is coordinate-consistent with the reference cache.
  - Output is always a plain VCF (not bgzipped); VEP accepts plain VCF.

Output consumed by:
  analysis/annotator.py  →  VEP annotation
"""

import logging
import shlex
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

_NORM_TIMEOUT = 7200   # 2 hours; bcftools norm on a WGS VCF is fast but allow margin


# ---------------------------------------------------------------------------
# Tool / resource resolution
# ---------------------------------------------------------------------------

def _get_bcftools(paths_config: dict) -> str:
    """
    Return the bcftools executable path from paths_config, or locate on PATH.

    Raises RuntimeError with install guidance if not found anywhere.
    """
    entry = paths_config.get("tools", {}).get("bcftools", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        logger.debug("bcftools from paths_config: %s", entry["path"])
        return entry["path"]
    path = shutil.which("bcftools")
    if path:
        logger.debug("bcftools located on PATH: %s", path)
        return path
    raise RuntimeError(
        "bcftools is required for VCF normalization but was not found.\n"
        "\n"
        "  Install it and ensure it is on your PATH:\n"
        "    conda install -c bioconda bcftools\n"
        "\n"
        "  After installation, re-run setup:\n"
        "    python cli.py --setup"
    )


def _get_reference_fasta(paths_config: dict) -> str:
    """
    Return the GRCh38 reference FASTA path from paths_config.

    Raises RuntimeError if not configured or not found.
    """
    entry = paths_config.get("resources", {}).get("reference_fasta", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        return entry["path"]
    raise RuntimeError(
        "GRCh38 reference FASTA is required for VCF normalization but is not configured.\n"
        "\n"
        "  Add the path to config/paths.yaml under resources.reference_fasta.path\n"
        "  and re-run setup:\n"
        "    python cli.py --setup"
    )


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

def normalize_vcf(input_vcf: str, output_vcf: str, paths_config: dict, threads: int = 4) -> str:
    """
    Normalize a raw VCF using bcftools norm.

    Splits multi-allelic records into separate biallelic records and
    left-aligns indels against the GRCh38 reference FASTA.

    Parameters
    ----------
    input_vcf:    Path to the raw (unannotated) input VCF.
    output_vcf:   Path to write the normalized output VCF.  The parent
                  directory is created if it does not exist.
    paths_config: Parsed config/paths.yaml.

    Returns
    -------
    Path to the normalized output VCF (same as output_vcf).

    Raises
    ------
    RuntimeError: If bcftools or the reference FASTA are missing, if
                  normalization times out, or if bcftools exits non-zero.
    """
    bcftools = _get_bcftools(paths_config)
    fasta    = _get_reference_fasta(paths_config)

    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        bcftools, "norm",
        "--threads", str(threads),
        "-m", "-any",       # split all multi-allelic into biallelic records
        "-f", fasta,        # left-align and normalise indels against reference
        "-o", output_vcf,   # write output here (plain VCF)
        input_vcf,
    ]

    logger.info("Normalizing VCF with bcftools norm:\n  %s", " ".join(shlex.quote(t) for t in cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=_NORM_TIMEOUT,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"bcftools not found at runtime: {bcftools}\n"
            f"  Re-run setup to re-validate:  python cli.py --setup"
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"bcftools norm timed out after {_NORM_TIMEOUT // 3600} hour(s).\n"
            f"\n"
            f"  Input: {input_vcf}\n"
            f"\n"
            f"  If this is a very large VCF, consider splitting by chromosome first:\n"
            f"    bcftools view -r chr1 {input_vcf} -Oz -o chr1.vcf.gz\n"
            f"    tabix -p vcf chr1.vcf.gz"
        )

    if result.stderr:
        logger.debug("bcftools norm stderr:\n%s", result.stderr.strip())

    if result.returncode != 0:
        stderr_tail = "\n    ".join(
            result.stderr.strip().splitlines()[-10:]
        ) or "(no output)"
        raise RuntimeError(
            f"bcftools norm failed (exit code {result.returncode}).\n"
            f"\n"
            f"  Output:\n"
            f"    {stderr_tail}\n"
            f"\n"
            f"  Common causes:\n"
            f"    - Reference FASTA chromosome names do not match the VCF\n"
            f"      (e.g. 'chr1' vs '1') — check with: bcftools view -h {input_vcf}\n"
            f"    - Input VCF is corrupt or truncated\n"
            f"\n"
            f"  Re-run setup to validate all paths:\n"
            f"    python cli.py --setup"
        )

    if not Path(output_vcf).exists():
        raise RuntimeError(
            f"bcftools norm exited successfully but output file was not produced.\n"
            f"  Expected: {output_vcf}\n"
            f"  Check that the output directory is writable and has sufficient disk space."
        )

    logger.info("Normalization complete. Output: %s", output_vcf)
    return output_vcf
