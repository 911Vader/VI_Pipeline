"""
validation.py — Validate GUI form inputs before launching the pipeline.

Pure functions (no tkinter state) so they can be unit-tested without a display.
"""

from pathlib import Path

_FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
_VCF_EXTS   = (".vcf", ".vcf.gz", ".bcf")
_BAM_EXTS   = (".bam", ".cram")


def _is_fastq(path: str) -> bool:
    p = path.lower()
    return any(p.endswith(e) for e in _FASTQ_EXTS)


def _is_vcf(path: str) -> bool:
    p = path.lower()
    return any(p.endswith(e) for e in _VCF_EXTS)


def _is_bam(path: str) -> bool:
    p = path.lower()
    return any(p.endswith(e) for e in _BAM_EXTS)


def validate_inputs(inputs: dict) -> tuple:
    """
    Check that all required fields are present and that file paths exist.

    Parameters
    ----------
    inputs:
        Dict with string values for keys: input, input2, gene_list,
        sample_id, case_id.  Missing keys default to "".
        output_dir is derived automatically as runs/{sample_id}/.

    Returns
    -------
    (True, "")           — all required fields are valid.
    (False, message)     — a required field is missing or invalid.
    """
    input_file = inputs.get("input", "").strip()
    gene_list  = inputs.get("gene_list", "").strip()
    sample_id  = inputs.get("sample_id", "").strip()
    input2     = inputs.get("input2", "").strip()

    if not input_file:
        return False, "Input file is required."
    if not Path(input_file).exists():
        return False, f"Input file not found:\n{input_file}"
    if not _is_fastq(input_file) and not _is_vcf(input_file) and not _is_bam(input_file):
        return False, (
            "Unrecognised input file type.\n"
            "Expected a VCF (.vcf, .vcf.gz, .bcf), "
            "FASTQ (.fastq, .fastq.gz, .fq, .fq.gz), "
            "or BAM/CRAM (.bam, .cram) file."
        )

    if input2:
        if not _is_fastq(input_file):
            return False, (
                "R2 FASTQ is only valid for FASTQ input.\n"
                "Remove the R2 file or change the input to a FASTQ file."
            )
        if not Path(input2).exists():
            return False, f"R2 FASTQ file not found:\n{input2}"

    if not gene_list:
        return False, "Gene list is required."
    if not Path(gene_list).exists():
        return False, f"Gene list file not found:\n{gene_list}"

    if not sample_id:
        return False, "Sample ID is required."

    return True, ""
