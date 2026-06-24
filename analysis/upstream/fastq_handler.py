"""
fastq_handler.py -- FASTQ to sorted BAM to VCF for VI_Pipeline.

Pipeline steps (all via subprocess):
  1.  Validate required tools (bwa, samtools) and paths_config resources.
  1b. [paired-end only] fastp pair synchronisation — discards orphaned reads so
      BWA never encounters mismatched read names across R1 / R2.
  2.  BWA MEM alignment  ->  samtools sort  ->  sorted BAM.
  3.  samtools index     ->  BAM index (.bai).
  3b. Duplicate marking/removal (position-based, coordinate-sorted BAM):
        Tool preference: sambamba  (fastest, fully multithreaded)
                      → picard     (standard, detailed metrics)
                      → samtools   (always available, no extra install needed)
      Controlled via settings.yaml deepvariant.markdup_tool and
      deepvariant.remove_duplicates.  Output: {sample_id}.dedup.bam.
  4.  DeepVariant (Docker) -> VCF.GZ with high-accuracy variant calls.

Design decisions:
  - bwa and samtools are located via shutil.which() at call time.
  - Docker is read from paths_config (validated at setup time).
  - DeepVariant is the variant caller.  It runs via the official Docker
    image (google/deepvariant:<version>).  The image is pulled automatically
    on first use if not already cached (~1 GB download).
  - Volume mounts: the output directory is mounted as /workdir (rw); the
    reference FASTA directory is mounted as /ref (ro) if it differs from
    the output directory.  No other host paths enter the container.
  - The reference FASTA and BWA index are read from paths_config.
  - All intermediate files are written under output_dir / sample_id to
    avoid collisions between parallel runs.

Output consumed by:
  analysis/annotator.py  ->  VEP annotation
  analysis/parser.py     ->  variant parsing
"""

import gzip
import logging
import os
import re
import threading
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

# Timeouts per step (seconds). 0 = no limit.
_FASTP_TIMEOUT      = 0
_ALIGN_TIMEOUT      = 0
_SORT_TIMEOUT       = 0
_MARKDUP_TIMEOUT    = 0
_DV_TIMEOUT_DEFAULT = 0

_DV_IMAGE = "google/deepvariant"
_DV_VERSION_DEFAULT = "1.10.0"
_DV_MODEL_DEFAULT = "WGS"
_DV_SHARDS_DEFAULT = 4
_DV_USE_GPU_DEFAULT = "auto"


# ---------------------------------------------------------------------------
# Tool / resource validation
# ---------------------------------------------------------------------------

def _require_tool(name: str) -> str:
    """
    Locate a required tool on PATH via shutil.which().

    Returns the absolute path.
    Raises RuntimeError with install guidance if absent.
    """
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"'{name}' is required for FASTQ processing but was not found on PATH.\n"
            f"\n"
            f"  Install it and ensure it is on your PATH, for example:\n"
            f"    conda install -c bioconda {name}\n"
            f"  or follow the tool's official installation instructions.\n"
            f"\n"
            f"  After installation, re-run setup:\n"
            f"    python cli.py --setup"
        )
    logger.debug("Tool located: %s -> %s", name, path)
    return path


def _get_docker(paths_config: dict) -> str:
    """Return the Docker executable path from paths_config or PATH."""
    entry = paths_config.get("tools", {}).get("docker", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        return entry["path"]
    return _require_tool("docker")


def _get_bcftools(paths_config: dict) -> str:
    """Return the bcftools executable path from paths_config or PATH."""
    entry = paths_config.get("tools", {}).get("bcftools", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        return entry["path"]
    return _require_tool("bcftools")


def _get_reference_fasta(paths_config: dict) -> str:
    """
    Return the GRCh38 reference FASTA path from paths_config.

    Raises RuntimeError if not configured or not found.
    """
    entry = paths_config.get("resources", {}).get("reference_fasta", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        return entry["path"]
    raise RuntimeError(
        "GRCh38 reference FASTA is required for FASTQ alignment but is not configured.\n"
        "\n"
        "  Add the path to config/paths.yaml under resources.reference_fasta.path\n"
        "  and re-run setup:\n"
        "    python cli.py --setup"
    )


def _get_bwa_index_prefix(paths_config: dict, fasta_path: str) -> str:
    """
    Return the BWA index prefix.

    Tries paths_config.resources.bwa_index_prefix first; falls back to
    checking whether the standard BWA index files (.amb, .ann, .bwt, .pac,
    .sa) exist alongside the reference FASTA.

    Raises RuntimeError if no index is found.
    """
    entry = paths_config.get("resources", {}).get("bwa_index_prefix", {})
    if isinstance(entry, dict) and entry.get("status") == "found" and entry.get("path"):
        prefix = entry["path"]
        logger.debug("Using bwa_index_prefix from paths_config: %s", prefix)
        return prefix

    required_suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    missing = [s for s in required_suffixes if not Path(fasta_path + s).exists()]
    if not missing:
        logger.debug("BWA index found alongside reference FASTA: %s", fasta_path)
        return fasta_path

    raise RuntimeError(
        f"BWA index not found for reference FASTA: {fasta_path}\n"
        f"\n"
        f"  Missing index files: {', '.join(missing)}\n"
        f"\n"
        f"  To build the index (this takes ~1 hour for GRCh38):\n"
        f"    bwa index {fasta_path}\n"
        f"\n"
        f"  Alternatively, add the path to a pre-built BWA index prefix to\n"
        f"  config/paths.yaml under resources.bwa_index_prefix.path, then\n"
        f"  re-run setup:\n"
        f"    python cli.py --setup"
    )


# ---------------------------------------------------------------------------
# Sub-steps
# ---------------------------------------------------------------------------

def _run_markdup(
    samtools_exe: str,
    input_bam: str,
    output_bam: str,
    metrics_file: str,
    remove_duplicates: bool,
    preferred_tool: str,
    threads: int,
    tmp_prefix: str,
    timeout: int = _MARKDUP_TIMEOUT,
) -> str:
    """
    Mark (and optionally remove) PCR/optical duplicates from a coordinate-sorted BAM.

    Tool selection order (first found wins):
      1. sambamba  — fastest; fully multithreaded duplicate detection
      2. picard    — standard; detailed metrics file
      3. samtools  — always available; position-based (-m t, template start/end)

    All three use coordinate-sorted input and produce a coordinate-sorted output.
    Returns the path to the output dedup BAM.
    """
    preferred = preferred_tool.strip().lower()

    # --- Try sambamba ---
    if preferred in ("sambamba", "auto"):
        sambamba_exe = shutil.which("sambamba")
        if sambamba_exe:
            return _markdup_sambamba(
                sambamba_exe, input_bam, output_bam,
                remove_duplicates=remove_duplicates,
                threads=threads,
                tmp_dir=str(Path(tmp_prefix).parent),
                timeout=timeout,
            )
        if preferred == "sambamba":
            raise RuntimeError(
                "sambamba is configured as the markdup tool but was not found on PATH.\n"
                "  Install:  conda install -c bioconda sambamba\n"
                "  Or change deepvariant.markdup_tool in config/settings.yaml to 'picard' or 'samtools'."
            )

    # --- Try picard ---
    if preferred in ("picard", "auto"):
        picard_exe = shutil.which("picard")
        if picard_exe:
            return _markdup_picard(
                picard_exe, input_bam, output_bam,
                metrics_file=metrics_file,
                remove_duplicates=remove_duplicates,
                timeout=timeout,
            )
        if preferred == "picard":
            raise RuntimeError(
                "picard is configured as the markdup tool but was not found on PATH.\n"
                "  Install:  conda install -c bioconda picard\n"
                "  Or change deepvariant.markdup_tool in config/settings.yaml to 'sambamba' or 'samtools'."
            )

    # --- samtools fallback (always available) ---
    return _markdup_samtools(
        samtools_exe, input_bam, output_bam,
        metrics_file=metrics_file,
        remove_duplicates=remove_duplicates,
        threads=threads,
        tmp_prefix=tmp_prefix,
        timeout=timeout,
    )


def _markdup_sambamba(
    sambamba_exe: str,
    input_bam: str,
    output_bam: str,
    remove_duplicates: bool,
    threads: int,
    tmp_dir: str,
    timeout: int,
) -> str:
    """Duplicate marking with sambamba markdup (fastest — fully multithreaded)."""
    cmd = [
        sambamba_exe, "markdup",
        "--nthreads", str(threads),
        "--tmpdir", tmp_dir,
    ]
    if remove_duplicates:
        cmd.append("--remove-duplicates")
    cmd += [input_bam, output_bam]

    logger.info(
        "sambamba markdup (threads=%d, remove_duplicates=%s): %s → %s",
        threads, remove_duplicates, input_bam, output_bam,
    )
    _run_markdup_cmd(cmd, "sambamba markdup", timeout)
    if not Path(output_bam).exists():
        raise RuntimeError(f"sambamba markdup finished but output BAM not found: {output_bam}")
    logger.info("sambamba markdup complete: %s", output_bam)
    return output_bam


def _markdup_picard(
    picard_exe: str,
    input_bam: str,
    output_bam: str,
    metrics_file: str,
    remove_duplicates: bool,
    timeout: int,
) -> str:
    """Duplicate marking with Picard MarkDuplicates (standard — detailed metrics)."""
    remove_str = "true" if remove_duplicates else "false"
    cmd = [
        picard_exe, "MarkDuplicates",
        f"INPUT={input_bam}",
        f"OUTPUT={output_bam}",
        f"METRICS_FILE={metrics_file}",
        f"REMOVE_DUPLICATES={remove_str}",
        "VALIDATION_STRINGENCY=LENIENT",
        "TMP_DIR=/tmp",
    ]
    logger.info(
        "Picard MarkDuplicates (remove_duplicates=%s): %s → %s",
        remove_duplicates, input_bam, output_bam,
    )
    _run_markdup_cmd(cmd, "Picard MarkDuplicates", timeout)
    if not Path(output_bam).exists():
        raise RuntimeError(f"Picard MarkDuplicates finished but output BAM not found: {output_bam}")
    logger.info("Picard MarkDuplicates complete: %s  metrics: %s", output_bam, metrics_file)
    return output_bam


def _markdup_samtools(
    samtools_exe: str,
    input_bam: str,
    output_bam: str,
    metrics_file: str,
    remove_duplicates: bool,
    threads: int,
    tmp_prefix: str,
    timeout: int,
) -> str:
    """
    Duplicate marking with samtools markdup.

    Uses -m t (template start/end position-based detection, default) which
    considers both read positions in a pair — equivalent to Picard's coordinate
    mode and more accurate than sequence-start-only (-m s).
    -S also marks supplementary alignments of detected duplicates.
    """
    cmd = [
        samtools_exe, "markdup",
        "-@", str(threads),
        "-m", "t",
        "-S",
        "-f", metrics_file,
        "-T", tmp_prefix,
    ]
    if remove_duplicates:
        cmd.append("-r")
    cmd += [input_bam, output_bam]

    logger.info(
        "samtools markdup (threads=%d, remove_duplicates=%s, mode=template): %s → %s",
        threads, remove_duplicates, input_bam, output_bam,
    )
    _run_markdup_cmd(cmd, "samtools markdup", timeout)
    if not Path(output_bam).exists():
        raise RuntimeError(f"samtools markdup finished but output BAM not found: {output_bam}")
    logger.info("samtools markdup complete: %s  metrics: %s", output_bam, metrics_file)
    return output_bam


def _run_markdup_cmd(cmd: List[str], label: str, timeout: int) -> None:
    """Run a markdup subprocess, raising RuntimeError on non-zero exit."""
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout if timeout else None,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"{label} executable not found: {exc}") from exc
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"{label} timed out after {timeout // 3600} hours.\n"
            f"  Increase deepvariant.markdup_timeout_hours in config/settings.yaml, "
            f"or set to 0 to disable."
        )
    if result.returncode != 0:
        tail = (result.stderr or result.stdout).strip().splitlines()[-15:]
        raise RuntimeError(
            f"{label} failed (exit code {result.returncode}).\n"
            f"  Last output:\n    " + "\n    ".join(tail)
        )


def _run_seqkit_pair(
    seqkit_exe: str,
    fastq_r1: str,
    fastq_r2: str,
    out_r1: str,
    out_r2: str,
    timeout: int = 0,
) -> tuple:
    """
    Synchronize paired-end FASTQs with seqkit pair, matching reads by name.

    Handles both mid-file orphaned reads AND files with different total read
    counts (e.g. truncated downloads).  Unpaired reads are silently discarded.

    seqkit writes output files into a temp directory using the input basenames,
    which are then moved to out_r1 / out_r2.

    Returns (out_r1, out_r2) paths on success.
    """
    tmp_dir = Path(out_r1).parent / "_seqkit_pair_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        seqkit_exe, "pair",
        "--read1", fastq_r1,
        "--read2", fastq_r2,
        "--out-dir", str(tmp_dir),
        "--force",
    ]
    logger.info(
        "seqkit pair synchronisation (name-matched): %s + %s -> %s + %s",
        fastq_r1, fastq_r2, out_r1, out_r2,
    )
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout if timeout else None,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"seqkit executable not found at runtime: {seqkit_exe}\n  {exc}") from exc

    if result.returncode != 0:
        tail = result.stderr.strip().splitlines()[-15:]
        raise RuntimeError(
            f"seqkit pair synchronisation failed (exit code {result.returncode}).\n"
            f"  Last output:\n    " + "\n    ".join(tail)
        )

    # Move seqkit outputs (named after input basenames) to the desired paths
    r1_tmp = tmp_dir / Path(fastq_r1).name
    r2_tmp = tmp_dir / Path(fastq_r2).name
    shutil.move(str(r1_tmp), out_r1)
    shutil.move(str(r2_tmp), out_r2)
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    logger.info("seqkit pair synchronisation complete: %s + %s", out_r1, out_r2)
    return out_r1, out_r2


def _run_fastp_repair(
    fastp_exe: str,
    fastq_r1: str,
    fastq_r2: str,
    out_r1: str,
    out_r2: str,
    threads: int,
    timeout: int = _FASTP_TIMEOUT,
) -> tuple:
    """
    Synchronize paired-end FASTQs with fastp, discarding orphaned reads.

    NOTE: fastp requires R1 and R2 to have equal total read counts.  If the
    files have different counts (e.g. truncated download), use seqkit pair
    instead.  This function is kept as a fallback when seqkit is unavailable.

    Returns (out_r1, out_r2) paths on success.
    """
    cmd = [
        fastp_exe,
        "--in1",  fastq_r1,
        "--in2",  fastq_r2,
        "--out1", out_r1,
        "--out2", out_r2,
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--disable_adapter_trimming",
        "--thread", str(min(16, threads)),  # fastp caps at 16
        "--json", "/dev/null",
        "--html", "/dev/null",
    ]
    logger.info(
        "fastp pair synchronisation (threads=%d): %s + %s -> %s + %s",
        min(16, threads), fastq_r1, fastq_r2, out_r1, out_r2,
    )
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout if timeout else None,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"fastp executable not found at runtime: {fastp_exe}\n  {exc}") from exc

    if result.returncode != 0:
        tail = result.stderr.strip().splitlines()[-15:]
        raise RuntimeError(
            f"fastp pair synchronisation failed (exit code {result.returncode}).\n"
            f"  Last output:\n    " + "\n    ".join(tail)
        )
    logger.info("fastp synchronisation complete: %s  %s", out_r1, out_r2)
    return out_r1, out_r2


def _run_bwa_align(
    bwa_exe: str,
    samtools_exe: str,
    bwa_index: str,
    fastq_r1: str,
    fastq_r2: Optional[str],
    sample_id: str,
    sorted_bam: str,
    threads: int,
    timeout: int = _ALIGN_TIMEOUT,
) -> str:
    """
    Align FASTQ reads with BWA MEM and pipe directly through samtools sort.

    Command (paired-end):
        bwa mem -t {threads} -R "@RG\\tID:{id}\\tSM:{id}\\tPL:ILLUMINA" \\
            {index} R1.fastq.gz R2.fastq.gz | samtools sort -@ {threads} -o {sorted_bam}

    Returns the path to the sorted BAM.
    """
    rg_tag = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA\\tLB:{sample_id}"

    # Use pigz for parallel decompression of gzipped FASTQs — BWA's built-in
    # gzip reader is single-threaded and starves alignment threads on large files.
    pigz_exe = shutil.which("pigz")

    def _is_gz(path: str) -> bool:
        return path.endswith(".gz") or path.endswith(".gzip")

    if pigz_exe and (_is_gz(fastq_r1) or (fastq_r2 and _is_gz(fastq_r2))):
        # Decompress threads: quarter of alignment threads, min 8, max 32 per file.
        # Two pigz processes run in parallel (R1 + R2), so total decomp = 2× this.
        decomp_threads = max(8, min(32, threads // 4))
        if fastq_r2 and _is_gz(fastq_r1) and _is_gz(fastq_r2):
            r1_input = f"<(pigz -dc -p {decomp_threads} {shlex.quote(fastq_r1)})"
            r2_input = f"<(pigz -dc -p {decomp_threads} {shlex.quote(fastq_r2)})"
        elif _is_gz(fastq_r1):
            r1_input = f"<(pigz -dc -p {decomp_threads} {shlex.quote(fastq_r1)})"
            r2_input = shlex.quote(fastq_r2) if fastq_r2 else None
        else:
            r1_input = shlex.quote(fastq_r1)
            r2_input = f"<(pigz -dc -p {decomp_threads} {shlex.quote(fastq_r2)})"
        inputs = [r1_input] + ([r2_input] if r2_input else [])
        # fixmate -m adds the ms (mate score) tag required by samtools markdup.
        # BWA MEM outputs reads in query-name order, which is what fixmate needs.
        bwa_shell_cmd = (
            f"set -o pipefail; "
            f"{shlex.quote(bwa_exe)} mem -t {threads} -R {shlex.quote(rg_tag)} "
            f"{shlex.quote(bwa_index)} {' '.join(inputs)} | "
            f"{shlex.quote(samtools_exe)} fixmate -m -O bam - - | "
            f"{shlex.quote(samtools_exe)} sort -@ {threads} -o {shlex.quote(sorted_bam)}"
        )
        use_shell = True
        logger.info(
            "BWA alignment with pigz decompression (bwa_threads=%d, decomp_threads=%d, timeout=%s)",
            threads, decomp_threads,
            f"{timeout // 3600}h" if timeout else "disabled",
        )
    else:
        bwa_cmd: List[str] = [
            bwa_exe, "mem",
            "-t", str(threads),
            "-R", rg_tag,
            bwa_index,
            fastq_r1,
        ]
        if fastq_r2:
            bwa_cmd.append(fastq_r2)
        bwa_shell_cmd = None
        use_shell = False
        logger.info(
            "BWA alignment (threads=%d, timeout=%s): pigz not found — using single-threaded decompression",
            threads,
            f"{timeout // 3600}h" if timeout else "disabled",
        )

    fixmate_cmd: List[str] = [samtools_exe, "fixmate", "-m", "-O", "bam", "-", "-"]
    sort_cmd: List[str] = [
        samtools_exe, "sort",
        "-@", str(threads),
        "-o", sorted_bam,
    ]

    try:
        if use_shell:
            sort_result = subprocess.run(
                ["bash", "-c", bwa_shell_cmd],
                capture_output=True,
                text=True,
                timeout=timeout if timeout else None,
            )
            bwa_stderr = ""
            bwa_returncode = sort_result.returncode
        else:
            bwa_proc = subprocess.Popen(
                bwa_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            # Drain BWA stderr in a background thread so the OS pipe buffer
            # never fills and deadlocks bwa_proc.wait() on long runs.
            _bwa_stderr_buf: list = []
            _stderr_thread = threading.Thread(
                target=lambda: _bwa_stderr_buf.append(bwa_proc.stderr.read()),
                daemon=True,
            )
            _stderr_thread.start()
            fixmate_proc = subprocess.Popen(
                fixmate_cmd,
                stdin=bwa_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            bwa_proc.stdout.close()
            sort_result = subprocess.run(
                sort_cmd,
                stdin=fixmate_proc.stdout,
                capture_output=True,
                text=True,
                timeout=timeout if timeout else None,
            )
            fixmate_proc.stdout.close()
            fixmate_proc.wait(timeout=30)
            bwa_proc.wait(timeout=30)
            _stderr_thread.join()
            bwa_stderr = _bwa_stderr_buf[0].decode("utf-8", errors="replace") if _bwa_stderr_buf else ""
            bwa_returncode = bwa_proc.returncode
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"BWA executable not found at runtime: {bwa_exe}\n  {exc}"
        ) from exc
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"BWA alignment timed out after {timeout // 3600} hours.\n"
            f"  Increase bwa_timeout_hours in config/settings.yaml, or set it to 0 to disable.\n"
            f"  Also consider increasing bwa_threads (currently {threads}) to speed up alignment."
        )

    if bwa_returncode != 0:
        tail = (sort_result.stderr if use_shell else bwa_stderr).strip().splitlines()[-15:]
        raise RuntimeError(
            f"BWA alignment failed (exit code {bwa_returncode}).\n"
            f"  Last output:\n    " + "\n    ".join(tail)
        )
    if sort_result.returncode != 0:
        tail = sort_result.stderr.strip().splitlines()[-10:]
        raise RuntimeError(
            f"samtools sort failed (exit code {sort_result.returncode}).\n"
            f"  Last output:\n    " + "\n    ".join(tail)
        )

    if not Path(sorted_bam).exists():
        raise RuntimeError(
            f"BWA/samtools sort completed but sorted BAM was not produced.\n"
            f"  Expected: {sorted_bam}"
        )

    logger.info("Alignment complete: %s", sorted_bam)
    return sorted_bam


def _bam_is_intact(samtools_exe: str, bam_path: str) -> bool:
    """
    Return True if bam_path passes samtools quickcheck.

    quickcheck verifies the EOF block, catching files truncated mid-write
    by a crash.  Returns False on any failure so the step is re-run.
    """
    try:
        result = subprocess.run(
            [samtools_exe, "quickcheck", bam_path],
            capture_output=True,
        )
        return result.returncode == 0
    except Exception:
        return False


def _vcf_has_records(vcf_path: str) -> bool:
    """
    Return True if vcf_path contains at least one non-header data record.

    Guards against resuming from a VCF written by a process that crashed
    after creating the file and writing headers but before any records.
    """
    try:
        lower = vcf_path.lower()
        open_fn = gzip.open if lower.endswith((".gz", ".bgz")) else open
        with open_fn(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line and not line.startswith("#"):
                    return True
        return False
    except Exception:
        return False


def _run_samtools_index(samtools_exe: str, bam_path: str, threads: int = min(os.cpu_count() or 8, 8)) -> None:
    """Index a sorted BAM with samtools index."""
    cmd = [samtools_exe, "index", "-@", str(threads), bam_path]
    logger.info("Indexing BAM: %s", bam_path)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools index failed (exit code {result.returncode}).\n"
            f"  stderr: {result.stderr.strip()}"
        )
    logger.debug("BAM index written: %s.bai", bam_path)


def _run_flagstat(
    samtools_exe: str,
    bam_path: str,
    out_dir: str,
    sample_id: str,
    threads: int = 8,
) -> dict:
    """
    Run samtools flagstat and write results to {sample_id}.flagstat.txt.

    Parses key QC metrics and logs warnings for problems that would cause
    DeepVariant to waste hours on a bad BAM.  Never raises — failures are
    logged and an empty dict is returned so the pipeline continues.

    Returns a dict with "total_reads" and "mapped_reads" on success.
    """
    out_file = str(Path(out_dir) / f"{sample_id}.flagstat.txt")
    cmd = [samtools_exe, "flagstat", "-@", str(threads), bam_path]

    logger.info("BAM QC — samtools flagstat: %s", bam_path)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        logger.warning("samtools flagstat timed out after 5 minutes — skipping BAM QC")
        return {}
    except Exception as exc:
        logger.warning("samtools flagstat could not run — skipping BAM QC: %s", exc)
        return {}

    if result.returncode != 0:
        logger.warning(
            "samtools flagstat returned exit code %d — skipping BAM QC", result.returncode,
        )
        return {}

    try:
        Path(out_file).write_text(result.stdout)
        logger.info("Flagstat written: %s", out_file)
    except Exception:
        pass

    total = 0
    mapped = 0
    for line in result.stdout.splitlines():
        m_total = re.match(r"^(\d+)\s+\+\s+\d+\s+in total\b", line)
        if m_total:
            total = int(m_total.group(1))
            continue
        # "mapped" line but NOT "primary mapped"
        m_mapped = re.match(r"^(\d+)\s+\+\s+\d+\s+mapped\b", line)
        if m_mapped and "primary" not in line:
            mapped = int(m_mapped.group(1))

    if total == 0:
        logger.error(
            "BAM QC: BAM contains no reads — DeepVariant will produce no output: %s",
            bam_path,
        )
    elif mapped == 0:
        logger.error(
            "BAM QC: 0 reads mapped — DeepVariant will produce no variants: %s",
            bam_path,
        )
    else:
        pct = mapped / total * 100
        if pct < 70.0:
            logger.warning(
                "BAM QC: low mapping rate %.1f%% (%d / %d reads mapped) — "
                "check alignment quality before proceeding with DeepVariant: %s",
                pct, mapped, total, bam_path,
            )
        else:
            logger.info("BAM QC: %.1f%% mapped (%d / %d reads)", pct, mapped, total)

    return {"total_reads": total, "mapped_reads": mapped}


def _detect_gpu() -> bool:
    """Return True if an NVIDIA GPU is present and accessible."""
    if not shutil.which("nvidia-smi"):
        return False
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10,
        )
        return result.returncode == 0 and bool(result.stdout.strip())
    except Exception:
        return False


def _run_deepvariant(
    docker_exe: str,
    fasta_path: str,
    bam_path: str,
    out_dir: str,
    sample_id: str,
    num_shards: int,
    model_type: str,
    version: str,
    timeout: int,
    visual_report: bool = False,
    use_gpu: bool = False,
    disable_small_model: bool = False,
) -> str:
    """
    Call variants with DeepVariant running inside the official Docker image.

    Volume layout inside the container:
      /workdir  <-- out_dir (rw): contains the sorted BAM; output VCF written here
      /ref      <-- fasta parent directory (ro): FASTA + .fai index
                    (omitted when fasta lives inside out_dir -- /workdir is used instead)

    When use_gpu=True the GPU-enabled image tag (<version>-gpu) is used and
    --gpus all is passed to Docker.  call_variants runs on the GPU; make_examples
    is always CPU-only regardless.

    The output file on the host is:
      {out_dir}/{sample_id}.dv.vcf.gz  (bgzipped VCF)

    Returns the host path to the output VCF.GZ.
    """
    fasta  = Path(fasta_path).resolve()
    out    = Path(out_dir).resolve()
    bam    = Path(bam_path).resolve()
    out_vcf_name = f"{sample_id}.dv.vcf.gz"

    # --- Build Docker volume mounts ---
    # /workdir — output directory (receives the VCF; also contains the BAM when
    #            called from the FASTQ path where the BAM was written there).
    # /ref     — FASTA parent directory (added when distinct from /workdir).
    # /bam     — BAM parent directory (added when the BAM is not inside /workdir
    #            or /ref, e.g. when a pre-existing BAM from elsewhere is supplied).
    fasta_dir = fasta.parent
    bam_dir   = bam.parent

    vol_args: List[str] = ["-v", f"{out}:/workdir"]

    if fasta_dir == out:
        ref_in_container = f"/workdir/{fasta.name}"
    else:
        vol_args += ["-v", f"{fasta_dir}:/ref:ro"]
        ref_in_container = f"/ref/{fasta.name}"

    if bam_dir == out:
        bam_in_container = f"/workdir/{bam.name}"
    elif bam_dir == fasta_dir:
        bam_in_container = f"/ref/{bam.name}"
    else:
        vol_args += ["-v", f"{bam_dir}:/bam:ro"]
        bam_in_container = f"/bam/{bam.name}"

    out_vcf_in_container = f"/workdir/{out_vcf_name}"

    # Run as the current user so output files are owned correctly.
    uid_gid = f"{os.getuid()}:{os.getgid()}"

    image_tag = f"{version}-gpu" if use_gpu else version
    gpu_args: List[str] = ["--gpus", "all"] if use_gpu else []

    cmd: List[str] = [
        docker_exe, "run", "--rm",
        "--user", uid_gid,
    ] + gpu_args + vol_args + [
        f"{_DV_IMAGE}:{image_tag}",
        "/opt/deepvariant/bin/run_deepvariant",
        f"--model_type={model_type}",
        f"--ref={ref_in_container}",
        f"--reads={bam_in_container}",
        f"--output_vcf={out_vcf_in_container}",
        f"--num_shards={num_shards}",
    ]
    if not visual_report:
        cmd.append("--novcf_stats_report")
    if disable_small_model:
        cmd.append("--disable_small_model")

    logger.info(
        "DeepVariant variant calling (%s):\n  %s",
        "GPU" if use_gpu else "CPU",
        " ".join(shlex.quote(t) for t in cmd),
    )
    logger.info(
        "DeepVariant is running -- this typically takes 1-8 hours for WGS (less with GPU). "
        "If the image is not cached locally, Docker will pull it first (~1 GB)."
    )

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout if timeout else None,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"Docker executable not found: {docker_exe}\n"
            f"  Ensure Docker is installed and on your PATH, then re-run setup:\n"
            f"    python cli.py --setup"
        )
    except subprocess.TimeoutExpired:
        hours = timeout / 3600
        raise RuntimeError(
            f"DeepVariant timed out after {hours:.1f} hour(s).\n"
            f"\n"
            f"  To fix, increase the timeout in config/settings.yaml:\n"
            f"    deepvariant:\n"
            f"      timeout_hours: {int(hours) + 6}\n"
            f"\n"
            f"  Or increase num_shards to use more CPU cores:\n"
            f"    deepvariant:\n"
            f"      num_shards: {num_shards + 4}"
        )

    if result.returncode != 0:
        stderr_lines = result.stderr.strip().splitlines()
        tail = "\n    ".join(stderr_lines[-20:]) if stderr_lines else "(no output)"
        raise RuntimeError(
            f"DeepVariant exited with code {result.returncode}.\n"
            f"\n"
            f"  Last output:\n    {tail}\n"
            f"\n"
            f"  Common causes:\n"
            f"    - Docker image could not be pulled (check internet connection)\n"
            f"    - Insufficient memory (WGS requires ~16 GB RAM)\n"
            f"    - Reference FASTA index (.fai) missing -- run: samtools faidx {fasta_path}\n"
            f"    - BAM is not coordinate-sorted or index (.bai) is missing\n"
            f"    - Chromosome names in BAM do not match the reference (e.g. '1' vs 'chr1')"
        )

    actual_output = str(out / out_vcf_name)
    if not Path(actual_output).exists():
        raise RuntimeError(
            f"DeepVariant exited successfully but output VCF was not produced.\n"
            f"  Expected: {actual_output}\n"
            f"  Check that {out} is writable and has sufficient disk space."
        )

    logger.info("DeepVariant variant calling complete: %s", actual_output)
    return actual_output


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_fastq_to_vcf(
    fastq_r1: str,
    fastq_r2: Optional[str],
    sample_id: str,
    output_dir: str,
    paths_config: dict,
    settings: dict,
    stage_callbacks=None,
) -> str:
    """
    Run the FASTQ -> sorted BAM -> VCF pipeline using BWA + DeepVariant.

    Parameters
    ----------
    fastq_r1:     Path to the R1 (or only) FASTQ file (.fastq / .fastq.gz).
    fastq_r2:     Path to the R2 FASTQ file for paired-end input, or None.
    sample_id:    Sample identifier used for read-group tags and output naming.
    output_dir:   Directory where all intermediate and final files are written.
    paths_config: Parsed config/paths.yaml.
    settings:     Parsed config/settings.yaml.

    Returns
    -------
    Path to the raw (unannotated) VCF.GZ file produced by DeepVariant.

    Raises
    ------
    RuntimeError: If any required tool or resource is missing, or if any
                  pipeline step exits non-zero.
    """
    # --- Validate and locate tools ---
    bwa_exe      = _require_tool("bwa")
    samtools_exe = _require_tool("samtools")
    docker_exe   = _get_docker(paths_config)
    fasta_path   = _get_reference_fasta(paths_config)
    bwa_index    = _get_bwa_index_prefix(paths_config, fasta_path)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # DeepVariant settings (with fallbacks)
    dv_cfg     = settings.get("deepvariant", {})
    dv_threads = int(dv_cfg.get("num_shards",
                                settings.get("vep", {}).get("fork", _DV_SHARDS_DEFAULT)))
    version    = str(dv_cfg.get("version",    _DV_VERSION_DEFAULT))
    model_type = str(dv_cfg.get("model_type", _DV_MODEL_DEFAULT))
    timeout_h  = float(dv_cfg.get("timeout_hours", _DV_TIMEOUT_DEFAULT / 3600))
    timeout    = int(timeout_h * 3600) if timeout_h > 0 else None

    gpu_setting = str(dv_cfg.get("use_gpu", _DV_USE_GPU_DEFAULT)).lower()
    if gpu_setting == "auto":
        use_gpu = _detect_gpu()
        logger.info("GPU auto-detection: %s", "GPU found — using GPU image" if use_gpu else "no GPU — using CPU")
    else:
        use_gpu = gpu_setting in ("true", "1", "yes")

    disable_small_model = bool(dv_cfg.get("disable_small_model", False))

    # BWA-specific settings — bwa_threads defaults to dv_threads if not set separately.
    # bwa_timeout_hours: 0 = no timeout.
    bwa_threads   = int(dv_cfg.get("bwa_threads", dv_threads))
    bwa_timeout_h = float(dv_cfg.get("bwa_timeout_hours", _ALIGN_TIMEOUT / 3600))
    bwa_timeout   = int(bwa_timeout_h * 3600) if bwa_timeout_h > 0 else 0

    mode = "paired-end" if fastq_r2 else "single-end"
    logger.info(
        "FASTQ -> VCF  sample=%s  mode=%s  caller=DeepVariant-%s  shards=%d  "
        "bwa_threads=%d  bwa_timeout=%s",
        sample_id, mode, version, dv_threads, bwa_threads,
        f"{bwa_timeout // 3600}h" if bwa_timeout else "disabled",
    )

    # --- Step 0: Synchronize paired-end FASTQs ---
    # SRA/ENA datasets often have R1/R2 files that are desynced — reads missing
    # from one file shift all subsequent read names out of lock-step with BWA.
    # Some datasets also have different total read counts (truncated downloads).
    # seqkit pair matches reads by name and discards unpaired reads, handling
    # both cases.  fastp is kept as a fallback but only works when total read
    # counts are equal.
    if fastq_r2:
        seqkit_exe = shutil.which("seqkit")
        fastp_exe  = shutil.which("fastp")
        sync_tool  = seqkit_exe or fastp_exe

        if sync_tool:
            r1_synced = str(out / f"{sample_id}.r1.synced.fastq.gz")
            r2_synced = str(out / f"{sample_id}.r2.synced.fastq.gz")
            if (Path(r1_synced).exists() and Path(r2_synced).exists()
                    and Path(r1_synced).stat().st_size > 0
                    and Path(r2_synced).stat().st_size > 0):
                logger.info(
                    "Resuming: synced FASTQs already exist — skipping sync: %s / %s",
                    r1_synced, r2_synced,
                )
                fastq_r1 = r1_synced
                fastq_r2 = r2_synced
            elif seqkit_exe:
                fastq_r1, fastq_r2 = _run_seqkit_pair(
                    seqkit_exe=seqkit_exe,
                    fastq_r1=fastq_r1,
                    fastq_r2=fastq_r2,
                    out_r1=r1_synced,
                    out_r2=r2_synced,
                )
            else:
                logger.warning(
                    "seqkit not found — using fastp for sync (requires equal read counts). "
                    "Install seqkit for full support:  conda install -c bioconda seqkit"
                )
                fastq_r1, fastq_r2 = _run_fastp_repair(
                    fastp_exe=fastp_exe,
                    fastq_r1=fastq_r1,
                    fastq_r2=fastq_r2,
                    out_r1=r1_synced,
                    out_r2=r2_synced,
                    threads=bwa_threads,
                )
        else:
            logger.warning(
                "Neither seqkit nor fastp found — skipping FASTQ pair synchronisation. "
                "BWA may fail with 'paired reads have different names' if reads are desynced. "
                "Install seqkit:  conda install -c bioconda seqkit"
            )

    def _fire(idx: int) -> None:
        """Fire stage_callbacks[idx] if it exists."""
        if stage_callbacks and idx < len(stage_callbacks) and stage_callbacks[idx]:
            stage_callbacks[idx]()

    # --- Step 1: Align and sort ---
    _fire(0)  # → "Aligning reads  (BWA)"
    sorted_bam = str(out / f"{sample_id}.sorted.bam")
    bam_index  = sorted_bam + ".bai"
    if (Path(sorted_bam).exists() and Path(bam_index).exists()
            and _bam_is_intact(samtools_exe, sorted_bam)):
        logger.info("Resuming: sorted BAM already exists — skipping alignment: %s", sorted_bam)
    else:
        _run_bwa_align(
            bwa_exe=bwa_exe,
            samtools_exe=samtools_exe,
            bwa_index=bwa_index,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            sample_id=sample_id,
            sorted_bam=sorted_bam,
            threads=bwa_threads,
            timeout=bwa_timeout,
        )
        _run_samtools_index(samtools_exe, sorted_bam, threads=bwa_threads)

    # --- Step 2: Mark and remove duplicates ---
    # Uses coordinate-sorted (position-based) duplicate detection.
    # Tool auto-selected: sambamba (fastest) → picard → samtools (always available).
    _fire(1)  # → "Marking duplicates  (sambamba)"
    markdup_tool    = str(dv_cfg.get("markdup_tool", "sambamba")).strip()
    remove_dups     = bool(dv_cfg.get("remove_duplicates", True))
    markdup_timeout_h = float(dv_cfg.get("markdup_timeout_hours", _MARKDUP_TIMEOUT / 3600))
    markdup_timeout = int(markdup_timeout_h * 3600) if markdup_timeout_h > 0 else 0
    # sambamba markdup creates one temp file per thread — cap at 32 to avoid
    # hitting OS open-file limits on very-high-core machines.  Markdup is
    # I/O-bound so gains nothing above ~32 threads.
    markdup_threads = min(32, bwa_threads)

    if markdup_tool:
        dedup_bam = str(out / f"{sample_id}.dedup.bam")
        dedup_bai = dedup_bam + ".bai"

        if (Path(dedup_bam).exists() and Path(dedup_bai).exists()
                and _bam_is_intact(samtools_exe, dedup_bam)):
            logger.info(
                "Resuming: dedup BAM already exists — skipping duplicate marking: %s",
                dedup_bam,
            )
            bam_for_calling = dedup_bam
        else:
            metrics_file = str(out / f"{sample_id}.markdup.metrics.txt")
            tmp_prefix   = str(out / f"{sample_id}.markdup.tmp")
            bam_for_calling = _run_markdup(
                samtools_exe=samtools_exe,
                input_bam=sorted_bam,
                output_bam=dedup_bam,
                metrics_file=metrics_file,
                remove_duplicates=remove_dups,
                preferred_tool=markdup_tool,
                threads=markdup_threads,
                tmp_prefix=tmp_prefix,
                timeout=markdup_timeout,
            )
            _run_samtools_index(samtools_exe, bam_for_calling, threads=bwa_threads)
    else:
        logger.info("Duplicate marking disabled (markdup_tool is empty) — using sorted BAM directly.")
        bam_for_calling = sorted_bam

    _fire(2)  # → "Calling variants  (DeepVariant)"

    # --- BAM QC: flagstat before the expensive DeepVariant step ---
    _run_flagstat(samtools_exe, bam_for_calling, str(out), sample_id, threads=bwa_threads)

    # --- Step 3: Call variants with DeepVariant ---
    expected_vcf = str(out / f"{sample_id}.dv.vcf.gz")
    if Path(expected_vcf).exists() and _vcf_has_records(expected_vcf):
        logger.info("Resuming: DeepVariant VCF already exists — skipping variant calling: %s", expected_vcf)
        raw_vcf = expected_vcf
    else:
        raw_vcf = _run_deepvariant(
            docker_exe=docker_exe,
            fasta_path=fasta_path,
            bam_path=bam_for_calling,
            out_dir=str(out),
            sample_id=sample_id,
            num_shards=dv_threads,
            model_type=model_type,
            version=version,
            timeout=timeout,
            visual_report=bool(dv_cfg.get("visual_report", False)),
            use_gpu=use_gpu,
            disable_small_model=disable_small_model,
        )

    logger.info(
        "FASTQ processing complete: %s -> %s -> %s -> %s",
        fastq_r1, sorted_bam, bam_for_calling, raw_vcf,
    )
    return raw_vcf
