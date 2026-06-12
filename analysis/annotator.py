"""
annotator.py — Build and execute the VEP CLI annotation command.

Version 1 design decisions:
  - Command-line VEP only. No VEP web API is used or contacted.
  - Genome assembly: GRCh38 only.
  - Cache mode: --offline --cache. All annotation is local.
  - Output format: plain VCF (--vcf), suitable for parser.py.
  - Required resources: VEP executable, VEP cache dir, reference FASTA.
  - Optional plugins: loaded only when the resource file is present and
    confirmed in paths.yaml. Missing optional resources are logged clearly
    as warnings, not errors.
  - ClinVar: loaded via VEP --custom if configured.
  - gnomAD AF: sourced from the VEP merged cache (--af_gnomad / --max_af),
    not from a separate gnomAD VCF, to avoid requiring the user to download
    a multi-hundred-GB file just for AF lookups.

Output consumed by:
  analysis/parser.py — parses the CSQ INFO field from the annotated VCF.
"""

import gzip
import logging
import os
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths-config accessors
# ---------------------------------------------------------------------------

def _get_required_path(
    paths_config: dict,
    section: str,
    key: str,
    label: str,
) -> str:
    """
    Extract a required validated path from paths_config.

    Raises RuntimeError with a beginner-friendly message if the path is
    absent, null, or has any status other than "found".
    """
    entry = paths_config.get(section, {}).get(key, {})
    if not isinstance(entry, dict):
        _missing_required(label, section, key, status="not_in_config")

    status = entry.get("status", "not_checked")
    path = entry.get("path")

    if status != "found" or not path:
        _missing_required(label, section, key, status=status)

    return path


def _missing_required(label: str, section: str, key: str, status: str) -> None:
    raise RuntimeError(
        f"{label} is required for annotation but is not available.\n"
        f"  Config entry : {section}.{key}\n"
        f"  Status       : {status}\n"
        f"\n"
        f"  To fix:\n"
        f"    1. Run setup:  python cli.py --setup\n"
        f"    2. If the path was not detected automatically, add it to\n"
        f"       config/paths.yaml under {section}.{key}.path\n"
        f"    3. Re-run setup to validate the path."
    )


def _get_optional_path(paths_config: dict, section: str, key: str) -> Optional[str]:
    """
    Extract an optional resource path from paths_config.

    Returns the path if status == "found", None otherwise.
    Logs a WARNING when a path is set but status is not "found" (e.g. the
    file was configured but has since gone missing or failed validation).
    Never raises.
    """
    entry = paths_config.get(section, {}).get(key, {})
    if not isinstance(entry, dict):
        return None
    status = entry.get("status", "not_checked")
    path = entry.get("path")
    if path and status != "found":
        logger.warning(
            "Resource %s.%s is configured (path: %s) but has status '%s' — "
            "it will be skipped. Re-run setup to re-validate: "
            "python cli.py --setup --force",
            section, key, path, status,
        )
    if path and status == "found":
        return path
    return None


# ---------------------------------------------------------------------------
# Plugin and custom-annotation builders
# ---------------------------------------------------------------------------

def _build_plugin_args(paths_config: dict, settings: dict) -> List[str]:
    """
    Build VEP --plugin arguments for all optional annotation resources.

    Each plugin is only added when its required resource file(s) are present
    in paths_config with status "found".  Missing plugins are logged as INFO
    (not configured) or WARNING (partially configured — e.g. one of a pair).

    Returns a flat list of argument tokens to append to the VEP command.
    """
    args: List[str] = []

    # --- CADD (requires both SNV and indel files) ---
    cadd_snv = _get_optional_path(paths_config, "resources", "cadd_snv_file")
    cadd_indel = _get_optional_path(paths_config, "resources", "cadd_indel_file")
    if cadd_snv and cadd_indel:
        args += ["--plugin", f"CADD,{cadd_snv},{cadd_indel}"]
        logger.info("Plugin enabled: CADD (SNV=%s, indel=%s)", cadd_snv, cadd_indel)
    elif cadd_snv or cadd_indel:
        present = cadd_snv or cadd_indel
        missing_key = "cadd_indel_file" if cadd_snv else "cadd_snv_file"
        logger.warning(
            "CADD plugin skipped: both SNV and indel score files are required. "
            "Found: %s. Missing: resources.%s in config/paths.yaml.",
            present, missing_key,
        )
    else:
        logger.info(
            "Plugin not configured: CADD — skipped "
            "(set resources.cadd_snv_file + resources.cadd_indel_file in config/paths.yaml)"
        )

    # --- REVEL ---
    revel = _get_optional_path(paths_config, "resources", "revel_file")
    if revel:
        args += ["--plugin", f"REVEL,{revel}"]
        logger.info("Plugin enabled:        REVEL  %s", revel)
    else:
        logger.info(
            "Plugin not configured: REVEL — skipped "
            "(set resources.revel_file in config/paths.yaml to enable)"
        )

    # --- SpliceAI ---
    # Both SNV and indel files give full coverage.
    # SNV-only mode is supported when the indel file is unavailable
    # (Illumina's indel file requires BaseSpace authentication; the SNV file
    # is publicly available from the Ensembl FTP).
    spliceai_snv   = _get_optional_path(paths_config, "resources", "spliceai_snv_file")
    spliceai_indel = _get_optional_path(paths_config, "resources", "spliceai_indel_file")
    if spliceai_snv and spliceai_indel:
        args += ["--plugin", f"SpliceAI,snv={spliceai_snv},indel={spliceai_indel}"]
        logger.info("Plugin enabled:        SpliceAI  snv=%s  indel=%s",
                    spliceai_snv, spliceai_indel)
    elif spliceai_snv:
        # SNV-only mode: covers substitutions and most splice-site variants.
        # In-frame indel splice impact will not be scored.
        args += ["--plugin", f"SpliceAI,snv={spliceai_snv}"]
        logger.warning(
            "Plugin enabled (SNV-only): SpliceAI  snv=%s  "
            "indel file not configured — in-frame indel splice impact will not be scored. "
            "To enable full coverage, set resources.spliceai_indel_file in config/paths.yaml.",
            spliceai_snv,
        )
    elif spliceai_indel:
        # Indel-only: unusual but allow it rather than silently skipping.
        args += ["--plugin", f"SpliceAI,indel={spliceai_indel}"]
        logger.warning(
            "Plugin enabled (indel-only): SpliceAI  indel=%s  "
            "SNV file not configured. Set resources.spliceai_snv_file in config/paths.yaml.",
            spliceai_indel,
        )
    else:
        logger.info(
            "Plugin not configured: SpliceAI — skipped "
            "(set resources.spliceai_snv_file in config/paths.yaml to enable; "
            "spliceai_indel_file is optional)"
        )

    # --- AlphaMissense ---
    alphamissense = _get_optional_path(paths_config, "resources", "alphamissense_file")
    if alphamissense:
        args += ["--plugin", f"AlphaMissense,file={alphamissense}"]
        logger.info("Plugin enabled:        AlphaMissense  %s", alphamissense)
    else:
        logger.info(
            "Plugin not configured: AlphaMissense — skipped "
            "(set resources.alphamissense_file in config/paths.yaml to enable)"
        )

    # --- LOFTEE ---
    # loftee_path points to the cloned LOFTEE source directory (contains
    # LoF.pm + helper scripts like gerp_dist.pl, splice_data/, maxEntScan/).
    # This is separate from vep_plugins_dir (which holds all other .pm files).
    loftee_dir      = _get_optional_path(paths_config, "resources", "loftee_path")
    loftee_human    = _get_optional_path(paths_config, "resources", "loftee_human_ancestor_fa")
    loftee_sql      = _get_optional_path(paths_config, "resources", "loftee_conservation_file")
    loftee_gerp_bw  = _get_optional_path(paths_config, "resources", "loftee_gerp_bigwig")

    if loftee_dir:
        plugin_parts = [f"LoF", f"loftee_path:{loftee_dir}"]

        if loftee_human:
            plugin_parts.append(f"human_ancestor_fa:{loftee_human}")
        else:
            logger.warning(
                "LOFTEE: human_ancestor_fa not configured — ancestral allele filter "
                "will be skipped. Set resources.loftee_human_ancestor_fa in "
                "config/paths.yaml to enable."
            )

        if loftee_sql:
            plugin_parts.append(f"conservation_file:{loftee_sql}")
        else:
            logger.warning(
                "LOFTEE: conservation_file (loftee.sql) not configured — GERP "
                "conservation filter will be skipped. Set resources.loftee_conservation_file "
                "in config/paths.yaml to enable."
            )

        if loftee_gerp_bw:
            plugin_parts.append(f"gerp_bigwig:{loftee_gerp_bw}")
        else:
            logger.info(
                "LOFTEE: gerp_bigwig not configured — BigWig conservation filter "
                "will be skipped. Set resources.loftee_gerp_bigwig in config/paths.yaml."
            )

        args += ["--plugin", ",".join(plugin_parts)]
        logger.info(
            "Plugin enabled:        LOFTEE  loftee_path=%s  human_ancestor=%s  "
            "conservation_file=%s  gerp_bigwig=%s",
            loftee_dir,
            loftee_human or "not configured",
            loftee_sql or "not configured",
            loftee_gerp_bw or "not configured",
        )
    else:
        logger.info(
            "Plugin not configured: LOFTEE — skipped "
            "(set resources.vep_plugins_dir in config/paths.yaml; "
            "LOFTEE must be installed in that directory)"
        )

    return args


def _build_custom_args(paths_config: dict) -> List[str]:
    """
    Build VEP --custom arguments for file-based annotation sources.

    Currently: ClinVar only.
    gnomAD is sourced from the VEP merged cache (--af_gnomad), not a
    separate VCF, to avoid requiring a multi-hundred-GB download.

    Returns a flat list of argument tokens.
    """
    args: List[str] = []

    clinvar = _get_optional_path(paths_config, "resources", "clinvar_vcf")
    if clinvar:
        # ClinVar fields extracted:
        #   CLNSIG      — clinical significance (Pathogenic, Benign, etc.)
        #   CLNREVSTAT  — review status (review star count proxy)
        #   CLNDN       — disease name(s)
        args += [
            "--custom",
            f"{clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN",
        ]
        logger.info("Custom annotation enabled:        ClinVar  %s", clinvar)
    else:
        logger.info(
            "Custom annotation not configured: ClinVar — skipped "
            "(set resources.clinvar_vcf in config/paths.yaml to enable)"
        )

    return args


# ---------------------------------------------------------------------------
# Command builder
# ---------------------------------------------------------------------------

def build_vep_command(
    input_vcf: str,
    output_vcf: str,
    paths_config: dict,
    settings: dict,
) -> List[str]:
    """
    Build the complete VEP CLI command as a list of argument tokens.

    Does not touch the filesystem beyond validating paths in paths_config.
    Raises RuntimeError if any required resource is missing or misconfigured.

    Parameters
    ----------
    input_vcf:    Validated input VCF path (from vcf_validator).
    output_vcf:   Desired output path.  Plain VCF in version 1.
    paths_config: Parsed config/paths.yaml (from config_writer.load_paths()).
    settings:     Parsed config/settings.yaml.

    Returns
    -------
    List of command token strings suitable for subprocess.run(cmd).
    """
    # --- Required paths ---
    vep_exe = _get_required_path(
        paths_config, "tools", "vep", "VEP executable"
    )
    cache_dir = _get_required_path(
        paths_config, "resources", "vep_cache_dir", "VEP cache directory"
    )
    fasta = _get_required_path(
        paths_config, "resources", "reference_fasta", "GRCh38 reference FASTA"
    )

    # --- Optional paths ---
    plugins_dir = _get_optional_path(paths_config, "resources", "vep_plugins_dir")

    # --- Settings ---
    vep_cfg     = settings.get("vep", {})
    fork        = int(vep_cfg.get("fork", 4))
    buffer_size = int(vep_cfg.get("buffer_size", 20000))
    distance    = int(vep_cfg.get("distance", 5000))

    # --- Base command ---
    # Flags are ordered: I/O → cache → assembly → annotation → output options.
    cmd: List[str] = [
        vep_exe,

        # I/O
        "--input_file",      input_vcf,
        "--output_file",     output_vcf,
        "--vcf",                              # output format: VCF with CSQ INFO field

        # Cache and offline mode (no web traffic)
        "--cache",
        "--offline",
        "--assembly",        "GRCh38",
        "--dir_cache",       cache_dir,
        "--fasta",           fasta,

        # Performance
        "--fork",            str(fork),
        "--buffer_size",     str(buffer_size),

        # Gene / transcript annotation
        "--distance",        str(distance),
        "--symbol",                           # HGNC gene symbol
        "--biotype",                          # transcript biotype (protein_coding etc.)
        "--numbers",                          # exon/intron numbers
        "--canonical",                        # flag canonical transcript
        "--transcript_version",               # include transcript version in output

        # HGVS notation (required for reporting)
        "--hgvs",                             # HGVS c. and p.
        "--hgvsg",                            # HGVS g. (genomic)

        # Population frequency (from cache — no separate gnomAD VCF required)
        "--af",                               # global 1000 Genomes AF
        "--af_gnomad",                        # gnomAD exomes AF per population
        "--max_af",                           # report maximum AF across all populations

        # Pathogenicity predictions (built into VEP — no plugin required)
        "--sift",            "b",             # SIFT score + label (tolerated/deleterious)
        "--polyphen",        "b",             # PolyPhen score + label

        # Functional context
        "--domains",                          # overlapping protein domains
        "--uniprot",                          # UniProt accessions (for domain queries)

        # Output cleanliness
        "--no_stats",                         # suppress HTML statistics report
        "--no_progress",                      # suppress progress bar (cleaner log)
        "--force_overwrite",                  # overwrite output if it exists
    ]

    # --- Optional: plugins directory ---
    # All plugin .pm files and LOFTEE helper scripts live in vep_plugins_dir.
    # The GRCh38 LoF.pm and its helpers are symlinked in from loftee_grch38/.
    if plugins_dir:
        cmd += ["--dir_plugins", plugins_dir]
    else:
        logger.info(
            "vep_plugins_dir not configured — plugins requiring a separate dir "
            "may fail to load. Add the path to config/paths.yaml if using plugins."
        )

    # --- Optional plugin arguments ---
    cmd += _build_plugin_args(paths_config, settings)

    # --- Optional custom annotation arguments ---
    cmd += _build_custom_args(paths_config)

    # --- User-defined extra flags (from settings.yaml) ---
    extra = vep_cfg.get("extra_flags", []) or []
    if extra:
        logger.info("Appending extra_flags from settings.yaml: %s", extra)
        cmd += extra

    return cmd


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def _count_vcf_records(vcf_path: str) -> Optional[int]:
    """
    Count the number of data records in a VCF (plain or gzipped).

    Returns None if the file cannot be read.  Used only for logging — never
    raises so it cannot block annotation.
    """
    try:
        count = 0
        lower = vcf_path.lower()
        open_fn = gzip.open if lower.endswith((".gz", ".bgz")) else open
        with open_fn(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line and not line.startswith("#"):
                    count += 1
        return count
    except Exception:
        return None


def _resolve_timeout(settings: dict) -> Optional[int]:
    """
    Return the VEP timeout in seconds from settings, or None to disable.

    Reads  settings.vep.timeout_hours.
    0 or negative  → None (no timeout).
    Default if absent: 6 hours (21600 seconds).
    """
    hours = settings.get("vep", {}).get("timeout_hours", 6)
    try:
        hours = float(hours)
    except (TypeError, ValueError):
        hours = 6
    if hours <= 0:
        return None
    return int(hours * 3600)


def _output_is_complete(output_vcf: str) -> bool:
    """
    Return True if output_vcf already exists and contains a VEP CSQ header.

    Used to implement skip_if_annotated_exists: the file must not only exist
    but also have a valid CSQ annotation line so a partial/truncated output
    from a previous failed run is not mistakenly accepted.
    """
    if not Path(output_vcf).exists():
        return False
    try:
        lower = output_vcf.lower()
        open_fn = gzip.open if lower.endswith((".gz", ".bgz")) else open
        with open_fn(output_vcf, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                if line.startswith("##INFO=<ID=CSQ"):
                    return True
    except Exception:
        pass
    return False


def run_vep(command: List[str], output_vcf: str, settings: Optional[dict] = None) -> str:
    """
    Execute a pre-built VEP command via subprocess.

    Parameters
    ----------
    command:    Command list from build_vep_command().
    output_vcf: Expected output path; checked for existence after VEP exits.
    settings:   Parsed config/settings.yaml.  Used to read vep.timeout_hours.
                Defaults to {} (6-hour timeout) if not supplied.

    Returns
    -------
    Path to the output VCF.

    Raises
    ------
    RuntimeError: If the VEP executable is not found, times out, exits
                  non-zero, or does not produce the expected output file.
    """
    settings = settings or {}
    timeout_sec = _resolve_timeout(settings)
    timeout_label = (
        f"{timeout_sec // 3600}h {(timeout_sec % 3600) // 60}m"
        if timeout_sec else "disabled"
    )

    cmd_display = " ".join(shlex.quote(t) for t in command)
    logger.info("Executing VEP command:\n  %s", cmd_display)
    logger.info(
        "VEP output will be written to: %s  |  timeout: %s",
        output_vcf, timeout_label,
    )
    logger.info(
        "VEP is running — this can take several hours for a whole-genome VCF. "
        "Progress is written to VEP's own log (--no_progress suppresses the bar). "
        "The output file grows as VEP writes annotated records."
    )

    # Build the subprocess environment.
    # VEP (and LOFTEE) require Perl modules from both the ensembl-vep tree
    # and the user local Perl lib. We extend — not replace — the current env
    # so that PATH, LD_LIBRARY_PATH, etc. are all inherited.
    vep_dir   = str(Path(command[0]).parent)           # e.g. /home/admin7995/ensembl-vep
    local_lib = str(Path.home() / "perl5" / "lib" / "perl5")
    existing  = os.environ.get("PERL5LIB", "")
    parts     = [p for p in [vep_dir, local_lib, existing] if p]
    env = os.environ.copy()
    env["PERL5LIB"] = ":".join(parts)

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_sec,        # None = no limit
            env=env,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"VEP executable not found at runtime: {command[0]}\n"
            f"  The path was set during setup but the file is no longer accessible.\n"
            f"  Re-run setup to re-validate:  python cli.py --setup"
        )
    except subprocess.TimeoutExpired:
        hours_used = (timeout_sec or 0) / 3600
        raise RuntimeError(
            f"VEP annotation timed out after {hours_used:.1f} hour(s).\n"
            f"\n"
            f"  Output file (may be partial): {output_vcf}\n"
            f"\n"
            f"  To fix — increase the timeout in config/settings.yaml:\n"
            f"    vep:\n"
            f"      timeout_hours: {int(hours_used) + 6}   # current: {hours_used:.0f}\n"
            f"\n"
            f"  Other options:\n"
            f"    - Increase vep.fork to use more CPU threads (current: "
            f"{settings.get('vep', {}).get('fork', 4)})\n"
            f"    - Split the VCF by chromosome before annotating:\n"
            f"        bcftools view -r chr1 input.vcf.gz -Oz -o chr1.vcf.gz\n"
            f"        tabix -p vcf chr1.vcf.gz\n"
            f"    - Set timeout_hours: 0 to disable the timeout entirely"
        )

    # Log VEP's own output for the run log
    if result.stdout:
        logger.debug("VEP stdout:\n%s", result.stdout.strip())
    if result.stderr:
        # VEP writes progress/summary to stderr — log at DEBUG, not WARNING
        logger.debug("VEP stderr:\n%s", result.stderr.strip())

    if result.returncode != 0:
        stderr_lines = result.stderr.strip().splitlines()
        stderr_tail = "\n    ".join(stderr_lines[-15:]) if stderr_lines else "(no output)"
        raise RuntimeError(
            f"VEP exited with error code {result.returncode}.\n"
            f"\n"
            f"  VEP output (last lines):\n"
            f"    {stderr_tail}\n"
            f"\n"
            f"  Common causes:\n"
            f"    - VEP cache version does not match VEP executable version\n"
            f"    - Reference FASTA does not match the cache species/assembly\n"
            f"    - A plugin file is corrupt, missing, or in an unexpected format\n"
            f"    - Input VCF contains chromosome names not present in the cache\n"
            f"      (e.g. '1' vs 'chr1')\n"
            f"\n"
            f"  Re-run setup to validate all paths:\n"
            f"    python cli.py --setup"
        )

    if not Path(output_vcf).exists():
        raise RuntimeError(
            f"VEP exited successfully but output file was not produced.\n"
            f"  Expected: {output_vcf}\n"
            f"  Check that the output directory is writable and has sufficient "
            f"disk space."
        )

    logger.info("VEP annotation complete. Output: %s", output_vcf)
    return output_vcf


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def annotate(
    input_vcf: str,
    output_vcf: str,
    paths_config: dict,
    settings: dict,
) -> str:
    """
    Annotate a validated VCF using command-line VEP.

    Orchestrates command construction and execution.  The resulting annotated
    VCF is consumed by analysis/parser.py in the next pipeline step.

    Parameters
    ----------
    input_vcf:    Validated VCF path (from vcf_validator.validate_vcf()).
    output_vcf:   Output file path.  Will be created (or overwritten).
                  The parent directory is created if it does not exist.
    paths_config: Parsed config/paths.yaml (load via config_writer.load_paths()).
    settings:     Parsed config/settings.yaml.

    Returns
    -------
    Path to the annotated VCF output file.
    """
    logger.info("Annotation step: %s → %s", input_vcf, output_vcf)

    # --- Skip-if-complete guard ---
    skip = settings.get("vep", {}).get("skip_if_annotated_exists", True)
    if skip and _output_is_complete(output_vcf):
        logger.info(
            "Skipping VEP annotation — output already exists and contains a valid "
            "CSQ header: %s  (set vep.skip_if_annotated_exists: false to force re-run)",
            output_vcf,
        )
        return output_vcf

    # Ensure output directory exists before VEP tries to write
    out_dir = Path(output_vcf).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    command = build_vep_command(input_vcf, output_vcf, paths_config, settings)
    return run_vep(command, output_vcf, settings)
