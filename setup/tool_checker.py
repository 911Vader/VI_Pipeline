"""
tool_checker.py — Locate and validate required executables.

Searches PATH, common conda/system locations, and any user-supplied hints.
Does not assume any default install location.
"""

import shutil
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Minimum acceptable versions (major, minor).
# Set to (0, 0) to accept any version.
# ---------------------------------------------------------------------------
MIN_VERSIONS: Dict[str, Tuple[int, int]] = {
    "docker":    (20, 0),
    "conda":     (4,  0),
    "python":    (3, 10),
    "nextflow":  (23,  0),
    "vep":       (111, 0),
    "bcftools":  (1, 17),
    "bgzip":     (1, 17),
    "tabix":     (1, 17),
}

# Candidate names to try for each logical tool
# (handles e.g. "vep" vs "ensembl-vep")
TOOL_CANDIDATES: Dict[str, List[str]] = {
    "docker":   ["docker"],
    "conda":    ["conda"],
    "python":   ["python3", "python"],
    "nextflow": ["nextflow"],
    "vep":      ["vep", "variant_effect_predictor.pl"],
    "bcftools": ["bcftools"],
    "bgzip":    ["bgzip"],
    "tabix":    ["tabix"],
}


def _find_executable(candidates: List[str]) -> Optional[str]:
    """Return the absolute path of the first candidate found on PATH."""
    for name in candidates:
        found = shutil.which(name)
        if found:
            return found
    return None


def _run_version_command(exe: str, tool: str) -> Optional[str]:
    """
    Run the tool's version flag and return raw stdout+stderr.
    Returns None if the command fails or times out.
    """
    version_flags: Dict[str, List[str]] = {
        "docker":   ["--version"],
        "conda":    ["--version"],
        "python":   ["--version"],
        "nextflow": ["-version"],
        "vep":      ["--help"],          # VEP prints version in --help header
        "bcftools": ["--version"],
        "bgzip":    ["--version"],
        "tabix":    ["--version"],
    }
    flags = version_flags.get(tool, ["--version"])
    try:
        result = subprocess.run(
            [exe] + flags,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return (result.stdout + result.stderr).strip()
    except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError) as exc:
        logger.debug("Version command failed for %s: %s", tool, exc)
        return None


def _parse_version(raw: str, tool: str) -> Optional[Tuple[int, int]]:
    """
    Extract (major, minor) from raw version output.
    Returns None if parsing fails.
    """
    import re

    patterns: Dict[str, str] = {
        "docker":   r"Docker version (\d+)\.(\d+)",
        "conda":    r"conda (\d+)\.(\d+)",
        "python":   r"Python (\d+)\.(\d+)",
        "nextflow": r"version (\d+)\.(\d+)",
        "vep":      r"ensembl-vep\s*:\s*(\d+)",     # single major only
        "bcftools": r"bcftools (\d+)\.(\d+)",
        "bgzip":    r"bgzip \(htslib\) (\d+)\.(\d+)",
        "tabix":    r"tabix \(htslib\) (\d+)\.(\d+)",
    }

    pattern = patterns.get(tool)
    if not pattern:
        return None

    match = re.search(pattern, raw, re.IGNORECASE | re.MULTILINE)
    if not match:
        return None

    try:
        if len(match.groups()) == 1:
            return (int(match.group(1)), 0)
        return (int(match.group(1)), int(match.group(2)))
    except (ValueError, IndexError):
        return None


def _version_ok(found: Tuple[int, int], minimum: Tuple[int, int]) -> bool:
    return found >= minimum


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class ToolCheckResult:
    """Result for a single tool check."""

    def __init__(self, tool: str):
        self.tool = tool
        self.path: Optional[str] = None
        self.version_raw: Optional[str] = None
        self.version_parsed: Optional[Tuple[int, int]] = None
        self.version_ok: bool = False
        self.status: str = "not_checked"   # found | missing | version_fail
        self.message: str = ""

    def as_dict(self) -> dict:
        return {
            "path": self.path,
            "version": (
                f"{self.version_parsed[0]}.{self.version_parsed[1]}"
                if self.version_parsed else None
            ),
            "status": self.status,
        }


def check_tool(tool: str, hint_path: Optional[str] = None) -> ToolCheckResult:
    """
    Locate and validate a single tool.

    Parameters
    ----------
    tool:       Logical tool name (key in TOOL_CANDIDATES).
    hint_path:  Optional explicit path supplied by user/config.
    """
    result = ToolCheckResult(tool)
    candidates = []

    if hint_path:
        candidates.append(hint_path)
    candidates.extend(TOOL_CANDIDATES.get(tool, [tool]))

    # Find executable
    exe = None
    if hint_path and Path(hint_path).is_file():
        exe = hint_path
    else:
        exe = _find_executable(candidates)

    if not exe:
        result.status = "missing"
        result.message = (
            f"'{tool}' was not found on PATH or at any expected location.\n"
            f"  Tried: {', '.join(candidates)}\n"
            f"  Fix: install {tool} and ensure it is on your PATH, or supply "
            f"the path in config/paths.yaml."
        )
        return result

    result.path = exe

    # Get version
    raw = _run_version_command(exe, tool)
    if raw is None:
        result.status = "missing"
        result.message = (
            f"'{tool}' was found at {exe} but could not be executed. "
            f"Check permissions and dependencies."
        )
        return result

    result.version_raw = raw
    result.version_parsed = _parse_version(raw, tool)
    minimum = MIN_VERSIONS.get(tool, (0, 0))

    if result.version_parsed is None:
        # Found and runnable, but version could not be parsed — warn, don't fail.
        result.status = "found"
        result.version_ok = True
        result.message = (
            f"'{tool}' found at {exe}. "
            f"Version could not be parsed from output — proceeding with caution."
        )
    elif not _version_ok(result.version_parsed, minimum):
        result.status = "version_fail"
        result.version_ok = False
        found_str = f"{result.version_parsed[0]}.{result.version_parsed[1]}"
        min_str = f"{minimum[0]}.{minimum[1]}"
        result.message = (
            f"'{tool}' at {exe} is version {found_str}, "
            f"but version >= {min_str} is required. Please upgrade."
        )
    else:
        result.status = "found"
        result.version_ok = True
        v = result.version_parsed
        result.message = f"'{tool}' OK — {exe} (v{v[0]}.{v[1]})"

    return result


def check_all_tools(hint_paths: Optional[Dict[str, str]] = None) -> Dict[str, ToolCheckResult]:
    """
    Check all required tools.

    Parameters
    ----------
    hint_paths: Mapping of tool_name → explicit_path (from existing paths.yaml).

    Returns
    -------
    Dict of tool_name → ToolCheckResult.
    """
    hints = hint_paths or {}
    results = {}
    for tool in TOOL_CANDIDATES:
        results[tool] = check_tool(tool, hint_path=hints.get(tool))
    return results
