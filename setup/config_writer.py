"""
config_writer.py — Write validated paths to config/paths.yaml.

Only called after all required checks have passed.
Writes exactly what was validated — no defaults, no invented paths.
"""

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

import yaml

from setup.tool_checker import ToolCheckResult
from setup.resource_checker import ResourceCheckResult, RESOURCE_SPEC

logger = logging.getLogger(__name__)

PATHS_YAML = Path(__file__).parent.parent / "config" / "paths.yaml"


def _tool_entry(result: ToolCheckResult) -> dict:
    v = result.version_parsed
    return {
        "path": result.path,
        "version": f"{v[0]}.{v[1]}" if v else None,
        "status": result.status,
    }


def _resource_entry(result: ResourceCheckResult) -> dict:
    """
    Serialise a ResourceCheckResult to a paths.yaml entry dict.

    For optional resources that were not provided, the notes field is
    populated with a brief hint about what file is needed, so users can
    read paths.yaml and know exactly what to download and where to put it.
    """
    notes = None
    if result.status == "not_provided":
        spec = RESOURCE_SPEC.get(result.resource, {})
        if not spec.get("required", False):
            raw_hint = spec.get("user_hint", "")
            if raw_hint:
                # Keep only the first line for the notes field — full details
                # are shown during setup output.
                notes = raw_hint.splitlines()[0]
    return {
        "path": result.path,
        "status": result.status,
        "notes": notes,
    }


def write_paths(
    tool_results: Dict[str, ToolCheckResult],
    resource_results: Dict[str, ResourceCheckResult],
    setup_version: str,
    success: bool,
) -> None:
    """
    Serialise validated tool and resource paths to config/paths.yaml.

    Parameters
    ----------
    tool_results:     Output of check_all_tools().
    resource_results: Output of check_all_resources().
    setup_version:    Pipeline version string.
    success:          True if all required checks passed.
    """
    timestamp = datetime.now(timezone.utc).isoformat()

    doc = {
        "setup_completed": success,
        "setup_timestamp": timestamp,
        "setup_version": setup_version,
        "tools": {name: _tool_entry(r) for name, r in tool_results.items()},
        "resources": {name: _resource_entry(r) for name, r in resource_results.items()},
    }

    header = (
        "# =============================================================================\n"
        "# VI_Pipeline - Validated Tool and Resource Paths\n"
        "#\n"
        "# tools:     Auto-detected by setup. Do not edit this section directly.\n"
        "#\n"
        "# resources: Required resources are auto-detected.\n"
        "#            Optional resources (status: not_provided) must be configured\n"
        "#            manually — the 'notes' field shows what file is needed.\n"
        "#\n"
        "# To add an optional resource:\n"
        "#   1. Set its path:  resources.<name>.path: /your/path/to/file\n"
        "#   2. Re-run setup:  python cli.py --setup --force\n"
        "#      Setup will validate the path and update the status field.\n"
        "#\n"
        f"# Generated: {timestamp}\n"
        "# =============================================================================\n\n"
    )

    PATHS_YAML.parent.mkdir(parents=True, exist_ok=True)
    with open(PATHS_YAML, "w") as fh:
        fh.write(header)
        yaml.dump(doc, fh, default_flow_style=False, sort_keys=False)

    logger.info("paths.yaml written to %s (setup_completed=%s)", PATHS_YAML, success)


def load_paths() -> Optional[dict]:
    """
    Load config/paths.yaml.

    Returns the parsed dict, or None if the file does not exist.
    Does NOT validate contents — call the checker modules for that.
    """
    if not PATHS_YAML.exists():
        return None
    with open(PATHS_YAML) as fh:
        return yaml.safe_load(fh)


def assert_setup_complete() -> dict:
    """
    Load paths.yaml and raise a clear RuntimeError if setup has not been completed.

    Called by analysis and reporting modules before doing any work.
    """
    data = load_paths()
    if data is None:
        raise RuntimeError(
            "config/paths.yaml not found.\n"
            "  Run setup first:   python cli.py --setup\n"
            "  This validates all required tools and resources before analysis."
        )
    if not data.get("setup_completed", False):
        raise RuntimeError(
            "Setup did not complete successfully.\n"
            "  Re-run:   python cli.py --setup\n"
            "  Review the errors reported and fix them before running analysis."
        )
    return data
