"""
setup_runner.py — Orchestrates the full setup check sequence.

Called by cli.py --setup. Prints beginner-friendly output to stdout.
Writes paths.yaml on success (or partial success with warnings).
"""

import logging
import sys
from typing import Dict, Optional

import yaml

from setup.tool_checker import check_all_tools, ToolCheckResult
from setup.resource_checker import check_all_resources, ResourceCheckResult, RESOURCE_SPEC
from setup.config_writer import write_paths, PATHS_YAML

logger = logging.getLogger(__name__)

PIPELINE_VERSION = "1.0.0"

# Tools that MUST be found for setup to succeed
REQUIRED_TOOLS = {"docker", "python", "vep", "bcftools", "bgzip", "tabix"}
# Tools that are expected but not fatal if absent in VCF-first mode
OPTIONAL_TOOLS = {"conda", "nextflow"}


def _print_section(title: str) -> None:
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def _print_ok(msg: str) -> None:
    print(f"  [OK]  {msg}")


def _print_warn(msg: str) -> None:
    print(f"  [WARN] {msg}")


def _print_fail(msg: str) -> None:
    print(f"  [FAIL] {msg}")


def _load_existing_paths() -> dict:
    """Load existing paths.yaml for user-supplied resource paths, if present."""
    if PATHS_YAML.exists():
        with open(PATHS_YAML) as fh:
            data = yaml.safe_load(fh) or {}
        return data
    return {}


def _extract_resource_paths(existing: dict) -> Dict[str, Optional[str]]:
    """Pull resource paths out of existing paths.yaml (user may have filled them in)."""
    resources_block = existing.get("resources", {})
    paths = {}
    for key, val in resources_block.items():
        if isinstance(val, dict):
            paths[key] = val.get("path")
        else:
            paths[key] = None
    return paths


def run_setup(force: bool = False) -> bool:
    """
    Run the full setup check sequence.

    Parameters
    ----------
    force: If True, re-check even if paths.yaml already shows setup_completed.

    Returns
    -------
    True if setup succeeded (all required checks passed).
    """
    _print_section("VI_Pipeline Setup Check")
    print(f"  Pipeline version : {PIPELINE_VERSION}")
    print(f"  Config target    : {PATHS_YAML}")

    existing = _load_existing_paths()
    if existing.get("setup_completed") and not force:
        print(
            "\n  Setup was already completed successfully.\n"
            "  Use --force to re-run all checks.\n"
        )
        return True

    # -----------------------------------------------------------------------
    # 1. Tool checks
    # -----------------------------------------------------------------------
    _print_section("1 / 2  Checking required tools")

    existing_tool_hints = {}
    tools_block = existing.get("tools", {})
    for name, val in tools_block.items():
        if isinstance(val, dict) and val.get("path"):
            existing_tool_hints[name] = val["path"]

    tool_results = check_all_tools(hint_paths=existing_tool_hints)
    tools_ok = True

    for name, result in tool_results.items():
        if result.status == "found":
            _print_ok(result.message)
        elif name in OPTIONAL_TOOLS:
            _print_warn(result.message)
        else:
            _print_fail(result.message)
            tools_ok = False

    # -----------------------------------------------------------------------
    # 2. Resource checks
    # -----------------------------------------------------------------------
    _print_section("2 / 2  Checking annotation resources")

    resource_paths = _extract_resource_paths(existing)
    resource_results = check_all_resources(resource_paths)
    resources_ok = True
    not_provided_optional = []

    for name, result in resource_results.items():
        spec = RESOURCE_SPEC.get(name, {})
        required = spec.get("required", False)
        for msg in result.messages:
            if result.status == "found":
                _print_ok(msg)
            elif result.status == "not_provided" and not required:
                _print_warn(msg)
            elif result.status == "not_provided" and required:
                _print_fail(msg)
                resources_ok = False
            elif result.status == "missing" and required:
                _print_fail(msg)
                resources_ok = False
            else:
                _print_warn(msg)
        if result.status == "not_provided" and not required:
            not_provided_optional.append(name)

    if not_provided_optional:
        print(
            f"\n  {len(not_provided_optional)} optional resource(s) not configured "
            f"(pipeline will run without them):\n"
            f"    {', '.join(not_provided_optional)}\n"
            f"\n"
            f"  To configure any of these:\n"
            f"    1. Open config/paths.yaml\n"
            f"    2. Set the path:  resources.<name>.path: /your/path/to/file\n"
            f"       (The 'notes' field in paths.yaml shows what file is needed)\n"
            f"    3. Re-run:        python cli.py --setup --force\n"
        )

    # -----------------------------------------------------------------------
    # 3. Write paths.yaml
    # -----------------------------------------------------------------------
    success = tools_ok and resources_ok
    write_paths(tool_results, resource_results, PIPELINE_VERSION, success)

    _print_section("Setup Summary")
    if success:
        print("  All required checks passed.")
        print(f"  Paths saved to: {PATHS_YAML}")
        print("\n  You can now run analysis:")
        print("    python cli.py --input sample.vcf.gz --gene-list config/gene_list.txt\n")
    else:
        print("  Setup FAILED. Resolve the errors above before running analysis.")
        print("\n  Common fixes:")
        print("    - Install missing tools and ensure they are on your PATH")
        print("    - Add resource file paths to config/paths.yaml under 'resources'")
        print("    - Re-run: python cli.py --setup\n")

    return success
