"""
run_log.py — Structured audit log for pipeline runs.

Captures: run timestamp, input file, sample/case ID, steps completed,
output file paths, summary counts, warnings, and errors.

Usage:
    log = RunLog("results/run_SAMPLE.log.json", metadata, paths_config)
    log.add_step("parse_vcf", "completed", detail="1234 variants parsed")
    log.add_warning("SpliceAI scores not available")
    log.add_output_file("xlsx", "results/SAMPLE.xlsx")
    log.finalize()          # writes JSON file, returns path

    # or as a context manager (finalize() called on exit):
    with RunLog(...) as log:
        log.add_step(...)
"""

import datetime
import json
import logging
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from reporting import PIPELINE_VERSION

logger = logging.getLogger(__name__)


class RunLog:
    """
    Accumulates structured log entries for one pipeline run and writes them
    as a JSON file when finalize() is called.
    """

    def __init__(
        self,
        output_path: str,
        metadata: dict,
        paths_config: Optional[dict] = None,
    ) -> None:
        self._output_path   = str(output_path)
        self._metadata      = dict(metadata)
        self._paths_config  = paths_config or {}
        self._started_at    = _now()
        self._completed_at: Optional[str] = None
        self._steps:        List[Dict[str, Any]] = []
        self._warnings:     List[str] = []
        self._errors:       List[str] = []
        self._output_files: Dict[str, str] = {}
        self._summary:      Dict[str, Any] = {}
        self._run_id        = uuid.uuid4().hex[:8]

    # ------------------------------------------------------------------
    # Accumulation methods
    # ------------------------------------------------------------------

    def add_step(
        self,
        name: str,
        status: str = "completed",
        detail: Optional[str] = None,
    ) -> None:
        """Record a named pipeline step with optional detail string."""
        entry: Dict[str, Any] = {
            "step":      name,
            "status":    status,
            "timestamp": _now(),
        }
        if detail:
            entry["detail"] = detail
        self._steps.append(entry)
        logger.debug("Step [%s]: %s — %s", status, name, detail or "")

    def add_warning(self, message: str) -> None:
        """Record a non-fatal warning."""
        self._warnings.append(message)
        logger.warning("RunLog: %s", message)

    def add_error(self, message: str) -> None:
        """Record an error (pipeline may still complete)."""
        self._errors.append(message)
        logger.error("RunLog: %s", message)

    def add_output_file(self, key: str, path: str) -> None:
        """Register an output file path under a descriptive key."""
        self._output_files[key] = str(path)

    def set_summary(self, summary: Dict[str, Any]) -> None:
        """Store a summary dict (e.g. variant counts) for the run record."""
        self._summary = summary

    # ------------------------------------------------------------------
    # Finalisation
    # ------------------------------------------------------------------

    def finalize(self) -> str:
        """Write the accumulated log to JSON and return the output path."""
        self._completed_at = _now()
        out = Path(self._output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(self._build_dict(), fh, indent=2, default=str)
        logger.info("Run log written: %s", out)
        return str(out)

    def to_dict(self) -> dict:
        """Return the log as a plain dict (without writing to disk)."""
        return self._build_dict()

    def _build_dict(self) -> dict:
        return {
            "run_id":           self._run_id,
            "pipeline_version": PIPELINE_VERSION,
            "started_at":       self._started_at,
            "completed_at":     self._completed_at,
            "metadata":         self._metadata,
            "steps_completed":  self._steps,
            "output_files":     self._output_files,
            "summary":          self._summary,
            "warnings":         self._warnings,
            "errors":           self._errors,
        }

    # ------------------------------------------------------------------
    # Context-manager support
    # ------------------------------------------------------------------

    def __enter__(self) -> "RunLog":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        if exc_type is not None:
            self.add_error(f"{exc_type.__name__}: {exc_val}")
        self.finalize()
        return False   # do not suppress exceptions


# ---------------------------------------------------------------------------
# Module-level factory (preserves stub signature)
# ---------------------------------------------------------------------------

def init_run_log(
    output_path: str,
    metadata: dict,
    paths_config: Optional[dict] = None,
) -> RunLog:
    """Initialise and return a RunLog instance ready to accumulate entries."""
    return RunLog(output_path, metadata, paths_config)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _now() -> str:
    """Return the current UTC time as an ISO-8601 string."""
    return datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z"
