"""
tool_concordance.py — Pairwise in-silico tool concordance analysis.

For each variant, each tool with a score is classified as:
  - "damaging"  — score above its pathogenic threshold
  - "benign"    — score below its benign threshold
  - "gray"      — score in the uncertain zone (or absent)

Two tools are "concordant" for a given variant when both have a non-gray
call AND they agree on direction (both damaging, or both benign).

Pairwise concordance rate = concordant variants / variants where both tools
                             had a non-gray call.

Tools assessed:
  CADD         — pathogenic threshold from settings (default ≥20)
  REVEL        — pathogenic ≥0.75, benign ≤0.15
  AlphaMissense— pathogenic ≥0.564, benign ≤0.340
  SpliceAI     — high ≥0.5 (damaging), low ≤0.1 (benign)
  ClinVar      — P/LP vs B/LB (review-status ≥1 star required)
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from analysis.parser import Variant

logger = logging.getLogger(__name__)

TOOLS = ["CADD", "REVEL", "AlphaMissense", "SpliceAI", "ClinVar"]

_DAMAGING = "damaging"
_BENIGN   = "benign"
_GRAY     = "gray"


# ---------------------------------------------------------------------------
# Per-tool call
# ---------------------------------------------------------------------------

def _call_cadd(v: Variant, t: dict) -> str:
    if v.cadd_phred is None:
        return _GRAY
    return _DAMAGING if v.cadd_phred >= t["cadd_pathogenic"] else _BENIGN


def _call_revel(v: Variant, t: dict) -> str:
    if v.revel_score is None:
        return _GRAY
    if v.revel_score >= t["revel_pathogenic"]:
        return _DAMAGING
    if v.revel_score <= t["revel_benign"]:
        return _BENIGN
    return _GRAY


def _call_alphamissense(v: Variant, t: dict) -> str:
    if v.alphamissense_score is None:
        return _GRAY
    if v.alphamissense_score >= t["am_pathogenic"]:
        return _DAMAGING
    if v.alphamissense_score <= t["am_benign"]:
        return _BENIGN
    return _GRAY


def _call_spliceai(v: Variant, t: dict) -> str:
    if v.spliceai_max_score is None:
        return _GRAY
    if v.spliceai_max_score >= t["spliceai_high"]:
        return _DAMAGING
    if v.spliceai_max_score <= t["spliceai_low"]:
        return _BENIGN
    return _GRAY


def _call_clinvar(v: Variant) -> str:
    """ClinVar call: P/LP (≥1-star) = damaging; B/LB = benign; else gray."""
    if not v.clinvar_significance:
        return _GRAY
    sig = v.clinvar_significance.lower().replace("_", " ")
    # Require at least 1-star review (criteria provided, single submitter)
    rs  = (v.clinvar_review_status or "").lower()
    if "no_assertion" in rs.replace(" ", "_") or (
        "criteria" not in rs and "expert" not in rs and "guideline" not in rs
    ):
        return _GRAY
    if "conflicting" in sig or "uncertain" in sig:
        return _GRAY
    if sig.startswith("pathogenic") or sig.startswith("likely pathogenic"):
        return _DAMAGING
    if "benign" in sig:
        return _BENIGN
    return _GRAY


def tool_calls(v: Variant, t: dict) -> Dict[str, str]:
    """Return a {tool_name: call} dict for one variant."""
    return {
        "CADD":          _call_cadd(v, t),
        "REVEL":         _call_revel(v, t),
        "AlphaMissense": _call_alphamissense(v, t),
        "SpliceAI":      _call_spliceai(v, t),
        "ClinVar":       _call_clinvar(v),
    }


# ---------------------------------------------------------------------------
# Pairwise concordance
# ---------------------------------------------------------------------------

@dataclass
class ToolConcordanceStats:
    tool_a: str
    tool_b: str
    n_both_called: int  = 0   # variants where both tools had non-gray call
    n_concordant: int   = 0   # of those, how many agreed on direction
    n_discordant: int   = 0

    @property
    def concordance_rate(self) -> Optional[float]:
        if self.n_both_called == 0:
            return None
        return self.n_concordant / self.n_both_called

    @property
    def pair_label(self) -> str:
        return f"{self.tool_a} vs {self.tool_b}"


@dataclass
class ToolConcordanceReport:
    """Aggregated tool concordance data across a variant list."""
    pairwise: List[ToolConcordanceStats] = field(default_factory=list)
    n_majority_discordant: int = 0   # variants where in-silico tools mostly disagree
    n_clinvar_conflicts: int   = 0   # ClinVar direction conflicts with pipeline

    def matrix(self) -> Dict[Tuple[str, str], Optional[float]]:
        """Return {(tool_a, tool_b): rate} for easy matrix rendering."""
        return {(s.tool_a, s.tool_b): s.concordance_rate for s in self.pairwise}


def compute_tool_concordance(
    variants: List[Variant],
    settings: dict,
) -> ToolConcordanceReport:
    """
    Compute pairwise tool concordance and aggregate flags across all variants.
    """
    t = _thresholds(settings)

    # Initialise pairwise stats for all pairs
    pairs: Dict[Tuple[str, str], ToolConcordanceStats] = {}
    for i, a in enumerate(TOOLS):
        for b in TOOLS[i + 1:]:
            pairs[(a, b)] = ToolConcordanceStats(tool_a=a, tool_b=b)

    n_majority_discordant = 0
    n_clinvar_conflicts   = 0

    for v in variants:
        calls = tool_calls(v, t)

        # Pairwise concordance
        for (a, b), stat in pairs.items():
            ca, cb = calls[a], calls[b]
            if ca == _GRAY or cb == _GRAY:
                continue
            stat.n_both_called += 1
            if ca == cb:
                stat.n_concordant += 1
            else:
                stat.n_discordant += 1

        # Majority discordance among in-silico tools (CADD, REVEL, AM, SpliceAI)
        insilico_calls = [calls[t_] for t_ in ("CADD", "REVEL", "AlphaMissense", "SpliceAI")
                          if calls[t_] != _GRAY]
        if len(insilico_calls) >= 2:
            n_dam = sum(1 for c in insilico_calls if c == _DAMAGING)
            n_ben = sum(1 for c in insilico_calls if c == _BENIGN)
            if n_dam > 0 and n_ben > 0 and abs(n_dam - n_ben) <= 1:
                n_majority_discordant += 1

        # ClinVar vs pipeline direction conflict
        clinvar_call = calls["ClinVar"]
        pipe_cls = v.acmg_classification or ""
        if clinvar_call != _GRAY and pipe_cls:
            pipe_path = pipe_cls in {"Pathogenic", "Likely Pathogenic"}
            pipe_ben  = pipe_cls in {"Likely Benign", "Benign"}
            if (clinvar_call == _DAMAGING and pipe_ben) or \
               (clinvar_call == _BENIGN   and pipe_path):
                n_clinvar_conflicts += 1

    report = ToolConcordanceReport(
        pairwise             = list(pairs.values()),
        n_majority_discordant= n_majority_discordant,
        n_clinvar_conflicts  = n_clinvar_conflicts,
    )

    for stat in report.pairwise:
        if stat.n_both_called > 0:
            logger.info(
                "Tool concordance %s: %.1f%% (%d/%d)",
                stat.pair_label,
                (stat.concordance_rate or 0) * 100,
                stat.n_concordant,
                stat.n_both_called,
            )

    return report


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _thresholds(settings: dict) -> dict:
    t = settings.get("acmg", {}).get("insilico_thresholds", {})
    return {
        "cadd_pathogenic":  float(t.get("cadd_phred_pathogenic",   20)),
        "revel_pathogenic": float(t.get("revel_pathogenic",        0.75)),
        "revel_benign":     float(t.get("revel_benign",            0.15)),
        "spliceai_high":    float(t.get("spliceai_high",           0.5)),
        "spliceai_low":     float(t.get("spliceai_low",            0.1)),
        "am_pathogenic":    float(t.get("alphamissense_pathogenic", 0.564)),
        "am_benign":        float(t.get("alphamissense_benign",     0.340)),
    }
