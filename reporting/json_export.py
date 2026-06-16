"""
json_export.py — Export ranked candidates to machine-readable JSON.

Output structure:
  {
    "pipeline_version": "1.0.0",
    "generated_at":     "<ISO timestamp>",
    "disclaimer":       "<mandatory text>",
    "metadata":         { sample_id, case_id, ... },
    "summary":          { total, retained, excluded, acmg_counts, ... },
    "retained_variants": [ { rank, priority_score, gene, ... }, ... ],
    "excluded_variants": [ { ... }, ... ]
  }

The disclaimer is always embedded verbatim — it must not be altered.
"""

import datetime
import json
import logging
from collections import Counter
from pathlib import Path
from typing import List, Optional

from analysis.parser import Variant
from reporting import PIPELINE_VERSION

logger = logging.getLogger(__name__)

#: Mandatory disclaimer — must appear verbatim in all outputs.
MANDATORY_DISCLAIMER = (
    "Top genomic candidate(s) requiring expert review. "
    "Variant pathogenicity alone does not prove cause of death."
)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_json(
    variants: List[Variant],
    output_path: str,
    metadata: dict,
    settings: Optional[dict] = None,
) -> str:
    """
    Write all ranked variants to a structured JSON file.

    Both retained and excluded variants are included in separate lists so
    downstream tools have access to the full call set.

    Parameters
    ----------
    variants:    Ranked list from ranker.rank_variants() (sorted in place).
    output_path: Destination file path (.json).
    metadata:    Dict with sample_id, case_id, and other run metadata.
    settings:    Parsed settings.yaml (used for summary only; optional).

    Returns
    -------
    Path to the written JSON file (str).
    """
    settings = settings or {}
    retained = [v for v in variants if not v.excluded]
    excluded = [v for v in variants if v.excluded]

    payload = {
        "pipeline_version": PIPELINE_VERSION,
        "generated_at":     _now(),
        "disclaimer":       MANDATORY_DISCLAIMER,
        "metadata":         metadata,
        "summary":          _build_summary(variants),
        "retained_variants": [
            _variant_to_dict(v, rank=i + 1) for i, v in enumerate(retained)
        ],
        "excluded_variants": [
            _variant_to_dict(v, rank=None) for v in excluded
        ],
    }

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)

    logger.info(
        "JSON export written: %s (%d retained, %d excluded)",
        out, len(retained), len(excluded),
    )
    return str(out)


# ---------------------------------------------------------------------------
# Serialisation helpers
# ---------------------------------------------------------------------------

def _variant_to_dict(v: Variant, rank: Optional[int]) -> dict:
    """Serialise one Variant to a JSON-compatible flat dict."""
    return {
        # Ranking
        "rank":                  rank,
        # Identity
        "gene_symbol":           v.gene_symbol or None,
        "chromosome":            v.chrom,
        "position":              v.pos,
        "ref":                   v.ref,
        "alt":                   v.alt,
        "variant_id":            v.variant_id or None,
        # Transcript annotation
        "transcript_id":         v.transcript_id or None,
        "hgvs_c":                v.hgvs_c,
        "hgvs_p":                v.hgvs_p,
        "consequence":           v.consequence or None,
        "zygosity":              v.zygosity,
        # ACMG classification
        "acmg_classification":   v.acmg_classification,
        "acmg_evidence_codes":   v.acmg_evidence_codes,
        # ClinVar
        "clinvar_significance":  v.clinvar_significance,
        "clinvar_review_status": v.clinvar_review_status,
        "clinvar_disease":       v.clinvar_disease,
        # Population frequency
        "gnomad_af":             v.gnomad_af,
        # In-silico scores
        "cadd_phred":            v.cadd_phred,
        "revel_score":           v.revel_score,
        "spliceai_max_score":    v.spliceai_max_score,
        "alphamissense_score":   v.alphamissense_score,
        # Context
        "in_gene_list":          v.in_gene_list,
        # Filtering
        "excluded":              v.excluded,
        "exclusion_reason":      v.exclusion_reason,
    }


def _build_summary(variants: List[Variant]) -> dict:
    """Compute summary counts for the JSON payload."""
    retained = [v for v in variants if not v.excluded]
    excluded = [v for v in variants if v.excluded]
    class_counts = Counter(
        v.acmg_classification for v in retained
        if v.acmg_classification
    )
    return {
        "total_variants":            len(variants),
        "retained":                  len(retained),
        "excluded":                  len(excluded),
        "acmg_classification_counts": dict(class_counts),
        "top_acmg_class":            retained[0].acmg_classification if retained else None,
        "top_gene":                  retained[0].gene_symbol if retained else None,
    }


def _now() -> str:
    return datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z"
