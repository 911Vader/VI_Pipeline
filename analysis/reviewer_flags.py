"""
reviewer_flags.py — Generate reviewer-facing flags for each Variant.

Reviewer flags are SEPARATE from ACMG classification.  They are advisory
annotations that help a human expert prioritize and interpret variants but
do NOT change the ACMG classification.

Flags are generated after ACMG classification so they can reference
acmg_evidence_codes and acmg_classification.

Output:
  variant.reviewer_flags — list of human-readable flag strings
"""

import logging
from typing import List, Optional, Set

from analysis.parser import Variant

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Thresholds (with defaults matching settings.yaml)
# ---------------------------------------------------------------------------

def _thresholds(settings: dict) -> dict:
    t   = settings.get("acmg", {}).get("insilico_thresholds", {})
    rev = settings.get("reviewer_flags", {})
    return {
        "cadd_high":         float(rev.get("cadd_high",         t.get("cadd_phred_pathogenic", 20))),
        "cadd_very_high":    float(rev.get("cadd_very_high",    30)),
        "revel_high":        float(rev.get("revel_high",        t.get("revel_pathogenic", 0.75))),
        "spliceai_high":     float(rev.get("spliceai_high",     t.get("spliceai_high", 0.5))),
        "spliceai_moderate": float(rev.get("spliceai_moderate", 0.2)),
        "max_af":            float(settings.get("filtering", {}).get("max_af", 0.01)),
        "ba1_af":            float(settings.get("acmg", {}).get("ba1_af_threshold", 0.05)),
        "low_depth":         int(rev.get("low_depth",   20)),
        "low_gq":            int(rev.get("low_gq",      20)),
        "low_vaf":           float(rev.get("low_vaf",   0.20)),
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def annotate_reviewer_flags(
    variants: List[Variant],
    settings: dict,
    gene_list: Optional[List[str]] = None,
) -> List[Variant]:
    """
    Populate variant.reviewer_flags for every variant in the list.

    Called after ACMG classification and before ranking so that flags are
    available for reporting.  Modifies variants in place.
    """
    gene_set: Set[str] = set()
    if gene_list:
        gene_set = {g.strip().upper() for g in gene_list if g.strip()}

    t = _thresholds(settings)

    for variant in variants:
        variant.reviewer_flags = _build_flags(variant, t, gene_set)

    logger.info("Reviewer flags annotated for %d variants", len(variants))
    return variants


# ---------------------------------------------------------------------------
# Flag generation
# ---------------------------------------------------------------------------

def _build_flags(variant: Variant, t: dict, gene_set: Set[str]) -> List[str]:
    flags: List[str] = []

    # --- Gene list ---
    if gene_set:
        if variant.in_gene_list:
            flags.append("In supplied gene list")
        else:
            flags.append("Outside supplied gene list")

    # --- ClinGen gene-disease validity (mandatory flag when confidence is low) ---
    if variant.clingen_validity:
        validity_norm = variant.clingen_validity.lower().strip()
        if validity_norm == "refuted":
            flags.append(
                f"ClinGen validity: REFUTED — gene-disease relationship not supported; "
                f"classification unreliable (manual expert review required)"
            )
        elif validity_norm == "disputed":
            flags.append(
                f"ClinGen validity: DISPUTED — conflicting evidence for gene-disease "
                f"relationship; interpret with caution (manual expert review required)"
            )
        elif validity_norm == "limited":
            flags.append(
                f"ClinGen validity: LIMITED — insufficient evidence for gene-disease "
                f"relationship; classification confidence reduced to Low"
            )

    # --- Population frequency ---
    af = variant.gnomad_af
    if af is None:
        flags.append("Absent from gnomAD — novel or ultra-rare variant")
    elif af < 0.0001:
        flags.append(f"Very rare in gnomAD (AF={af:.2e})")
    elif af < t["max_af"]:
        flags.append(f"Rare in gnomAD (AF={af:.4f})")
    elif af < t["ba1_af"]:
        flags.append(f"Intermediate AF (AF={af:.4f})")
    else:
        flags.append(f"Common variant in gnomAD (AF={af:.4f})")

    # --- CADD ---
    if variant.cadd_phred is not None:
        if variant.cadd_phred >= t["cadd_very_high"]:
            flags.append(f"Very high CADD (phred={variant.cadd_phred:.1f})")
        elif variant.cadd_phred >= t["cadd_high"]:
            flags.append(f"High CADD (phred={variant.cadd_phred:.1f})")

    # --- REVEL ---
    if variant.revel_score is not None and variant.revel_score >= t["revel_high"]:
        flags.append(f"High REVEL ({variant.revel_score:.3f})")

    # --- SpliceAI ---
    if variant.spliceai_max_score is not None:
        if variant.spliceai_max_score >= t["spliceai_high"]:
            flags.append(f"High SpliceAI ({variant.spliceai_max_score:.3f}) — likely splice-altering")
        elif variant.spliceai_max_score >= t["spliceai_moderate"]:
            flags.append(f"Moderate SpliceAI ({variant.spliceai_max_score:.3f})")

    # --- LOFTEE ---
    if variant.loftee_lof == "HC":
        flags.append("LOFTEE high-confidence LoF")
    elif variant.loftee_lof == "LC":
        lc_note = f" ({variant.loftee_filter})" if variant.loftee_filter else ""
        flags.append(f"LOFTEE low-confidence LoF{lc_note}")

    # LoF consequence without LOFTEE HC — explicit flag
    _LOF_CSQ = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
                "splice_donor_variant", "start_lost"}
    if variant.consequence:
        csq_parts = {c.strip() for c in variant.consequence.split("&")}
        if csq_parts & _LOF_CSQ and variant.loftee_lof != "HC":
            flags.append("LoF consequence — LOFTEE HC not confirmed (PVS1 not applied)")

    # --- ClinVar ---
    if variant.clinvar_significance:
        sig_l = variant.clinvar_significance.lower().replace("_", " ")
        if "conflicting" in sig_l:
            flags.append("ClinVar: conflicting interpretations")
        elif sig_l.startswith("pathogenic") or sig_l.startswith("likely pathogenic"):
            review = f" — {variant.clinvar_review_status}" if variant.clinvar_review_status else ""
            flags.append(f"ClinVar: {variant.clinvar_significance}{review}")
        elif "benign" in sig_l:
            flags.append(f"ClinVar: {variant.clinvar_significance}")

    # --- QC ---
    if variant.depth is not None and variant.depth < t["low_depth"]:
        flags.append(f"Low read depth (DP={variant.depth})")
    if variant.genotype_quality is not None and variant.genotype_quality < t["low_gq"]:
        flags.append(f"Low genotype quality (GQ={variant.genotype_quality})")
    if variant.vaf is not None and variant.vaf < t["low_vaf"]:
        flags.append(f"Low VAF ({variant.vaf:.2f}) — possible artefact or mosaic")
    if variant.filter_status and variant.filter_status not in ("PASS", "."):
        flags.append(f"FILTER={variant.filter_status}")

    # --- ACMG context flags ---
    if variant.acmg_classification == "VUS" and af is not None and af < 0.001:
        flags.append("Rare VUS — warrants expert review")
    if "PS3" in (variant.acmg_evidence_codes or []):
        flags.append("Strong ClinVar evidence (PS3)")
    if "PM2" in (variant.acmg_evidence_codes or []) and variant.gnomad_af is None:
        flags.append(
            "PM2: absent from gnomAD — population database coverage limitations; "
            "variant may be in under-represented population; manual review required"
        )
    if variant.classification_confidence == "Low":
        flags.append(f"Classification confidence: Low — interpret with caution")

    # --- Deprioritization ---
    if variant.deprioritized and variant.deprioritization_reason:
        flags.append(f"Deprioritized: {variant.deprioritization_reason.split('(')[0].strip()}")

    return flags
