"""
evidence_codes.py — ACMG/AMP evidence code definitions and combination rules.

References:
  Richards et al. (2015) Genetics in Medicine 17:405-424
  ClinGen sequence variant interpretation working group recommendations

This file defines:
  - The full catalogue of evidence codes (v1 subset active, rest defined but inactive)
  - Evidence weight categories
  - Classification combination rules per Table 5 of Richards 2015
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List


class EvidenceDirection(Enum):
    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class EvidenceStrength(Enum):
    VERY_STRONG = "very_strong"    # PVS
    STRONG = "strong"              # PS / BS
    MODERATE = "moderate"          # PM
    SUPPORTING = "supporting"      # PP / BP
    STAND_ALONE = "stand_alone"    # BA


@dataclass(frozen=True)
class EvidenceCode:
    code: str
    direction: EvidenceDirection
    strength: EvidenceStrength
    description: str
    implemented_v1: bool           # True = evaluated in version 1


# ---------------------------------------------------------------------------
# Full catalogue — implemented_v1=True for the v1 active subset only
# ---------------------------------------------------------------------------
EVIDENCE_CATALOGUE: Dict[str, EvidenceCode] = {
    # --- Pathogenic Very Strong ---
    "PVS1": EvidenceCode(
        "PVS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.VERY_STRONG,
        "Null variant (nonsense, frameshift, canonical ±1/2 splice, initiation codon, "
        "exon deletion) in gene where LoF is disease mechanism",
        implemented_v1=True,
    ),

    # --- Pathogenic Strong ---
    "PS1": EvidenceCode(
        "PS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG,
        "Same amino acid change as established pathogenic variant (different nucleotide)",
        implemented_v1=False,   # Requires per-amino-acid ClinVar lookup — not available in v1
    ),
    "PS2": EvidenceCode(
        "PS2", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG,
        "De novo variant (confirmed paternity/maternity) in patient with disease",
        implemented_v1=False,   # Requires trio data — not available in v1
    ),
    "PS3": EvidenceCode(
        "PS3", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG,
        "Well-established functional studies show damaging effect (ClinVar P/LP proxy in v1)",
        implemented_v1=True,
    ),
    "PS4": EvidenceCode(
        "PS4", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG,
        "Prevalence of variant significantly increased in affected vs controls",
        implemented_v1=False,   # Requires case-control data
    ),

    # --- Pathogenic Moderate ---
    "PM1": EvidenceCode(
        "PM1", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Located in mutational hotspot or critical functional domain (no benign variation)",
        implemented_v1=False,   # Requires curated hotspot/functional-domain database — not available in v1
    ),
    "PM2": EvidenceCode(
        "PM2", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Absent or extremely rare in population databases (gnomAD AF < threshold)",
        implemented_v1=True,
    ),
    "PM3": EvidenceCode(
        "PM3", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Detected in trans with pathogenic variant for recessive disorder",
        implemented_v1=False,   # Requires phase information
    ),
    "PM4": EvidenceCode(
        "PM4", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Protein length change due to in-frame indel or stop-loss",
        implemented_v1=True,
    ),
    "PM5": EvidenceCode(
        "PM5", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Novel missense at codon where different missense is known pathogenic",
        implemented_v1=False,   # Requires per-codon ClinVar lookup — not available in v1
    ),
    "PM6": EvidenceCode(
        "PM6", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE,
        "Assumed de novo (paternity/maternity not confirmed)",
        implemented_v1=False,
    ),

    # --- Pathogenic Supporting ---
    "PP1": EvidenceCode(
        "PP1", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING,
        "Co-segregation with disease in affected family members",
        implemented_v1=False,   # Requires family data
    ),
    "PP2": EvidenceCode(
        "PP2", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING,
        "Missense variant in gene with low rate of benign missense variation",
        implemented_v1=False,   # Requires per-gene missense-constraint data — not available in v1
    ),
    "PP3": EvidenceCode(
        "PP3", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING,
        "Multiple in-silico tools predict damaging effect (CADD, REVEL, SpliceAI)",
        implemented_v1=True,
    ),
    "PP4": EvidenceCode(
        "PP4", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING,
        "Patient phenotype highly specific to disease with single genetic aetiology",
        implemented_v1=False,
    ),
    "PP5": EvidenceCode(
        "PP5", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING,
        "Reputable source reports as pathogenic (not used as primary criterion in v1)",
        implemented_v1=False,
    ),

    # --- Benign Stand-Alone ---
    "BA1": EvidenceCode(
        "BA1", EvidenceDirection.BENIGN, EvidenceStrength.STAND_ALONE,
        "AF > 5% in population databases (gnomAD) — Benign stand-alone",
        implemented_v1=True,
    ),

    # --- Benign Strong ---
    "BS1": EvidenceCode(
        "BS1", EvidenceDirection.BENIGN, EvidenceStrength.STRONG,
        "AF greater than expected for disorder",
        implemented_v1=True,
    ),
    "BS2": EvidenceCode(
        "BS2", EvidenceDirection.BENIGN, EvidenceStrength.STRONG,
        "Observed in healthy adults for fully penetrant condition",
        implemented_v1=False,   # Requires healthy-adult zygosity counts — not available in v1
    ),
    "BS3": EvidenceCode(
        "BS3", EvidenceDirection.BENIGN, EvidenceStrength.STRONG,
        "Well-established functional studies show no damaging effect",
        implemented_v1=False,
    ),
    "BS4": EvidenceCode(
        "BS4", EvidenceDirection.BENIGN, EvidenceStrength.STRONG,
        "Lack of segregation in affected family members",
        implemented_v1=False,
    ),

    # --- Benign Supporting ---
    "BP1": EvidenceCode(
        "BP1", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Missense in gene where only truncating variants cause disease",
        implemented_v1=False,
    ),
    "BP2": EvidenceCode(
        "BP2", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Observed in trans with pathogenic variant for dominant disorder",
        implemented_v1=False,
    ),
    "BP3": EvidenceCode(
        "BP3", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "In-frame indel in repetitive region without known function",
        implemented_v1=False,
    ),
    "BP4": EvidenceCode(
        "BP4", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Multiple in-silico tools predict benign effect",
        implemented_v1=True,
    ),
    "BP5": EvidenceCode(
        "BP5", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Variant found in case with alternate molecular basis for disease",
        implemented_v1=False,
    ),
    "BP6": EvidenceCode(
        "BP6", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Reputable source reports variant as benign",
        implemented_v1=False,
    ),
    "BP7": EvidenceCode(
        "BP7", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING,
        "Synonymous variant with no predicted splice impact",
        implemented_v1=True,
    ),
}


# ---------------------------------------------------------------------------
# Classification combination rules
# Based on Table 5, Richards et al. 2015
# Returns one of: Pathogenic | Likely Pathogenic | VUS | Likely Benign | Benign
# ---------------------------------------------------------------------------

ClassificationCategory = str   # Type alias for readability

PATHOGENIC       = "Pathogenic"
LIKELY_PATHOGENIC = "Likely Pathogenic"
VUS              = "VUS"
LIKELY_BENIGN    = "Likely Benign"
BENIGN           = "Benign"


def _count(
    codes: List[str],
    strength: EvidenceStrength,
    direction: EvidenceDirection,
) -> int:
    """Count triggered evidence codes matching a given strength and direction."""
    return sum(
        1 for c in codes
        if c in EVIDENCE_CATALOGUE
        and EVIDENCE_CATALOGUE[c].strength == strength
        and EVIDENCE_CATALOGUE[c].direction == direction
    )


def classify(evidence_codes: List[str]) -> ClassificationCategory:
    """
    Combine triggered evidence codes into an ACMG/AMP classification category.

    Implements Table 5 from Richards et al. 2015 Genetics in Medicine.
    Each combination rule is a direct transcription of the published table.

    Conflicting evidence handling (v1 conservative approach):
      - BA1 (stand-alone Benign) overrides all pathogenic evidence.
      - 2× Strong Benign overrides all pathogenic evidence.
      - Conflicting Likely Pathogenic + Likely Benign → VUS.
      - Pathogenic evidence takes precedence over Likely Benign.

    Parameters
    ----------
    evidence_codes: List of triggered code strings, e.g. ["PVS1", "PM2", "PP3"].
                    Unknown / not-in-catalogue codes are silently ignored.

    Returns
    -------
    One of: "Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"
    """
    if not evidence_codes:
        return VUS

    # --- Count by category ---
    pvs = _count(evidence_codes, EvidenceStrength.VERY_STRONG, EvidenceDirection.PATHOGENIC)
    ps  = _count(evidence_codes, EvidenceStrength.STRONG,      EvidenceDirection.PATHOGENIC)
    pm  = _count(evidence_codes, EvidenceStrength.MODERATE,    EvidenceDirection.PATHOGENIC)
    pp  = _count(evidence_codes, EvidenceStrength.SUPPORTING,  EvidenceDirection.PATHOGENIC)
    ba  = _count(evidence_codes, EvidenceStrength.STAND_ALONE, EvidenceDirection.BENIGN)
    bs  = _count(evidence_codes, EvidenceStrength.STRONG,      EvidenceDirection.BENIGN)
    bp  = _count(evidence_codes, EvidenceStrength.SUPPORTING,  EvidenceDirection.BENIGN)

    # --- Benign (stand-alone BA1 overrides all pathogenic evidence) ---
    if ba >= 1:
        return BENIGN
    if bs >= 2:
        return BENIGN

    # --- Pathogenic (Table 5 rules, checked exhaustively) ---
    is_pathogenic = (
        # PVS1 + strong/moderate/supporting combinations
        (pvs >= 1 and ps >= 1)                        or  # PVS1 + ≥1 PS
        (pvs >= 1 and pm >= 2)                        or  # PVS1 + ≥2 PM
        (pvs >= 1 and pm >= 1 and pp >= 1)            or  # PVS1 + 1 PM + ≥1 PP
        (pvs >= 1 and pp >= 2)                        or  # PVS1 + ≥2 PP
        # Strong combinations
        (ps >= 2)                                     or  # ≥2 PS
        (ps >= 1 and pm >= 3)                         or  # 1 PS + ≥3 PM
        (ps >= 1 and pm >= 2 and pp >= 2)             or  # 1 PS + ≥2 PM + ≥2 PP
        (ps >= 1 and pm >= 1 and pp >= 4)                 # 1 PS + 1 PM + ≥4 PP
    )

    # --- Likely Pathogenic (Table 5 rules) ---
    # Note: these are evaluated AFTER Pathogenic, so stricter cases are already
    # captured above.  The conditions below correctly identify the LP-only cases.
    is_likely_pathogenic = (
        (pvs >= 1 and pm >= 1)                        or  # PVS1 + 1 PM
        (ps >= 1 and pm >= 1)                         or  # 1 PS + 1-2 PM (≥3 PM → P above)
        (ps >= 1 and pp >= 2)                         or  # 1 PS + ≥2 PP
        (pm >= 3)                                     or  # ≥3 PM
        (pm >= 2 and pp >= 2)                         or  # 2 PM + ≥2 PP
        (pm >= 1 and pp >= 4)                             # 1 PM + ≥4 PP
    )

    # --- Likely Benign (Table 5 rules) ---
    is_likely_benign = (
        (bs >= 1 and bp >= 1)                         or  # 1 BS + 1 BP
        (bp >= 2)                                         # ≥2 BP
    )

    # --- Resolution (most severe classification wins unless conflicting at same level) ---
    if is_pathogenic:
        return PATHOGENIC
    if is_likely_pathogenic:
        # Conflicting LP + LB → VUS (conservative)
        if is_likely_benign:
            return VUS
        return LIKELY_PATHOGENIC
    if is_likely_benign:
        return LIKELY_BENIGN
    return VUS
