"""
Tests for analysis/acmg/evidence_codes.py

Validates the evidence code catalogue is internally consistent.
Run: pytest tests/analysis/acmg/test_evidence_codes.py -v
"""

import pytest
from analysis.acmg.evidence_codes import (
    EVIDENCE_CATALOGUE,
    EvidenceDirection,
    EvidenceStrength,
    PATHOGENIC, LIKELY_PATHOGENIC, VUS, LIKELY_BENIGN, BENIGN,
)


class TestEvidenceCatalogue:
    def test_all_expected_codes_present(self):
        expected = {
            "PVS1",
            "PS1", "PS2", "PS3", "PS4",
            "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
            "PP1", "PP2", "PP3", "PP4", "PP5",
            "BA1",
            "BS1", "BS2", "BS3", "BS4",
            "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
        }
        assert expected.issubset(set(EVIDENCE_CATALOGUE.keys()))

    def test_pvs1_is_pathogenic_very_strong(self):
        code = EVIDENCE_CATALOGUE["PVS1"]
        assert code.direction == EvidenceDirection.PATHOGENIC
        assert code.strength == EvidenceStrength.VERY_STRONG
        assert code.implemented_v1 is True

    def test_ba1_is_benign_stand_alone(self):
        code = EVIDENCE_CATALOGUE["BA1"]
        assert code.direction == EvidenceDirection.BENIGN
        assert code.strength == EvidenceStrength.STAND_ALONE
        assert code.implemented_v1 is True

    def test_ps2_not_implemented_in_v1(self):
        # PS2 requires trio data — must not be active in v1
        assert EVIDENCE_CATALOGUE["PS2"].implemented_v1 is False

    def test_all_pathogenic_codes_have_pathogenic_direction(self):
        pathogenic_prefixes = ("PVS", "PS", "PM", "PP")
        for code, entry in EVIDENCE_CATALOGUE.items():
            if any(code.startswith(p) for p in pathogenic_prefixes):
                assert entry.direction == EvidenceDirection.PATHOGENIC, (
                    f"{code} should be PATHOGENIC direction"
                )

    def test_all_benign_codes_have_benign_direction(self):
        benign_prefixes = ("BA", "BS", "BP")
        for code, entry in EVIDENCE_CATALOGUE.items():
            if any(code.startswith(p) for p in benign_prefixes):
                assert entry.direction == EvidenceDirection.BENIGN, (
                    f"{code} should be BENIGN direction"
                )

    def test_v1_active_subset_is_not_empty(self):
        active = [k for k, v in EVIDENCE_CATALOGUE.items() if v.implemented_v1]
        assert len(active) >= 9, "v1 should activate at least 9 criteria"

    def test_each_code_has_description(self):
        for code, entry in EVIDENCE_CATALOGUE.items():
            assert entry.description, f"{code} is missing a description"
