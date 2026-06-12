"""
Tests for analysis/parser.py

All tests use synthetic VCF content written to pytest tmp_path fixtures.
No real VEP output or genomic data is required.

Run: pytest tests/analysis/test_parser.py -v
"""

import gzip
import pytest
from pathlib import Path
from typing import Dict, Optional

from analysis.parser import (
    # Public
    parse_vep_output,
    Variant,
    TranscriptAnnotation,
    CONSEQUENCE_SEVERITY,
    # Private helpers (tested directly for unit coverage)
    _extract_csq_fields,
    _has_sample_column,
    _parse_info_str,
    _parse_csq_value,
    _parse_single_csq_entry,
    _select_primary_transcript,
    _check_multiallelic,
    _extract_zygosity,
    _parse_score_pred,
    _safe_float,
    _safe_str,
    _consequence_severity_rank,
    _extract_spliceai_max,
    _extract_revel,
    _extract_clinvar,
)

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

# A representative CSQ format string (minimal for most tests)
MIN_CSQ_FORMAT = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE"
    "|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position"
    "|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS"
    "|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen"
)

# Extended CSQ format including all optional plugin fields
FULL_CSQ_FORMAT = (
    MIN_CSQ_FORMAT
    + "|AF|gnomAD_AF|MAX_AF"
    + "|CADD_PHRED|REVEL"
    + "|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL"
    + "|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL"
    + "|am_pathogenicity"
)


def _csq_entry(
    fmt: str,
    allele: str = "T",
    consequence: str = "missense_variant",
    impact: str = "MODERATE",
    symbol: str = "BRCA1",
    gene: str = "ENSG00000012048",
    feature: str = "ENST00000357654",
    biotype: str = "protein_coding",
    hgvs_c: str = "c.5123A>T",
    hgvs_p: str = "p.Asp1708Val",
    canonical: str = "YES",
    sift: str = "0.02(deleterious)",
    polyphen: str = "0.95(probably_damaging)",
    gnomad_af: str = "0.001",
    max_af: str = "0.002",
    cadd_phred: str = "28.5",
    revel: str = "0.82",
    spliceai_ag: str = "0.05",
    spliceai_al: str = "0.02",
    spliceai_dg: str = "0.01",
    spliceai_dl: str = "0.01",
    am_path: str = "0.95",
    **overrides,
) -> str:
    """Build a CSQ pipe-delimited entry string for the given format."""
    defaults: Dict[str, str] = {
        "Allele": allele,
        "Consequence": consequence,
        "IMPACT": impact,
        "SYMBOL": symbol,
        "Gene": gene,
        "Feature_type": "Transcript",
        "Feature": feature,
        "BIOTYPE": biotype,
        "EXON": "10/23",
        "INTRON": "",
        "HGVSc": hgvs_c,
        "HGVSp": hgvs_p,
        "cDNA_position": "5123",
        "CDS_position": "5123",
        "Protein_position": "1708",
        "Amino_acids": "D/V",
        "Codons": "gAt/gTt",
        "Existing_variation": "",
        "DISTANCE": "",
        "STRAND": "1",
        "FLAGS": "",
        "VARIANT_CLASS": "SNV",
        "SYMBOL_SOURCE": "HGNC",
        "HGNC_ID": "1100",
        "CANONICAL": canonical,
        "SIFT": sift,
        "PolyPhen": polyphen,
        "AF": "",
        "gnomAD_AF": gnomad_af,
        "MAX_AF": max_af,
        "CADD_PHRED": cadd_phred,
        "REVEL": revel,
        "SpliceAI_pred_DS_AG": spliceai_ag,
        "SpliceAI_pred_DS_AL": spliceai_al,
        "SpliceAI_pred_DS_DG": spliceai_dg,
        "SpliceAI_pred_DS_DL": spliceai_dl,
        "am_pathogenicity": am_path,
        **overrides,
    }
    fields = fmt.split("|")
    return "|".join(defaults.get(f, "") for f in fields)


def _vcf_lines(
    records,
    csq_format: str = FULL_CSQ_FORMAT,
    has_sample: bool = True,
) -> str:
    """Build a complete VCF string from a list of data record strings."""
    header = [
        "##fileformat=VCFv4.2",
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description='
        f'"Consequence annotations from Ensembl VEP. Format: {csq_format}">',
    ]
    if has_sample:
        header.append(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        )
    else:
        header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(header + list(records)) + "\n"


def _record(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    filter_val: str = "PASS",
    csq: str = "",
    gt: str = "0/1",
    dp: str = "30",
    has_sample: bool = True,
) -> str:
    info = f"CSQ={csq}" if csq else "."
    if has_sample:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filter_val}\t{info}\tGT:DP\t{gt}:{dp}"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filter_val}\t{info}"


def _write_vcf(tmp_path: Path, content: str, filename: str = "test.vcf") -> Path:
    p = tmp_path / filename
    p.write_text(content, encoding="utf-8")
    return p


def _write_gz_vcf(tmp_path: Path, content: str, filename: str = "test.vcf.gz") -> Path:
    p = tmp_path / filename
    with gzip.open(str(p), "wt", encoding="utf-8") as fh:
        fh.write(content)
    return p


# ---------------------------------------------------------------------------
# _safe_float
# ---------------------------------------------------------------------------

class TestSafeFloat:
    def test_valid_float(self):
        assert _safe_float("0.95") == pytest.approx(0.95)

    def test_integer_string(self):
        assert _safe_float("25") == pytest.approx(25.0)

    def test_scientific_notation(self):
        assert _safe_float("1e-5") == pytest.approx(1e-5)

    def test_dot_returns_none(self):
        assert _safe_float(".") is None

    def test_empty_string_returns_none(self):
        assert _safe_float("") is None

    def test_na_returns_none(self):
        assert _safe_float("NA") is None

    def test_none_returns_none(self):
        assert _safe_float(None) is None

    def test_non_numeric_returns_none(self):
        assert _safe_float("deleterious") is None


# ---------------------------------------------------------------------------
# _safe_str
# ---------------------------------------------------------------------------

class TestSafeStr:
    def test_valid_string_returned(self):
        assert _safe_str("BRCA1") == "BRCA1"

    def test_empty_string_returns_none(self):
        assert _safe_str("") is None

    def test_dot_returns_none(self):
        assert _safe_str(".") is None

    def test_dash_returns_none(self):
        assert _safe_str("-") is None

    def test_none_input_returns_none(self):
        assert _safe_str(None) is None

    def test_whitespace_stripped(self):
        assert _safe_str("  BRCA1  ") == "BRCA1"


# ---------------------------------------------------------------------------
# _parse_score_pred (SIFT / PolyPhen)
# ---------------------------------------------------------------------------

class TestParseScorePred:
    def test_score_and_pred_parsed(self):
        score, pred = _parse_score_pred("0.02(deleterious)")
        assert score == pytest.approx(0.02)
        assert pred == "deleterious"

    def test_polyphen_probably_damaging(self):
        score, pred = _parse_score_pred("0.95(probably_damaging)")
        assert score == pytest.approx(0.95)
        assert pred == "probably_damaging"

    def test_empty_returns_none_none(self):
        assert _parse_score_pred("") == (None, None)

    def test_dot_returns_none_none(self):
        assert _parse_score_pred(".") == (None, None)

    def test_score_only_no_pred(self):
        score, pred = _parse_score_pred("0.75")
        assert score == pytest.approx(0.75)
        assert pred is None

    def test_tolerated(self):
        score, pred = _parse_score_pred("0.8(tolerated)")
        assert pred == "tolerated"


# ---------------------------------------------------------------------------
# _consequence_severity_rank
# ---------------------------------------------------------------------------

class TestConsequenceSeverityRank:
    def test_stop_gained_more_severe_than_missense(self):
        assert (
            _consequence_severity_rank(["stop_gained"])
            < _consequence_severity_rank(["missense_variant"])
        )

    def test_frameshift_more_severe_than_synonymous(self):
        assert (
            _consequence_severity_rank(["frameshift_variant"])
            < _consequence_severity_rank(["synonymous_variant"])
        )

    def test_compound_consequence_uses_most_severe(self):
        rank = _consequence_severity_rank(
            ["synonymous_variant", "splice_region_variant"]
        )
        assert rank == _consequence_severity_rank(["splice_region_variant"])

    def test_unknown_term_gets_default_rank(self):
        rank = _consequence_severity_rank(["made_up_variant_term"])
        assert rank == 9999

    def test_empty_list_gets_default_rank(self):
        assert _consequence_severity_rank([]) == 9999

    def test_splice_acceptor_more_severe_than_splice_region(self):
        assert (
            _consequence_severity_rank(["splice_acceptor_variant"])
            < _consequence_severity_rank(["splice_region_variant"])
        )


# ---------------------------------------------------------------------------
# _extract_csq_fields
# ---------------------------------------------------------------------------

class TestExtractCsqFields:
    def test_extracts_field_names_in_order(self):
        lines = [
            '##fileformat=VCFv4.2',
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description='
            f'"... Format: Allele|Consequence|IMPACT|SYMBOL">',
        ]
        fields = _extract_csq_fields(lines)
        assert fields == ["Allele", "Consequence", "IMPACT", "SYMBOL"]

    def test_whitespace_stripped_from_fields(self):
        lines = [
            '##INFO=<ID=CSQ,Number=.,Type=String,Description='
            '"x. Format: Allele | Consequence | IMPACT">',
        ]
        fields = _extract_csq_fields(lines)
        assert fields == ["Allele", "Consequence", "IMPACT"]

    def test_raises_when_no_csq_header(self):
        lines = ["##fileformat=VCFv4.2", "##FILTER=<ID=PASS>"]
        with pytest.raises(ValueError, match="CSQ annotation not found"):
            _extract_csq_fields(lines)

    def test_raises_with_helpful_message(self):
        with pytest.raises(ValueError, match="VEP"):
            _extract_csq_fields(["##fileformat=VCFv4.2"])

    def test_case_insensitive_match(self):
        lines = [
            '##info=<id=CSQ,number=.,type=String,description='
            '"x. format: Allele|Consequence">',
        ]
        fields = _extract_csq_fields(lines)
        assert "Allele" in fields


# ---------------------------------------------------------------------------
# _has_sample_column
# ---------------------------------------------------------------------------

class TestHasSampleColumn:
    def test_nine_cols_no_sample(self):
        line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        assert _has_sample_column(line) is False

    def test_ten_cols_has_sample(self):
        line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        assert _has_sample_column(line) is True


# ---------------------------------------------------------------------------
# _parse_info_str
# ---------------------------------------------------------------------------

class TestParseInfoStr:
    def test_key_value_pair(self):
        d = _parse_info_str("AF=0.01;DP=30")
        assert d["AF"] == "0.01"
        assert d["DP"] == "30"

    def test_flag_field(self):
        d = _parse_info_str("SOMATIC;AF=0.5")
        assert "SOMATIC" in d
        assert d["AF"] == "0.5"

    def test_csq_with_internal_commas(self):
        # CSQ value has commas inside — must not split prematurely
        csq = "T|missense_variant|MODERATE,T|synonymous_variant|LOW"
        d = _parse_info_str(f"CSQ={csq}")
        assert d["CSQ"] == csq

    def test_dot_returns_empty_dict(self):
        assert _parse_info_str(".") == {}

    def test_empty_string_returns_empty_dict(self):
        assert _parse_info_str("") == {}

    def test_value_with_equals_sign(self):
        # e.g. a base64-encoded field
        d = _parse_info_str("FOO=a=b=c")
        assert d["FOO"] == "a=b=c"


# ---------------------------------------------------------------------------
# _check_multiallelic
# ---------------------------------------------------------------------------

class TestCheckMultiallelic:
    def test_biallelic_does_not_raise(self):
        _check_multiallelic("chr1", 100, "A", "T")  # no exception

    def test_multiallelic_raises_value_error(self):
        with pytest.raises(ValueError, match="Multi-allelic"):
            _check_multiallelic("chr1", 100, "A", "T,G")

    def test_error_message_contains_bcftools_command(self):
        with pytest.raises(ValueError, match="bcftools norm"):
            _check_multiallelic("chr7", 500, "C", "A,G")

    def test_error_message_contains_position(self):
        with pytest.raises(ValueError, match="chr7:500"):
            _check_multiallelic("chr7", 500, "C", "A,G")

    def test_error_message_contains_fix_instructions(self):
        with pytest.raises(ValueError, match="split"):
            _check_multiallelic("chr1", 100, "A", "T,C")


# ---------------------------------------------------------------------------
# _extract_zygosity
# ---------------------------------------------------------------------------

class TestExtractZygosity:
    def test_heterozygous(self):
        assert _extract_zygosity("GT:DP", "0/1:30") == "heterozygous"

    def test_homozygous_alt(self):
        assert _extract_zygosity("GT:DP", "1/1:30") == "homozygous_alt"

    def test_homozygous_ref(self):
        assert _extract_zygosity("GT:DP", "0/0:30") == "homozygous_ref"

    def test_phased_heterozygous(self):
        assert _extract_zygosity("GT", "0|1") == "heterozygous"

    def test_missing_gt(self):
        assert _extract_zygosity("GT", "./.") == "unknown"

    def test_no_gt_in_format(self):
        assert _extract_zygosity("DP:GQ", "30:99") is None

    def test_empty_format_returns_none(self):
        assert _extract_zygosity("", "0/1:30") is None

    def test_gt_not_first_field(self):
        # GT can be anywhere in FORMAT
        assert _extract_zygosity("DP:GT", "30:0/1") == "heterozygous"


# ---------------------------------------------------------------------------
# _select_primary_transcript
# ---------------------------------------------------------------------------

class TestSelectPrimaryTranscript:
    def _ta(self, canonical=False, consequences=None, transcript_id="TX1"):
        c = consequences or ["missense_variant"]
        return TranscriptAnnotation(
            transcript_id=transcript_id,
            is_canonical=canonical,
            consequences=c,
            consequence="&".join(c),
            severity_rank=_consequence_severity_rank(c),
        )

    def test_returns_none_for_empty_list(self):
        assert _select_primary_transcript([]) is None

    def test_canonical_preferred_over_noncanonical(self):
        noncanon = self._ta(canonical=False, consequences=["stop_gained"], transcript_id="TX_nc")
        canon = self._ta(canonical=True, consequences=["missense_variant"], transcript_id="TX_c")
        primary = _select_primary_transcript([noncanon, canon])
        assert primary.transcript_id == "TX_c"

    def test_most_severe_canonical_selected_when_multiple_canonical(self):
        c1 = self._ta(canonical=True, consequences=["synonymous_variant"], transcript_id="TX1")
        c2 = self._ta(canonical=True, consequences=["stop_gained"], transcript_id="TX2")
        primary = _select_primary_transcript([c1, c2])
        assert primary.transcript_id == "TX2"

    def test_most_severe_noncanon_when_no_canonical(self):
        t1 = self._ta(canonical=False, consequences=["synonymous_variant"], transcript_id="TX1")
        t2 = self._ta(canonical=False, consequences=["frameshift_variant"], transcript_id="TX2")
        primary = _select_primary_transcript([t1, t2])
        assert primary.transcript_id == "TX2"

    def test_first_entry_used_when_all_severity_tied(self):
        t1 = self._ta(canonical=False, consequences=["intron_variant"], transcript_id="FIRST")
        t2 = self._ta(canonical=False, consequences=["intron_variant"], transcript_id="SECOND")
        primary = _select_primary_transcript([t1, t2])
        assert primary.transcript_id == "FIRST"

    def test_single_transcript_returned(self):
        t = self._ta(canonical=True, consequences=["missense_variant"])
        assert _select_primary_transcript([t]) is t

    def test_canonical_trumps_high_severity_noncanonical(self):
        """
        A canonical transcript with a low-severity consequence should still be
        preferred over a non-canonical with a high-severity consequence.
        """
        high_sev = self._ta(
            canonical=False, consequences=["transcript_ablation"], transcript_id="HIGH"
        )
        low_sev_canon = self._ta(
            canonical=True, consequences=["synonymous_variant"], transcript_id="CANON"
        )
        primary = _select_primary_transcript([high_sev, low_sev_canon])
        assert primary.transcript_id == "CANON"


# ---------------------------------------------------------------------------
# _extract_spliceai_max
# ---------------------------------------------------------------------------

class TestExtractSpliceaiMax:
    def test_max_of_four_scores(self):
        raw = {
            "SpliceAI_pred_DS_AG": "0.1",
            "SpliceAI_pred_DS_AL": "0.5",
            "SpliceAI_pred_DS_DG": "0.2",
            "SpliceAI_pred_DS_DL": "0.3",
        }
        assert _extract_spliceai_max(raw) == pytest.approx(0.5)

    def test_missing_all_fields_returns_none(self):
        assert _extract_spliceai_max({}) is None

    def test_partial_fields_uses_available(self):
        raw = {"SpliceAI_pred_DS_AG": "0.7"}
        assert _extract_spliceai_max(raw) == pytest.approx(0.7)

    def test_dot_values_return_none(self):
        raw = {k: "." for k in (
            "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
        )}
        assert _extract_spliceai_max(raw) is None


# ---------------------------------------------------------------------------
# _extract_revel
# ---------------------------------------------------------------------------

class TestExtractRevel:
    def test_single_revel_score(self):
        raw = {"REVEL": "0.82"}
        assert _extract_revel(raw) == pytest.approx(0.82)

    def test_revel_score_field_name(self):
        raw = {"REVEL_score": "0.75"}
        assert _extract_revel(raw) == pytest.approx(0.75)

    def test_semicolon_separated_values_takes_max(self):
        raw = {"REVEL_score": "0.6;0.85;0.7"}
        assert _extract_revel(raw) == pytest.approx(0.85)

    def test_missing_returns_none(self):
        assert _extract_revel({}) is None

    def test_dot_value_returns_none(self):
        assert _extract_revel({"REVEL": "."}) is None


# ---------------------------------------------------------------------------
# _extract_clinvar
# ---------------------------------------------------------------------------

class TestExtractClinvar:
    def test_all_fields_present(self):
        info = {
            "ClinVar_CLNSIG": "Pathogenic",
            "ClinVar_CLNREVSTAT": "criteria_provided,_multiple_submitters",
            "ClinVar_CLNDN": "Hereditary_breast_cancer",
        }
        sig, revstat, disease = _extract_clinvar(info)
        assert sig == "Pathogenic"
        assert revstat == "criteria_provided,_multiple_submitters"
        assert disease == "Hereditary_breast_cancer"

    def test_all_absent_returns_none_triple(self):
        sig, revstat, disease = _extract_clinvar({})
        assert sig is None
        assert revstat is None
        assert disease is None

    def test_partial_fields(self):
        sig, revstat, disease = _extract_clinvar({"ClinVar_CLNSIG": "Benign"})
        assert sig == "Benign"
        assert revstat is None
        assert disease is None


# ---------------------------------------------------------------------------
# _parse_csq_value
# ---------------------------------------------------------------------------

class TestParseCsqValue:
    def test_single_transcript_parsed(self):
        fmt = MIN_CSQ_FORMAT
        entry = _csq_entry(fmt, canonical="YES")
        transcripts = _parse_csq_value(entry, fmt.split("|"))
        assert len(transcripts) == 1
        assert transcripts[0].is_canonical is True
        assert transcripts[0].gene_symbol == "BRCA1"

    def test_multiple_transcripts_parsed(self):
        fmt = MIN_CSQ_FORMAT
        e1 = _csq_entry(fmt, feature="ENST00000001", canonical="YES")
        e2 = _csq_entry(fmt, feature="ENST00000002", canonical="")
        csq_val = f"{e1},{e2}"
        transcripts = _parse_csq_value(csq_val, fmt.split("|"))
        assert len(transcripts) == 2

    def test_empty_csq_value_returns_empty_list(self):
        assert _parse_csq_value("", MIN_CSQ_FORMAT.split("|")) == []

    def test_wrong_field_count_skips_entry(self):
        fields = MIN_CSQ_FORMAT.split("|")
        bad_entry = "only|two|fields"
        transcripts = _parse_csq_value(bad_entry, fields)
        assert len(transcripts) == 0

    def test_compound_consequence_split(self):
        fmt = MIN_CSQ_FORMAT
        entry = _csq_entry(fmt, consequence="missense_variant&splice_region_variant")
        transcripts = _parse_csq_value(entry, fmt.split("|"))
        assert "missense_variant" in transcripts[0].consequences
        assert "splice_region_variant" in transcripts[0].consequences


# ---------------------------------------------------------------------------
# _parse_single_csq_entry — optional fields
# ---------------------------------------------------------------------------

class TestParseSingleCsqEntryOptionalFields:
    def _make_entry(self, overrides=None):
        fmt = FULL_CSQ_FORMAT
        kw = overrides or {}
        entry = _csq_entry(fmt, **kw)
        fields = fmt.split("|")
        values = entry.split("|")
        return _parse_single_csq_entry(fields, values)

    def test_all_scores_populated(self):
        ta = self._make_entry()
        assert ta.cadd_phred == pytest.approx(28.5)
        assert ta.revel_score == pytest.approx(0.82)
        assert ta.spliceai_max_ds is not None
        assert ta.alphamissense_score == pytest.approx(0.95)

    def test_cadd_missing_returns_none(self):
        ta = self._make_entry({"cadd_phred": ""})
        assert ta.cadd_phred is None

    def test_revel_missing_returns_none(self):
        ta = self._make_entry({"revel": ""})
        assert ta.revel_score is None

    def test_spliceai_all_missing_returns_none(self):
        ta = self._make_entry({
            "spliceai_ag": "", "spliceai_al": "",
            "spliceai_dg": "", "spliceai_dl": "",
        })
        assert ta.spliceai_max_ds is None

    def test_alphamissense_missing_returns_none(self):
        ta = self._make_entry({"am_path": ""})
        assert ta.alphamissense_score is None

    def test_gnomad_af_populated(self):
        ta = self._make_entry()
        assert ta.af_gnomad == pytest.approx(0.001)
        assert ta.max_af == pytest.approx(0.002)

    def test_gnomad_af_missing_returns_none(self):
        ta = self._make_entry({"gnomad_af": "", "max_af": ""})
        assert ta.af_gnomad is None
        assert ta.max_af is None

    def test_sift_parsed(self):
        ta = self._make_entry()
        assert ta.sift_score == pytest.approx(0.02)
        assert ta.sift_pred == "deleterious"

    def test_polyphen_parsed(self):
        ta = self._make_entry()
        assert ta.polyphen_score == pytest.approx(0.95)
        assert ta.polyphen_pred == "probably_damaging"

    def test_sift_missing_returns_none(self):
        ta = self._make_entry({"sift": ""})
        assert ta.sift_score is None
        assert ta.sift_pred is None

    def test_is_canonical_true(self):
        ta = self._make_entry({"canonical": "YES"})
        assert ta.is_canonical is True

    def test_is_canonical_false_when_empty(self):
        ta = self._make_entry({"canonical": ""})
        assert ta.is_canonical is False

    def test_severity_rank_set(self):
        ta = self._make_entry({"consequence": "stop_gained"})
        assert ta.severity_rank < 9999
        assert ta.severity_rank == CONSEQUENCE_SEVERITY["stop_gained"]

    def test_raw_dict_contains_all_fields(self):
        ta = self._make_entry()
        assert "CADD_PHRED" in ta.raw
        assert "REVEL" in ta.raw


# ---------------------------------------------------------------------------
# parse_vep_output — integration (plain VCF)
# ---------------------------------------------------------------------------

class TestParseVepOutputPlainVcf:
    def test_single_variant_parsed(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        variants = parse_vep_output(str(vcf))
        assert len(variants) == 1

    def test_variant_position_correct(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(chrom="chr7", pos=1234, csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.chrom == "chr7"
        assert v.pos == 1234

    def test_ref_alt_correct(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(ref="G", alt="C", csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.ref == "G"
        assert v.alt == "C"

    def test_gene_symbol_from_primary_transcript(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT, symbol="BRCA2")
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.gene_symbol == "BRCA2"

    def test_hgvs_c_and_p_populated(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT, hgvs_c="c.100A>T", hgvs_p="p.Asn34Tyr")
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.hgvs_c == "c.100A>T"
        assert v.hgvs_p == "p.Asn34Tyr"

    def test_zygosity_heterozygous(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq, gt="0/1")]))
        v = parse_vep_output(str(vcf))[0]
        assert v.zygosity == "heterozygous"

    def test_zygosity_homozygous_alt(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq, gt="1/1")]))
        v = parse_vep_output(str(vcf))[0]
        assert v.zygosity == "homozygous_alt"

    def test_filter_status_preserved(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq, filter_val="LowQual")]))
        v = parse_vep_output(str(vcf))[0]
        assert v.filter_status == "LowQual"

    def test_gnomad_af_uses_max_af_when_present(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT, gnomad_af="0.001", max_af="0.005")
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.gnomad_af == pytest.approx(0.005)

    def test_in_silico_scores_on_variant(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT, cadd_phred="30.0", revel="0.9", am_path="0.85")
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.cadd_phred == pytest.approx(30.0)
        assert v.revel_score == pytest.approx(0.9)
        assert v.alphamissense_score == pytest.approx(0.85)

    def test_transcripts_list_preserved(self, tmp_path):
        fmt = FULL_CSQ_FORMAT
        e1 = _csq_entry(fmt, feature="TX1", canonical="YES")
        e2 = _csq_entry(fmt, feature="TX2", canonical="")
        csq = f"{e1},{e2}"
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert len(v.transcripts) == 2

    def test_canonical_selected_as_primary(self, tmp_path):
        fmt = FULL_CSQ_FORMAT
        noncanon = _csq_entry(fmt, feature="TX_nc", canonical="",
                              consequence="stop_gained")
        canon = _csq_entry(fmt, feature="TX_c", canonical="YES",
                           consequence="missense_variant")
        csq = f"{noncanon},{canon}"
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq=csq)]))
        v = parse_vep_output(str(vcf))[0]
        assert v.primary_transcript.transcript_id == "TX_c"

    def test_multiple_variants_parsed(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        records = [
            _record(chrom="chr1", pos=100, csq=csq),
            _record(chrom="chr1", pos=200, csq=csq),
            _record(chrom="chr2", pos=300, csq=csq),
        ]
        vcf = _write_vcf(tmp_path, _vcf_lines(records))
        variants = parse_vep_output(str(vcf))
        assert len(variants) == 3

    def test_variant_no_csq_entry_has_empty_transcripts(self, tmp_path):
        vcf = _write_vcf(tmp_path, _vcf_lines([_record(csq="")]))
        v = parse_vep_output(str(vcf))[0]
        assert v.transcripts == []
        assert v.primary_transcript is None
        assert v.gene_symbol == ""

    def test_clinvar_from_csq_field(self, tmp_path):
        # VEP --custom embeds ClinVar inside the CSQ per-transcript field,
        # not as a standalone INFO key.  The parser must read from CSQ raw dict.
        fmt = FULL_CSQ_FORMAT + "|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN"
        csq = _csq_entry(
            fmt,
            ClinVar_CLNSIG="Pathogenic",
            ClinVar_CLNREVSTAT="reviewed_by_expert_panel",
            ClinVar_CLNDN="Hereditary_breast_cancer",
        )
        record_line = (
            f"chr1\t100\t.\tA\tT\t50\tPASS\tCSQ={csq}\tGT:DP\t0/1:30"
        )
        vcf = _write_vcf(tmp_path, _vcf_lines([record_line], csq_format=fmt))
        v = parse_vep_output(str(vcf))[0]
        assert v.clinvar_significance == "Pathogenic"
        assert v.clinvar_review_status == "reviewed_by_expert_panel"
        assert v.clinvar_disease == "Hereditary_breast_cancer"

    def test_no_sample_column_zygosity_is_none(self, tmp_path):
        fmt = FULL_CSQ_FORMAT
        csq = _csq_entry(fmt)
        record_line = f"chr1\t100\t.\tA\tT\t50\tPASS\tCSQ={csq}"
        vcf = _write_vcf(tmp_path, _vcf_lines([record_line], has_sample=False))
        v = parse_vep_output(str(vcf))[0]
        assert v.zygosity is None


# ---------------------------------------------------------------------------
# parse_vep_output — bgzipped VCF
# ---------------------------------------------------------------------------

class TestParseVepOutputGzVcf:
    def test_gz_vcf_parsed(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        content = _vcf_lines([_record(csq=csq)])
        vcf = _write_gz_vcf(tmp_path, content)
        variants = parse_vep_output(str(vcf))
        assert len(variants) == 1
        assert variants[0].gene_symbol == "BRCA1"


# ---------------------------------------------------------------------------
# parse_vep_output — error conditions
# ---------------------------------------------------------------------------

class TestParseVepOutputErrors:
    def test_no_csq_header_raises(self, tmp_path):
        content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        vcf = _write_vcf(tmp_path, content)
        with pytest.raises(ValueError, match="CSQ annotation not found"):
            parse_vep_output(str(vcf))

    def test_multiallelic_raises(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        record_line = f"chr1\t100\t.\tA\tT,G\t50\tPASS\tCSQ={csq}\tGT:DP\t0/1:30"
        vcf = _write_vcf(tmp_path, _vcf_lines([record_line]))
        with pytest.raises(ValueError, match="Multi-allelic"):
            parse_vep_output(str(vcf))

    def test_multiallelic_error_contains_normalization_hint(self, tmp_path):
        csq = _csq_entry(FULL_CSQ_FORMAT)
        record_line = f"chr1\t100\t.\tA\tT,G\t50\tPASS\tCSQ={csq}\tGT:DP\t0/1:30"
        vcf = _write_vcf(tmp_path, _vcf_lines([record_line]))
        with pytest.raises(ValueError, match="bcftools norm"):
            parse_vep_output(str(vcf))

    def test_empty_vcf_raises(self, tmp_path):
        vcf = _write_vcf(tmp_path, "##fileformat=VCFv4.2\n")
        with pytest.raises(ValueError, match="#CHROM"):
            parse_vep_output(str(vcf))

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            parse_vep_output(str(tmp_path / "nonexistent.vcf"))
