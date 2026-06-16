"""
parser.py — Parse VEP-annotated VCF output into the internal Variant data model.

Parsing is done with Python's standard library only (no pysam dependency).

Design decisions (v1):
  - CSQ field names are read dynamically from the ##INFO=<ID=CSQ,...> header.
    Field order is never assumed — it varies by VEP version and active plugins.
  - All transcript annotations per variant are preserved in Variant.transcripts.
  - One primary transcript is selected per variant using this priority:
      1. Canonical transcript (CANONICAL=YES), most severe consequence if tied.
      2. Most severe consequence across all transcripts.
      3. First transcript annotation if all severities are equal.
  - Multi-allelic records (comma in ALT) raise ValueError with normalization
    instructions. They are never silently skipped.
  - Missing optional annotation fields (CADD, REVEL, SpliceAI, AlphaMissense,
    dbNSFP-derived scores, ClinVar) → None. No exception raised.
  - Handles both plain VCF (.vcf) and bgzipped VCF (.vcf.gz / .vcf.bgz).

Output of this module feeds:
  analysis/filter_engine.py  → filtering
  analysis/acmg/classifier.py → ACMG/AMP classification
"""

import gzip
import logging
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Consequence severity ranking
# Source: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
# Lower rank = more severe. Used to select the primary transcript.
# ---------------------------------------------------------------------------

CONSEQUENCE_SEVERITY: Dict[str, int] = {
    "transcript_ablation":                    1,
    "splice_acceptor_variant":                2,
    "splice_donor_variant":                   3,
    "stop_gained":                            4,
    "frameshift_variant":                     5,
    "stop_lost":                              6,
    "start_lost":                             7,
    "transcript_amplification":               8,
    "feature_elongation":                     9,
    "feature_truncation":                    10,
    "inframe_insertion":                     11,
    "inframe_deletion":                      12,
    "missense_variant":                      13,
    "protein_altering_variant":              14,
    "splice_donor_5th_base_variant":         15,
    "splice_region_variant":                 16,
    "splice_donor_region_variant":           17,
    "splice_polypyrimidine_tract_variant":   18,
    "incomplete_terminal_codon_variant":     19,
    "start_retained_variant":               20,
    "stop_retained_variant":                21,
    "synonymous_variant":                    22,
    "coding_sequence_variant":               23,
    "mature_miRNA_variant":                  24,
    "5_prime_UTR_variant":                   25,
    "3_prime_UTR_variant":                   26,
    "non_coding_transcript_exon_variant":    27,
    "intron_variant":                        28,
    "NMD_transcript_variant":               29,
    "non_coding_transcript_variant":         30,
    "coding_transcript_variant":             31,
    "upstream_gene_variant":                 32,
    "downstream_gene_variant":               33,
    "TFBS_ablation":                         34,
    "TFBS_amplification":                    35,
    "TF_binding_site_variant":               36,
    "regulatory_region_ablation":            37,
    "regulatory_region_amplification":       38,
    "regulatory_region_variant":             39,
    "intergenic_variant":                    40,
    "sequence_variant":                      41,
}

_DEFAULT_SEVERITY = 9999

# Regex for VEP's combined "score(prediction)" format used by --sift b / --polyphen b
_SCORE_PRED_RE = re.compile(r"^([\d.]+(?:[eE][+-]?\d+)?)\(([^)]+)\)$")

# Regex for extracting the CSQ Format string from the ##INFO header line
_CSQ_FORMAT_RE = re.compile(
    r'##INFO=<ID=CSQ[^>]*Format:\s*([^"]+)"',
    re.IGNORECASE,
)

# SpliceAI plugin delta-score field names (max of all four = overall score)
_SPLICEAI_DS_FIELDS = (
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
)

# REVEL field names (plugin name differs between REVEL plugin and dbNSFP)
_REVEL_FIELD_NAMES = ("REVEL", "REVEL_score")

# gnomAD AF field name candidates (vary by VEP version and flags used)
_GNOMAD_AF_FIELDS = ("MAX_AF", "gnomAD_AF", "gnomADe_AF", "gnomADg_AF")

# CADD phred score field name candidates
_CADD_PHRED_FIELDS = ("CADD_PHRED", "CADD_phred")


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class TranscriptAnnotation:
    """
    One VEP CSQ entry: the annotation for one transcript overlapping one allele.

    Optional score fields default to None when the corresponding plugin or
    resource was not available at annotation time.
    """
    # Core per-transcript fields
    allele: str = ""
    consequence: str = ""                   # raw string, e.g. "missense_variant&splice_region_variant"
    consequences: List[str] = field(default_factory=list)   # split on &
    impact: str = ""                        # HIGH / MODERATE / LOW / MODIFIER
    gene_symbol: str = ""
    gene_id: str = ""
    feature_type: str = ""
    transcript_id: str = ""
    biotype: str = ""
    exon: Optional[str] = None
    intron: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    cdna_position: Optional[str] = None
    cds_position: Optional[str] = None
    protein_position: Optional[str] = None
    amino_acids: Optional[str] = None
    codons: Optional[str] = None
    is_canonical: bool = False

    # Built-in VEP predictions (--sift b / --polyphen b)
    sift_score: Optional[float] = None
    sift_pred: Optional[str] = None
    polyphen_score: Optional[float] = None
    polyphen_pred: Optional[str] = None

    # Population frequency from VEP cache
    af_gnomad: Optional[float] = None      # gnomAD exomes overall AF
    max_af: Optional[float] = None         # maximum AF across all populations

    # In-silico scores from optional plugins
    cadd_phred: Optional[float] = None
    revel_score: Optional[float] = None
    spliceai_max_ds: Optional[float] = None        # max of the four delta scores
    alphamissense_score: Optional[float] = None

    # LOFTEE loss-of-function classification (from LoF plugin)
    loftee_lof: Optional[str] = None       # HC (high confidence) or LC (low confidence)
    loftee_filter: Optional[str] = None    # reason for LC, e.g. END_TRUNC, INCOMPLETE_CDS
    loftee_flags: Optional[str] = None     # additional flags, e.g. NON_CAN_SPLICE

    # Derived severity (lower = more severe; used for primary transcript selection)
    severity_rank: int = _DEFAULT_SEVERITY

    # All raw CSQ field values for this entry
    raw: Dict[str, str] = field(default_factory=dict)


@dataclass
class Variant:
    """
    Internal representation of one annotated, single-allele variant.

    Flat fields (gene_symbol, transcript_id, etc.) mirror the primary
    transcript for downstream convenience. Full transcript detail lives in
    `transcripts`.

    Fields populated by later pipeline stages (filter_engine, classifier,
    ranker) are present here with defaults so the dataclass is complete from
    the start.
    """
    # Genomic position
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    variant_id: str = ""

    # VCF quality / filter
    qual: Optional[str] = None
    filter_status: str = ""

    # Zygosity (from FORMAT/sample GT field; None if no sample column)
    zygosity: Optional[str] = None

    # All transcript annotations from CSQ (preserved in full)
    transcripts: List[TranscriptAnnotation] = field(default_factory=list)
    primary_transcript: Optional[TranscriptAnnotation] = None

    # Flat fields populated from primary_transcript
    gene_symbol: str = ""
    transcript_id: str = ""
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: str = ""

    # Variant-level scores (from primary transcript)
    gnomad_af: Optional[float] = None          # max_af if present, else af_gnomad

    cadd_phred: Optional[float] = None
    revel_score: Optional[float] = None
    spliceai_max_score: Optional[float] = None
    alphamissense_score: Optional[float] = None

    sift_score: Optional[float] = None
    sift_pred: Optional[str] = None
    polyphen_score: Optional[float] = None
    polyphen_pred: Optional[str] = None

    # LOFTEE (from primary transcript CSQ LoF fields)
    loftee_lof: Optional[str] = None          # HC or LC
    loftee_filter: Optional[str] = None       # e.g. END_TRUNC
    loftee_flags: Optional[str] = None        # e.g. NON_CAN_SPLICE

    # ClinGen gene disease validity (populated by analysis/clingen.py)
    clingen_validity: Optional[str] = None    # Definitive, Strong, Moderate, Limited, Disputed, Refuted

    # ClinVar (from INFO field via --custom, not from CSQ)
    clinvar_significance: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    clinvar_disease: Optional[str] = None

    # ACMG/AMP classification (populated by acmg/classifier.py)
    acmg_evidence_codes: List[str] = field(default_factory=list)
    acmg_classification: Optional[str] = None
    classification_confidence: Optional[str] = None  # High / Moderate / Low

    # Ordering context (populated by ranker.py)
    in_gene_list: bool = False

    # Filtering (populated by filter_engine.py)
    excluded: bool = False
    exclusion_reason: Optional[str] = None

    # QC metrics extracted from FORMAT/sample columns
    depth: Optional[int] = None               # DP — total read depth
    genotype_quality: Optional[int] = None    # GQ — genotype quality score
    allele_depth: Optional[str] = None        # AD — "ref_count,alt_count"
    vaf: Optional[float] = None              # variant allele fraction = alt_AD / (ref_AD + alt_AD)

    # Reviewer flags (populated by analysis/reviewer_flags.py after classification)
    reviewer_flags: List[str] = field(default_factory=list)

    # Gene relevance tag (populated by filter_engine.py)
    gene_relevance: str = ""

    # Soft deprioritization (distinct from hard exclusion)
    deprioritized: bool = False
    deprioritization_reason: Optional[str] = None

    # Audit trail: populated by filter_engine.py and classifier.py so every
    # pipeline decision is traceable without reading raw logs.
    decision_trace: dict = field(default_factory=dict)

    # ACMG rules that threw an exception during evaluation (populated by classifier.py).
    # An empty list means all rules completed normally (not that all returned False).
    acmg_skipped_rules: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def parse_vep_output(vep_vcf_path: str) -> List[Variant]:
    """
    Parse a VEP-annotated VCF file into a list of Variant objects.

    Parameters
    ----------
    vep_vcf_path: Path to a VEP-annotated VCF or VCF.GZ file.

    Returns
    -------
    List of Variant objects, one per biallelic VCF record (in order).

    Raises
    ------
    ValueError: If the VCF has no CSQ header, if a record is multi-allelic,
                or if the file is fundamentally malformed.
    """
    header_lines: List[str] = []
    column_line: Optional[str] = None
    csq_fields: Optional[List[str]] = None
    has_sample: bool = False
    variants: List[Variant] = []
    record_count = 0

    with _open_vcf(vep_vcf_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith("##"):
                header_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                column_line = line
                csq_fields = _extract_csq_fields(header_lines)
                has_sample = _has_sample_column(column_line)
                continue

            # Data record
            if csq_fields is None:
                raise ValueError(
                    f"VCF data record encountered before the #CHROM header line "
                    f"in {vep_vcf_path}.\n"
                    f"  The file may be malformed or truncated."
                )

            record_count += 1
            cols = line.split("\t")
            try:
                variant = _parse_record(cols, csq_fields, has_sample)
                variants.append(variant)
            except ValueError:
                raise   # multi-allelic and structural errors stop execution
            except Exception as exc:
                chrom = cols[0] if cols else "?"
                pos = cols[1] if len(cols) > 1 else "?"
                logger.warning(
                    "Skipping malformed VCF record at %s:%s — %s",
                    chrom, pos, exc,
                )

    if column_line is None:
        raise ValueError(
            f"No #CHROM header line found in {vep_vcf_path}.\n"
            f"  The file may be empty or not a valid VCF."
        )

    logger.info(
        "Parsed %d / %d records from %s",
        len(variants), record_count, vep_vcf_path,
    )
    return variants


# ---------------------------------------------------------------------------
# File handling
# ---------------------------------------------------------------------------

def _open_vcf(vcf_path: str):
    """Open a plain or bgzipped VCF file for text reading."""
    lower = vcf_path.lower()
    if lower.endswith(".gz") or lower.endswith(".bgz"):
        return gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace")
    return open(vcf_path, "r", encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def _extract_csq_fields(header_lines: List[str]) -> List[str]:
    """
    Extract the ordered CSQ field names from the VCF ##INFO header.

    Raises ValueError if no CSQ header line is found.
    """
    for line in header_lines:
        m = _CSQ_FORMAT_RE.search(line)
        if m:
            fields = [f.strip() for f in m.group(1).split("|")]
            logger.debug("CSQ field order (%d fields): %s", len(fields), fields)
            return fields

    raise ValueError(
        "VEP CSQ annotation not found in VCF header.\n"
        "  This VCF does not appear to have been annotated with VEP (--vcf mode).\n"
        "  Run the annotation step before parsing:\n"
        "    python cli.py --annotate ...\n"
        "  If VEP was run, confirm the --vcf flag was used and the correct "
        "output file is being parsed."
    )


def _has_sample_column(column_line: str) -> bool:
    """Return True if the VCF has at least one sample column (col index >= 9)."""
    return len(column_line.split("\t")) >= 10


# ---------------------------------------------------------------------------
# Record parsing
# ---------------------------------------------------------------------------

def _parse_record(
    cols: List[str],
    csq_fields: List[str],
    has_sample: bool,
) -> Variant:
    """
    Parse one VCF data record into a Variant.

    Raises ValueError for multi-allelic records and structurally invalid lines.
    """
    if len(cols) < 8:
        raise ValueError(
            f"VCF record has only {len(cols)} columns (minimum 8 required). "
            f"The record may be malformed."
        )

    chrom      = cols[0]
    pos_str    = cols[1]
    variant_id = cols[2]
    ref        = cols[3]
    alt        = cols[4]
    qual_str   = cols[5]
    filter_val = cols[6]
    info_str   = cols[7]

    try:
        pos = int(pos_str)
    except ValueError:
        raise ValueError(f"Non-integer POS field: {pos_str!r} at {chrom}.")

    # Hard stop on multi-allelic records
    _check_multiallelic(chrom, pos, ref, alt)

    info_dict = _parse_info_str(info_str)
    csq_value = info_dict.get("CSQ", "")
    transcripts = _parse_csq_value(csq_value, csq_fields)
    primary = _select_primary_transcript(transcripts)
    # ClinVar is embedded in the per-transcript CSQ field (VEP --custom injects it
    # there), not as a standalone INFO key. Read from the primary transcript's raw dict.
    clinvar_raw = primary.raw if primary else {}
    clinvar_sig, clinvar_revstat, clinvar_disease = _extract_clinvar(clinvar_raw)

    zygosity = None
    if has_sample and len(cols) >= 10:
        zygosity = _extract_zygosity(cols[8], cols[9])

    depth = genotype_quality = allele_depth = vaf = None
    if has_sample and len(cols) >= 10:
        qc = _extract_qc_metrics(cols[8], cols[9])
        depth            = qc["depth"]
        genotype_quality = qc["genotype_quality"]
        allele_depth     = qc["allele_depth"]
        vaf              = qc["vaf"]

    v = Variant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_id=variant_id if variant_id != "." else "",
        qual=qual_str if qual_str != "." else None,
        filter_status=filter_val,
        zygosity=zygosity,
        transcripts=transcripts,
        primary_transcript=primary,
        clinvar_significance=clinvar_sig,
        clinvar_review_status=clinvar_revstat,
        clinvar_disease=clinvar_disease,
        depth=depth,
        genotype_quality=genotype_quality,
        allele_depth=allele_depth,
        vaf=vaf,
    )

    if primary:
        v.gene_symbol  = primary.gene_symbol
        v.transcript_id = primary.transcript_id
        v.hgvs_c       = primary.hgvs_c
        v.hgvs_p       = primary.hgvs_p
        v.consequence  = primary.consequence
        # Use MAX_AF (most conservative) when available; fall back to gnomAD_AF
        v.gnomad_af          = primary.max_af if primary.max_af is not None else primary.af_gnomad
        v.cadd_phred         = primary.cadd_phred
        v.revel_score        = primary.revel_score
        # SpliceAI plugin writes delta scores to INFO fields (not CSQ) when
        # VEP runs with --vcf output. Read from info_dict first; fall back to
        # CSQ-embedded value for annotated VCFs produced by other toolchains.
        spliceai_info = _extract_spliceai_max(info_dict)
        v.spliceai_max_score = spliceai_info if spliceai_info is not None else primary.spliceai_max_ds
        v.alphamissense_score = primary.alphamissense_score
        v.sift_score         = primary.sift_score
        v.sift_pred          = primary.sift_pred
        v.polyphen_score     = primary.polyphen_score
        v.polyphen_pred      = primary.polyphen_pred
        v.loftee_lof         = primary.loftee_lof
        v.loftee_filter      = primary.loftee_filter
        v.loftee_flags       = primary.loftee_flags

    return v


def _check_multiallelic(chrom: str, pos: int, ref: str, alt: str) -> None:
    """
    Raise a clear ValueError if the record has more than one ALT allele.

    Multi-allelic sites must be split before this pipeline can process them.
    """
    if "," in alt:
        raise ValueError(
            f"Multi-allelic variant at {chrom}:{pos} "
            f"(REF={ref}, ALT={alt}).\n"
            f"\n"
            f"  This pipeline requires biallelic (decomposed) records.\n"
            f"  Split multi-allelic sites before running analysis:\n"
            f"\n"
            f"    bcftools norm --multiallelics -both \\\n"
            f"      --fasta-ref /path/to/GRCh38.fa \\\n"
            f"      input.vcf.gz \\\n"
            f"      | bgzip > input_split.vcf.gz\n"
            f"    tabix -p vcf input_split.vcf.gz\n"
            f"\n"
            f"  Then re-run with input_split.vcf.gz."
        )


def _parse_info_str(info_str: str) -> Dict[str, str]:
    """
    Parse a VCF INFO field string into a flat dict.

    Splitting on ';' is safe because VEP URL-encodes semicolons within field
    values, and the CSQ value uses ',' (not ';') as its internal separator.
    """
    result: Dict[str, str] = {}
    if not info_str or info_str == ".":
        return result
    for token in info_str.split(";"):
        if "=" in token:
            key, val = token.split("=", 1)   # split on first '=' only
            result[key.strip()] = val
        elif token.strip():
            result[token.strip()] = ""        # flag-style INFO field
    return result


# ---------------------------------------------------------------------------
# CSQ parsing
# ---------------------------------------------------------------------------

def _parse_csq_value(
    csq_value: str,
    csq_fields: List[str],
) -> List[TranscriptAnnotation]:
    """
    Parse the full CSQ INFO field value into a list of TranscriptAnnotation objects.

    Each transcript entry is comma-separated; fields within are pipe-separated.
    Malformed individual entries are skipped with a warning — they do not stop parsing.
    """
    if not csq_value:
        return []

    transcripts: List[TranscriptAnnotation] = []
    for raw_entry in csq_value.split(","):
        values = raw_entry.split("|")
        if len(values) != len(csq_fields):
            logger.warning(
                "CSQ entry has %d fields but header defines %d — skipping entry: %s...",
                len(values), len(csq_fields), raw_entry[:80],
            )
            continue
        try:
            ta = _parse_single_csq_entry(csq_fields, values)
            transcripts.append(ta)
        except Exception as exc:
            logger.warning(
                "Failed to parse CSQ entry (skipping): %s... Error: %s",
                raw_entry[:80], exc,
            )

    return transcripts


def _parse_single_csq_entry(
    csq_fields: List[str],
    values: List[str],
) -> TranscriptAnnotation:
    """Build one TranscriptAnnotation from a single CSQ pipe-delimited entry."""
    raw = dict(zip(csq_fields, values))

    consequence_str  = raw.get("Consequence", "")
    consequences     = [c for c in consequence_str.split("&") if c]
    is_canonical     = raw.get("CANONICAL", "").upper() == "YES"

    sift_score, sift_pred       = _parse_score_pred(raw.get("SIFT", ""))
    polyphen_score, polyphen_pred = _parse_score_pred(raw.get("PolyPhen", ""))

    # gnomAD AF — try multiple field name candidates
    af_gnomad = None
    for fn in ("gnomAD_AF", "gnomADe_AF", "gnomADg_AF"):
        af_gnomad = _safe_float(raw.get(fn, ""))
        if af_gnomad is not None:
            break
    max_af = _safe_float(raw.get("MAX_AF", ""))

    # Optional in-silico scores
    cadd_phred        = _extract_cadd(raw)
    revel_score       = _extract_revel(raw)
    spliceai_max      = _extract_spliceai_max(raw)
    alphamissense     = _safe_float(raw.get("am_pathogenicity", ""))

    # LOFTEE — LoF, LoF_filter, LoF_flags are per-transcript CSQ fields
    loftee_lof    = _safe_str(raw.get("LoF", ""))
    loftee_filter = _safe_str(raw.get("LoF_filter", ""))
    loftee_flags  = _safe_str(raw.get("LoF_flags", ""))

    severity = _consequence_severity_rank(consequences)

    return TranscriptAnnotation(
        allele           = raw.get("Allele", ""),
        consequence      = consequence_str,
        consequences     = consequences,
        impact           = raw.get("IMPACT", ""),
        gene_symbol      = raw.get("SYMBOL", "") or raw.get("Gene", ""),
        gene_id          = raw.get("Gene", ""),
        feature_type     = raw.get("Feature_type", ""),
        transcript_id    = raw.get("Feature", ""),
        biotype          = raw.get("BIOTYPE", ""),
        exon             = _safe_str(raw.get("EXON", "")),
        intron           = _safe_str(raw.get("INTRON", "")),
        hgvs_c           = _safe_str(raw.get("HGVSc", "")),
        hgvs_p           = _safe_str(raw.get("HGVSp", "")),
        cdna_position    = _safe_str(raw.get("cDNA_position", "")),
        cds_position     = _safe_str(raw.get("CDS_position", "")),
        protein_position = _safe_str(raw.get("Protein_position", "")),
        amino_acids      = _safe_str(raw.get("Amino_acids", "")),
        codons           = _safe_str(raw.get("Codons", "")),
        is_canonical     = is_canonical,
        sift_score       = sift_score,
        sift_pred        = sift_pred,
        polyphen_score   = polyphen_score,
        polyphen_pred    = polyphen_pred,
        af_gnomad        = af_gnomad,
        max_af           = max_af,
        cadd_phred       = cadd_phred,
        revel_score      = revel_score,
        spliceai_max_ds  = spliceai_max,
        alphamissense_score = alphamissense,
        loftee_lof           = loftee_lof,
        loftee_filter        = loftee_filter,
        loftee_flags         = loftee_flags,
        severity_rank        = severity,
        raw                  = raw,
    )


# ---------------------------------------------------------------------------
# Primary transcript selection
# ---------------------------------------------------------------------------

def _select_primary_transcript(
    transcripts: List[TranscriptAnnotation],
) -> Optional[TranscriptAnnotation]:
    """
    Choose one primary transcript annotation per variant using this priority:
      1. Canonical transcript (CANONICAL=YES), most severe consequence if tied.
      2. Most severe consequence across all transcripts.
      3. First entry if severity is tied across all.
    """
    if not transcripts:
        return None

    canonical = [t for t in transcripts if t.is_canonical]
    if canonical:
        return min(canonical, key=lambda t: t.severity_rank)

    return min(transcripts, key=lambda t: t.severity_rank)


# ---------------------------------------------------------------------------
# Variant-level annotation extraction
# ---------------------------------------------------------------------------

def _extract_clinvar(
    info_dict: Dict[str, str],
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Extract ClinVar fields added to INFO via VEP --custom.

    VEP --custom names these as ClinVar_CLNSIG, ClinVar_CLNREVSTAT, ClinVar_CLNDN.
    Returns (significance, review_status, disease) — any may be None.
    """
    return (
        _safe_str(info_dict.get("ClinVar_CLNSIG", "")),
        _safe_str(info_dict.get("ClinVar_CLNREVSTAT", "")),
        _safe_str(info_dict.get("ClinVar_CLNDN", "")),
    )


def _extract_qc_metrics(format_str: str, sample_str: str) -> dict:
    """Extract DP, GQ, AD and derived VAF from FORMAT and sample columns."""
    result: dict = {
        "depth": None, "genotype_quality": None,
        "allele_depth": None, "vaf": None,
    }
    if not format_str or not sample_str:
        return result
    fmt_fields    = format_str.split(":")
    sample_fields = sample_str.split(":")
    field_map     = dict(zip(fmt_fields, sample_fields))

    dp_raw = field_map.get("DP", "")
    if dp_raw and dp_raw != ".":
        try:
            result["depth"] = int(dp_raw)
        except ValueError:
            pass

    gq_raw = field_map.get("GQ", "")
    if gq_raw and gq_raw != ".":
        try:
            result["genotype_quality"] = int(gq_raw)
        except ValueError:
            pass

    ad_raw = field_map.get("AD", "")
    if ad_raw and ad_raw != ".":
        result["allele_depth"] = ad_raw
        try:
            parts = [int(x) for x in ad_raw.split(",")]
            if len(parts) >= 2:
                total = sum(parts)
                if total > 0:
                    result["vaf"] = round(parts[1] / total, 4)
        except ValueError:
            pass

    return result


def _extract_zygosity(format_str: str, sample_str: str) -> Optional[str]:
    """
    Derive zygosity from FORMAT and sample columns.

    Returns one of: 'heterozygous', 'homozygous_alt', 'homozygous_ref',
    'unknown', or None if GT is absent.
    """
    if not format_str or not sample_str:
        return None

    fmt_fields    = format_str.split(":")
    sample_fields = sample_str.split(":")

    if "GT" not in fmt_fields:
        return None

    gt_idx = fmt_fields.index("GT")
    if gt_idx >= len(sample_fields):
        return None

    gt = sample_fields[gt_idx]
    # Normalise phased (|) to unphased (/) and split alleles
    alleles = gt.replace("|", "/").split("/")

    if not alleles or any(a == "." for a in alleles):
        return "unknown"

    unique = set(alleles)
    if unique == {"0"}:
        return "homozygous_ref"
    if len(unique) == 1:
        return "homozygous_alt"
    if "0" in unique:
        return "heterozygous"
    # 1/2 or similar — shouldn't reach here after multi-allelic check
    return "heterozygous"


# ---------------------------------------------------------------------------
# Score extraction helpers
# ---------------------------------------------------------------------------

def _extract_cadd(raw: Dict[str, str]) -> Optional[float]:
    """Extract CADD phred score from CSQ raw dict, trying multiple field names."""
    for fn in _CADD_PHRED_FIELDS:
        v = _safe_float(raw.get(fn, ""))
        if v is not None:
            return v
    return None


def _extract_revel(raw: Dict[str, str]) -> Optional[float]:
    """
    Extract REVEL score from CSQ raw dict.
    Handles both the REVEL plugin (single value) and dbNSFP (semicolon-delimited
    multiple values — takes the maximum, i.e. most damaging).
    """
    for fn in _REVEL_FIELD_NAMES:
        val = raw.get(fn, "")
        if not val or val in (".", "-", ""):
            continue
        parts  = [p.strip() for p in val.split(";") if p.strip()]
        scores = [_safe_float(p) for p in parts]
        valid  = [s for s in scores if s is not None]
        if valid:
            return max(valid)
    return None


def _extract_spliceai_max(raw: Dict[str, str]) -> Optional[float]:
    """Return the maximum SpliceAI delta score across all four delta-score fields."""
    scores = [_safe_float(raw.get(fn, "")) for fn in _SPLICEAI_DS_FIELDS]
    valid  = [s for s in scores if s is not None]
    return max(valid) if valid else None


def _first_nonempty(value: str) -> str:
    """
    dbNSFP returns semicolon-delimited values when multiple transcripts are present.
    Return the first non-empty, non-dot token; empty string if none found.
    """
    if not value:
        return ""
    for token in value.split(";"):
        t = token.strip()
        if t and t != ".":
            return t
    return ""


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _consequence_severity_rank(consequences: List[str]) -> int:
    """Return the severity rank of the most severe consequence in the list."""
    if not consequences:
        return _DEFAULT_SEVERITY
    return min(CONSEQUENCE_SEVERITY.get(c.strip(), _DEFAULT_SEVERITY) for c in consequences)


def _parse_score_pred(raw: str) -> Tuple[Optional[float], Optional[str]]:
    """
    Parse VEP's combined 'score(prediction)' format used by --sift b / --polyphen b.

    Examples: '0.02(deleterious)' → (0.02, 'deleterious')
              '0.95(probably_damaging)' → (0.95, 'probably_damaging')
              '' or '.' → (None, None)
              '0.02' → (0.02, None)   # score only, no prediction
    """
    if not raw or raw.strip() in ("", ".", "-"):
        return None, None
    m = _SCORE_PRED_RE.match(raw.strip())
    if m:
        return _safe_float(m.group(1)), m.group(2)
    return _safe_float(raw), None


def _safe_float(val: Optional[str]) -> Optional[float]:
    """Parse a string to float; return None for empty, '.', 'NA', or unparseable."""
    if val is None:
        return None
    v = str(val).strip()
    if v in ("", ".", "NA", "NaN", "-", "N/A"):
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _safe_str(val: Optional[str]) -> Optional[str]:
    """Return None for empty/placeholder strings; strip whitespace otherwise."""
    if val is None:
        return None
    v = val.strip()
    return None if v in ("", ".", "-") else v
