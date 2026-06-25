#!/usr/bin/env bash
# =============================================================================
# vep_filter.sh — VEP annotation + bcftools filtering pipeline
#
# Steps:
#   1. Annotate input VCF with VEP (all plugins)
#   2. Filter: remove common variants (gnomAD AF >= 1%)
#   3. subsets/ folder:
#        gene_list_canonical.vcf.gz    — variants in any gene from the gene list (canonical transcript)
#        gene_list_noncanonical.vcf.gz — same genes, non-canonical transcripts only
#        loftee_hc.vcf.gz    — LOFTEE high-confidence LoF only
#        high_cadd.vcf.gz    — CADD phred >= 20
#        spliceai.vcf.gz     — SpliceAI delta score >= 0.5
#        clinvar_plp.vcf.gz  — ClinVar Pathogenic / Likely Pathogenic / Conflicting
#
# Usage:
#   bash vep_filter.sh --input sample.vcf.gz --out-dir results/
#
# Optional:
#   --gene-list  /path/to/genes.txt   (default: config/gene_list.txt)
#   --skip-vep                         skip VEP if annotated.vcf.gz already exists
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

INPUT_VCF=""
OUT_DIR=""
GENE_LIST="${SCRIPT_DIR}/config/gene_list.txt"
SKIP_VEP=false

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT_VCF="$2";  shift 2 ;;
        --out-dir)   OUT_DIR="$2";    shift 2 ;;
        --gene-list) GENE_LIST="$2";  shift 2 ;;
        --skip-vep)  SKIP_VEP=true;   shift ;;
        -h|--help)
            sed -n '2,22p' "$0" | tr -d '#'
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ -z "$INPUT_VCF" ]]; then
    echo "ERROR: --input is required." >&2; exit 1
fi
if [[ -z "$OUT_DIR" ]]; then
    echo "ERROR: --out-dir is required." >&2; exit 1
fi
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" >&2; exit 1
fi
if [[ ! -f "$GENE_LIST" ]]; then
    echo "ERROR: gene list not found: $GENE_LIST" >&2; exit 1
fi

SUBSETS_DIR="${OUT_DIR}/subsets"
mkdir -p "$OUT_DIR" "$SUBSETS_DIR"

# Normalise gene list to uppercase — bcftools =@ is case-sensitive but VEP
# emits HGNC symbols in uppercase, so a user-supplied list with mixed case
# (e.g. "Brca1") would silently miss matches without this step.
GENE_LIST_UC="$(mktemp)"
trap 'rm -f "$GENE_LIST_UC"' EXIT
tr '[:lower:]' '[:upper:]' < "$GENE_LIST" > "$GENE_LIST_UC"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

VEP=/home/admin7995/ensembl-vep/vep
CACHE=/home/admin7995/.vep
PLUGINS=/home/admin7995/.vep/Plugins
FASTA=/mnt/server_share/Reference_Material/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa

CADD_SNV=/mnt/server_share/sara/whole_genome_SNVs.tsv.gz
CADD_INDEL=/mnt/server_share/sara/gnomad.genomes.r4.0.indel.tsv.gz
REVEL=/mnt/server_share/Reference_Material/annotation/revel_all_chromosomes.tsv.gz
SPLICEAI_SNV=/mnt/server_share/Reference_Material/annotation/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz
SPLICEAI_INDEL=/mnt/server_share/Reference_Material/annotation/spliceai_scores.raw.indel.hg38.vcf.gz
ALPHAMISSENSE=/mnt/server_share/Reference_Material/annotation/AlphaMissense_hg38.tsv.gz
CLINVAR=/mnt/server_share/Reference_Material/annotation/clinvar.chr.vcf.gz
LOFTEE_DIR=/home/admin7995/.vep/loftee_grch38
LOFTEE_HUMAN=/mnt/server_share/Reference_Material/annotation/loftee_GRCh38/human_ancestor.fa.gz
LOFTEE_SQL=/mnt/server_share/Reference_Material/annotation/loftee_GRCh38/loftee.sql
LOFTEE_GERP=/mnt/server_share/Reference_Material/annotation/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------

ANNOTATED="${OUT_DIR}/annotated.vcf.gz"
RARE="${OUT_DIR}/rare_annotated.1pc.vcf.gz"

GENE_LIST_CANONICAL="${SUBSETS_DIR}/gene_list_canonical.vcf.gz"
GENE_LIST_NONCANON="${SUBSETS_DIR}/gene_list_noncanonical.vcf.gz"
LOFTEE_HC="${SUBSETS_DIR}/loftee_hc.vcf.gz"
HIGH_CADD="${SUBSETS_DIR}/high_cadd.vcf.gz"
SPLICEAI="${SUBSETS_DIR}/spliceai.vcf.gz"
CLINVAR_PLP="${SUBSETS_DIR}/clinvar_plp.vcf.gz"

log()   { echo "[$(date '+%H:%M:%S')] $*"; }
stage() { echo "##STAGE## $1/6 $2"; }
count() { bcftools view -H "$1" | wc -l; }

# Run a command in background, store PID and label for error checking
declare -a _bg_pids=()
declare -a _bg_labels=()
bg() {
    local label="$1"; shift
    "$@" &
    _bg_pids+=($!)
    _bg_labels+=("$label")
}
wait_all() {
    local failed=0
    for i in "${!_bg_pids[@]}"; do
        if ! wait "${_bg_pids[$i]}"; then
            log "ERROR: ${_bg_labels[$i]} failed"
            failed=1
        fi
    done
    _bg_pids=(); _bg_labels=()
    [[ $failed -eq 0 ]]
}

# ---------------------------------------------------------------------------
# Step 1 — VEP annotation
# ---------------------------------------------------------------------------

stage 1 "Annotating variants with VEP"
if $SKIP_VEP; then
    log "Step 1/6: skipping VEP — using input as annotated VCF"
    ANNOTATED="$INPUT_VCF"
else
    log "Step 1/6: VEP annotation → $ANNOTATED"
    "$VEP" \
        --input_file  "$INPUT_VCF" \
        --output_file "$ANNOTATED" --vcf \
        --cache --offline --assembly GRCh38 \
        --dir_cache   "$CACHE" \
        --dir_plugins "$PLUGINS" \
        --fasta       "$FASTA" \
        --fork 64 --buffer_size 20000 \
        --symbol --biotype --numbers --canonical --transcript_version \
        --hgvs --hgvsg \
        --af --af_gnomad --max_af \
        --sift b --polyphen b \
        --no_stats --no_progress --force_overwrite \
        --plugin "CADD,${CADD_SNV},${CADD_INDEL}" \
        --plugin "REVEL,${REVEL}" \
        --plugin "SpliceAI,snv=${SPLICEAI_SNV},indel=${SPLICEAI_INDEL}" \
        --plugin "AlphaMissense,file=${ALPHAMISSENSE}" \
        --plugin "LoF,loftee_path:${LOFTEE_DIR},human_ancestor_fa:${LOFTEE_HUMAN},conservation_file:${LOFTEE_SQL},gerp_bigwig:${LOFTEE_GERP}" \
        --custom "${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" \
        2> >(grep -v "^Smartmatch is deprecated" >&2)
    tabix "$ANNOTATED"
    log "Step 1/6 done."
fi

# ---------------------------------------------------------------------------
# Step 2 — Remove common variants (AF >= 1%)
# ---------------------------------------------------------------------------

stage 2 "Filtering common variants (AF ≥ 1%)"
log "Step 2/6: AF filter (<1%) + PASS only → $RARE"
bcftools view -f "PASS,." "$ANNOTATED" \
    | bcftools +split-vep \
        -c MAX_AF:Float \
        -e 'MAX_AF >= 0.01' \
        -Oz -o "$RARE"
tabix "$RARE"
RARE_N=$(count "$RARE")
log "Step 2/6 done — ${RARE_N} rare variants retained."

# ---------------------------------------------------------------------------
# Steps 3–6 — All 6 subsets run in parallel (all read from $RARE)
# ---------------------------------------------------------------------------

stage 3 "Producing all subsets in parallel"
log "Steps 3-6: launching all 6 subset filters in parallel..."

bg "gene_list_canonical" \
    bcftools +split-vep "$RARE" \
        -c SYMBOL:String,CANONICAL:String \
        -i "CANONICAL='YES' && SYMBOL=@${GENE_LIST_UC}" \
        -Oz -o "$GENE_LIST_CANONICAL"

bg "gene_list_noncanonical" \
    bcftools +split-vep "$RARE" \
        -c SYMBOL:String,CANONICAL:String \
        -i "CANONICAL!='YES' && SYMBOL=@${GENE_LIST_UC}" \
        -Oz -o "$GENE_LIST_NONCANON"

bg "loftee_hc" \
    bcftools +split-vep "$RARE" \
        -c LoF:String,CANONICAL:String \
        -i 'CANONICAL="YES" && LoF="HC"' \
        -Oz -o "$LOFTEE_HC"

bg "high_cadd" \
    bcftools +split-vep "$RARE" \
        -c CADD_PHRED:Float,CANONICAL:String \
        -i 'CANONICAL="YES" && CADD_PHRED >= 20' \
        -Oz -o "$HIGH_CADD"

bg "spliceai" \
    bcftools +split-vep "$RARE" \
        -c ClinVar_CLNSIG:String \
        -i 'SpliceAI_pred_DS_AG>=0.5 || SpliceAI_pred_DS_AL>=0.5 || SpliceAI_pred_DS_DG>=0.5 || SpliceAI_pred_DS_DL>=0.5' \
        -Oz -o "$SPLICEAI"

bg "clinvar_plp" \
    bcftools +split-vep "$RARE" \
        -c ClinVar_CLNSIG:String \
        -i 'ClinVar_CLNSIG~"[Pp]athogenic"' \
        -Oz -o "$CLINVAR_PLP"

wait_all
log "All 6 subsets complete — indexing in parallel..."

bg "tabix gene_list_canonical"    tabix "$GENE_LIST_CANONICAL"
bg "tabix gene_list_noncanonical" tabix "$GENE_LIST_NONCANON"
bg "tabix loftee_hc"              tabix "$LOFTEE_HC"
bg "tabix high_cadd"              tabix "$HIGH_CADD"
bg "tabix spliceai"               tabix "$SPLICEAI"
bg "tabix clinvar_plp"            tabix "$CLINVAR_PLP"
wait_all
log "Indexing complete."

stage 4 "Counting subset results"
stage 5 "Counting subset results"
stage 6 "Counting subset results"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

log ""
log "=== Done ==="
log "  Output dir:  $OUT_DIR"
log "  annotated.vcf.gz       — full VEP output"
log "  rare_annotated.1pc.vcf.gz — ${RARE_N} variants (AF < 1%)"
log ""
log "  subsets/"
log "    gene_list_canonical.vcf.gz    — $(count "$GENE_LIST_CANONICAL") variants in gene list (canonical)"
log "    gene_list_noncanonical.vcf.gz — $(count "$GENE_LIST_NONCANON") variants in gene list (non-canonical)"
log "    loftee_hc.vcf.gz     — $(count "$LOFTEE_HC") LOFTEE high-confidence LoF"
log "    high_cadd.vcf.gz     — $(count "$HIGH_CADD") variants CADD >= 20"
log "    spliceai.vcf.gz      — $(count "$SPLICEAI") splice-altering (SpliceAI >= 0.5)"
log "    clinvar_plp.vcf.gz   — $(count "$CLINVAR_PLP") ClinVar Pathogenic / Likely Pathogenic / Conflicting"
