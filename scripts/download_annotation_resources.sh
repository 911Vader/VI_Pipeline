#!/usr/bin/env bash
# =============================================================================
# download_annotation_resources.sh
#
# Downloads and prepares annotation resources required for full ACMG scoring:
#   - ClinVar  (PS3 — strong pathogenic evidence from expert-reviewed variants)
#   - REVEL    (PP3/BP4 — missense pathogenicity scores)
#   - SpliceAI (PP3/BP4 — splicing impact scores)
#
# After running this script, re-run setup:
#   python /home/admin7995/VI_Pipeline/cli.py --setup --force
#
# Storage required:
#   ClinVar      :  ~120 MB
#   REVEL        :  ~7   GB  (compressed)
#   SpliceAI SNV :  ~28  GB  (Ensembl MANE v1.4, Ensembl FTP — public)
#   SpliceAI indel: ~2   GB  (Illumina BaseSpace — free account required)
#   Total        :  ~37  GB
#
# All files are placed in:
#   /mnt/server_share/Reference_Material/annotation/
# =============================================================================

set -euo pipefail

DEST="/mnt/server_share/Reference_Material/annotation"
TABIX=$(which tabix 2>/dev/null || echo "/usr/local/bin/tabix")
BGZIP=$(which bgzip 2>/dev/null || echo "/usr/local/bin/bgzip")

# Colour output
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
fail() { echo -e "${RED}[FAIL]${NC} $*"; }

echo "============================================================"
echo "  VI_Pipeline — Annotation Resource Downloader"
echo "  Destination: $DEST"
echo "============================================================"
echo ""

mkdir -p "$DEST"

# =============================================================================
# 1. ClinVar
#    Source: NCBI FTP (no authentication required)
#    Files:  clinvar.vcf.gz + .tbi (index provided by NCBI)
# =============================================================================
echo "------------------------------------------------------------"
echo "  1/3  ClinVar"
echo "------------------------------------------------------------"

CLINVAR_VCF="$DEST/clinvar.vcf.gz"
CLINVAR_TBI="$DEST/clinvar.vcf.gz.tbi"

if [ -f "$CLINVAR_VCF" ] && [ -f "$CLINVAR_TBI" ]; then
    ok "ClinVar already present — skipping download"
    ok "  $CLINVAR_VCF"
else
    echo "  Downloading clinvar.vcf.gz (~120 MB) ..."
    wget -q --show-progress -P "$DEST" \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    echo "  Downloading clinvar.vcf.gz.tbi ..."
    wget -q --show-progress -P "$DEST" \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
    ok "ClinVar downloaded: $CLINVAR_VCF"
fi

# =============================================================================
# 2. REVEL
#    Source: Zenodo (doi.org/10.5281/zenodo.7072866 — REVEL v1.3)
#    File:   revel-v1.3_all_chromosomes.zip  (~7 GB)
#    Prep:   unzip → bgzip → tabix index
#
#    The tabix command for REVEL:
#      tabix -f -s 1 -b 2 -e 2 revel_with_transcript_ids  (column layout)
#    Column 1 = chrom (no 'chr' prefix), column 2 = position
# =============================================================================
echo ""
echo "------------------------------------------------------------"
echo "  2/3  REVEL"
echo "------------------------------------------------------------"

REVEL_GZ="$DEST/revel_all_chromosomes.tsv.gz"

if [ -f "$REVEL_GZ" ] && [ -f "${REVEL_GZ}.tbi" ]; then
    ok "REVEL already present — skipping download"
    ok "  $REVEL_GZ"
else
    echo "  Downloading REVEL v1.3 scores (~7 GB zip) ..."
    echo "  This will take several minutes on a slow connection."
    wget -q --show-progress -O "$DEST/revel_v1.3.zip" \
        "https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip"

    echo "  Unzipping ..."
    unzip -o "$DEST/revel_v1.3.zip" -d "$DEST/"

    # The unzipped file is tab-separated with a header; bgzip and index it
    REVEL_RAW=$(ls "$DEST"/revel_with_transcript_ids 2>/dev/null \
             || ls "$DEST"/revel-v1.3_all_chromosomes.tsv 2>/dev/null \
             || ls "$DEST"/revel_all_chromosomes.tsv 2>/dev/null \
             || echo "")

    if [ -z "$REVEL_RAW" ]; then
        fail "Could not find unzipped REVEL file in $DEST — check zip contents"
        echo "  Contents of $DEST:"
        ls "$DEST"
        exit 1
    fi

    # File is CSV: chr,hg19_pos,grch38_pos,ref,alt,...
    # Must: comment header with #, filter rows lacking grch38_pos, convert to TSV
    echo "  Converting CSV to TSV (commenting header, filtering empty grch38 rows) ..."
    awk -F',' 'NR==1{printf "#"; print; next} $3!="" && $3!="." {print}' "$REVEL_RAW" \
      | tr ',' '\t' \
      | "$BGZIP" -c > "$REVEL_GZ"
    rm -f "$REVEL_RAW" "$DEST/revel_v1.3.zip"

    echo "  Indexing with tabix (chr=col1, grch38_pos=col3) ..."
    "$TABIX" -f -s 1 -b 3 -e 3 "$REVEL_GZ"

    ok "REVEL prepared: $REVEL_GZ"
fi

# =============================================================================
# 3. SpliceAI
#    SNV  source: Ensembl FTP (Ensembl MANE v1.4, GRCh38) — no auth required
#    Indel source: Illumina BaseSpace — free account required; instructions
#                  printed below if file is not present.
# =============================================================================
echo ""
echo "------------------------------------------------------------"
echo "  3/3  SpliceAI (SNV + indel)"
echo "------------------------------------------------------------"

SPLICE_SNV="$DEST/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz"
SPLICE_INDEL="$DEST/spliceai_scores.raw.indel.hg38.vcf.gz"

if [ -f "$SPLICE_SNV" ] && [ -f "${SPLICE_SNV}.tbi" ]; then
    ok "SpliceAI SNV already present — skipping download"
    ok "  $SPLICE_SNV"
else
    echo "  Downloading SpliceAI SNV scores from Ensembl FTP (~28 GB) ..."
    echo "  This is a large file and will take a long time."
    wget -q --show-progress -P "$DEST" \
        "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz"
    wget -q --show-progress -P "$DEST" \
        "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz.tbi"
    ok "SpliceAI SNV downloaded: $SPLICE_SNV"
fi

if [ -f "$SPLICE_INDEL" ] && [ -f "${SPLICE_INDEL}.tbi" ]; then
    ok "SpliceAI indel already present — skipping"
    ok "  $SPLICE_INDEL"
else
    warn "SpliceAI indel requires a free Illumina BaseSpace account."
    warn "  1. Create a free account at https://basespace.illumina.com"
    warn "  2. Accept the dataset share: https://basespace.illumina.com/s/otSPW8hnhaZR"
    warn "  3. Download with the BaseSpace CLI:"
    warn "       bs dataset download --id <dataset_id> -o $DEST"
    warn "  4. Tabix-index: tabix -p vcf $SPLICE_INDEL"
    warn "  5. Re-run: python cli.py --setup --force"
    warn "  The pipeline runs in SNV-only splice mode until this file is added."
fi

# =============================================================================
# 4. Update paths.yaml
# =============================================================================
echo ""
echo "------------------------------------------------------------"
echo "  Updating config/paths.yaml"
echo "------------------------------------------------------------"

PATHS_YAML="/home/admin7995/VI_Pipeline/config/paths.yaml"

python3 - <<PYEOF
import yaml
from pathlib import Path

paths_file = Path("$PATHS_YAML")
dest = Path("$DEST")

with open(paths_file) as f:
    config = yaml.safe_load(f)

resources = config.setdefault("resources", {})

updates = {
    "clinvar_vcf":         str(dest / "clinvar.vcf.gz"),
    "revel_file":          str(dest / "revel_all_chromosomes.tsv.gz"),
    "spliceai_snv_file":   str(dest / "spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz"),
    "alphamissense_file":  str(dest / "AlphaMissense_hg38.tsv.gz"),
    # spliceai_indel_file requires Illumina BaseSpace — must be added manually
    # dbnsfp_file is handled by download_dbnsfp.sh separately
}

for key, path in updates.items():
    if Path(path).exists():
        resources.setdefault(key, {})["path"] = path
        resources[key]["status"] = "found"
        print(f"  [OK]  {key} → {path}")
    else:
        print(f"  [SKIP] {key} — file not found at {path}")

with open(paths_file, "w") as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
print(f"  Saved: {paths_file}")
PYEOF

# =============================================================================
# 5. Summary
# =============================================================================
echo ""
echo "============================================================"
echo "  Done. Next steps:"
echo ""
echo "  1. Verify setup:"
echo "       python /home/admin7995/VI_Pipeline/cli.py --setup --force"
echo ""
echo "  2. Delete the cached annotated VCF for any sample you"
echo "     want to re-annotate with ClinVar/REVEL/SpliceAI:"
echo "       rm /path/to/runs/SAMPLE/SAMPLE.vep.vcf"
echo ""
echo "  3. Re-run the pipeline — VEP will pick up all three"
echo "     resources automatically."
echo "============================================================"
