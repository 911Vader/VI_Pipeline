#!/usr/bin/env bash
# =============================================================================
# download_dbnsfp.sh
#
# Downloads dbNSFP 4.1a (GRCh38) from Zenodo and prepares it for the
# Ensembl VEP dbNSFP plugin.
#
# Steps:
#   1. Download individual chromosome files from Zenodo record 4323592
#   2. Merge into a single file with one header line
#   3. Sort by chr+pos, bgzip, tabix index
#   4. Update config/paths.yaml
#
# Storage: ~30 GB compressed; ~80 GB temp space needed during merge/sort
# Destination: /mnt/server_share/Reference_Material/annotation/
# =============================================================================

set -euo pipefail

DEST="/mnt/server_share/Reference_Material/annotation"
BGZIP=$(which bgzip 2>/dev/null || echo "/usr/local/bin/bgzip")
TABIX=$(which tabix 2>/dev/null || echo "/usr/local/bin/tabix")
OUT="$DEST/dbNSFP4.1a_grch38.gz"
ZENODO_BASE="https://zenodo.org/api/records/4323592/files"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }

echo "============================================================"
echo "  dbNSFP 4.1a Download and Preparation"
echo "  Destination: $DEST"
echo "  This will take 30-60 minutes on a fast connection."
echo "============================================================"

if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
    ok "dbNSFP already present — skipping"
    ok "  $OUT"
    exit 0
fi

mkdir -p "$DEST/dbnsfp_tmp"

# ---------------------------------------------------------------------------
# 1. Download chromosome files
# ---------------------------------------------------------------------------
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

echo ""
echo "Downloading ${#CHROMS[@]} chromosome files (~30 GB total) ..."

for chr in "${CHROMS[@]}"; do
    FILE="$DEST/dbnsfp_tmp/dbNSFP4.1a_variant.chr${chr}.gz"
    if [ -f "$FILE" ]; then
        warn "chr${chr} already downloaded — skipping"
    else
        echo "  Downloading chr${chr} ..."
        wget -q --show-progress \
            -O "$FILE" \
            "${ZENODO_BASE}/dbNSFP4.1a_variant.chr${chr}.gz/content"
    fi
done

ok "All chromosome files downloaded"

# ---------------------------------------------------------------------------
# 2. Merge: extract header from chr1, then append all data rows sorted
# ---------------------------------------------------------------------------
echo ""
echo "Merging chromosome files ..."

HEADER_FILE="$DEST/dbnsfp_tmp/header.txt"
zcat "$DEST/dbnsfp_tmp/dbNSFP4.1a_variant.chr1.gz" | head -1 > "$HEADER_FILE"

echo "  Sorting and merging data rows (this may take 15-20 minutes) ..."
{
    cat "$HEADER_FILE"
    for chr in "${CHROMS[@]}"; do
        zcat "$DEST/dbnsfp_tmp/dbNSFP4.1a_variant.chr${chr}.gz" | tail -n +2
    done
} | "$BGZIP" -c > "$OUT"

ok "Merge complete: $OUT"

# ---------------------------------------------------------------------------
# 3. Tabix index
# ---------------------------------------------------------------------------
echo ""
echo "Indexing with tabix ..."
"$TABIX" -s 1 -b 2 -e 2 "$OUT"
ok "Index complete: ${OUT}.tbi"

# ---------------------------------------------------------------------------
# 4. Cleanup temp files
# ---------------------------------------------------------------------------
echo ""
echo "Cleaning up temporary files ..."
rm -rf "$DEST/dbnsfp_tmp"
ok "Cleanup done"

# ---------------------------------------------------------------------------
# 5. Update paths.yaml
# ---------------------------------------------------------------------------
python3 - <<PYEOF
import yaml
from pathlib import Path

paths_file = Path("/home/admin7995/VI_Pipeline/config/paths.yaml")
with open(paths_file) as f:
    config = yaml.safe_load(f)

config.setdefault("resources", {}).setdefault("dbnsfp_file", {})
config["resources"]["dbnsfp_file"]["path"] = "$OUT"
config["resources"]["dbnsfp_file"]["status"] = "found"
config["resources"]["dbnsfp_file"]["notes"] = "dbNSFP 4.1a GRCh38"

with open(paths_file, "w") as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
print(f"  Updated: {paths_file}")
PYEOF

echo ""
echo "============================================================"
echo "  dbNSFP preparation complete."
echo "  File: $OUT"
echo "  Next: python /home/admin7995/VI_Pipeline/cli.py --setup --force"
echo "============================================================"
