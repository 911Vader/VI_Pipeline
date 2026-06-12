#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Thin shell wrapper for VI_Pipeline
#
# Usage:
#   ./run_pipeline.sh --setup
#   ./run_pipeline.sh --input sample.vcf.gz --gene-list config/gene_list.txt \
#                     --output results/ --sample-id S001 --case-id C001
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Require python3
if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 not found on PATH."
    echo "       Install Python 3.10+ and ensure it is on your PATH."
    exit 1
fi

PYTHON_VERSION=$(python3 -c "import sys; print(sys.version_info.major * 100 + sys.version_info.minor)")
if [ "$PYTHON_VERSION" -lt 310 ]; then
    echo "ERROR: Python 3.10 or later is required."
    echo "       Found: $(python3 --version)"
    exit 1
fi

exec python3 cli.py "$@"
