#!/bin/bash
# Build all three ablation study configurations
# and copy binaries to separate directories for easy comparison

set -e  # Exit on error

echo "=========================================="
echo "Building All KERO Ablation Configurations"
echo "=========================================="
echo ""

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "Script directory: ${SCRIPT_DIR}"
echo "Project root: ${PROJECT_ROOT}"
echo ""

# Create output directory for binaries in project root
ABLATION_BIN_DIR="${PROJECT_ROOT}/ablation_binaries"
rm -rf ${ABLATION_BIN_DIR}
mkdir -p ${ABLATION_BIN_DIR}

# Build Config 1: Row
echo ">>> Building Config 1: KERO-Row..."
bash "${SCRIPT_DIR}/build_ablation_row.sh"
cp "${PROJECT_ROOT}/bin/Release/kero-tools" "${ABLATION_BIN_DIR}/kero-tools-row"
echo "✓ Config 1 complete: ${ABLATION_BIN_DIR}/kero-tools-row"
echo ""

# Build Config 2: Columnar (no integer array compression)
echo ">>> Building Config 2: KERO-Col..."
bash "${SCRIPT_DIR}/build_ablation_col.sh"
cp "${PROJECT_ROOT}/bin/Release/kero-tools" "${ABLATION_BIN_DIR}/kero-tools-col"
echo "✓ Config 2 complete: ${ABLATION_BIN_DIR}/kero-tools-col"
echo ""

# Build Config 3: Columnar + integer array compression
echo ">>> Building Config 3: KERO-Final..."
bash "${SCRIPT_DIR}/build_ablation_final.sh"
cp "${PROJECT_ROOT}/bin/Release/kero-tools" "${ABLATION_BIN_DIR}/kero-tools-final"
echo "✓ Config 3 complete: ${ABLATION_BIN_DIR}/kero-tools-final"
echo ""

echo "=========================================="
echo "All builds complete!"
echo "=========================================="
echo ""
echo "Binaries available in: ${ABLATION_BIN_DIR}/"
echo "  - kero-tools-row     (Row-oriented, no compression)"
echo "  - kero-tools-col     (Columnar, no integer array compression)"
echo "  - kero-tools-final   (Columnar + integer array compression)"
echo ""
echo "Next steps:"
echo "1. Prepare your dataset"
echo "2. Run experiments with each binary"
echo "3. Compare file sizes and query performance"
echo ""
