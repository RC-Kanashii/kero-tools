#!/bin/bash
# Build script for KERO-Row (Ablation Config 1)
# Row-oriented storage, no compression

set -e  # Exit on error

echo "======================================"
echo "Building KERO-Row (Config 1)"
echo "Row-oriented + No Compression"
echo "======================================"

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "Project root: ${PROJECT_ROOT}"
cd "${PROJECT_ROOT}"

# Create build directory in project root
BUILD_DIR="${PROJECT_ROOT}/build_ablation_row"
rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Configure with CMake
cmake -DCMAKE_BUILD_TYPE=Release \
      -DKERO_STORAGE_MODE=ROW \
      "${PROJECT_ROOT}"

# Build
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo ""
echo "======================================"
echo "Build complete!"
echo "Binary location: ${PROJECT_ROOT}/bin/Release/kero-tools"
echo "======================================"
echo ""
echo "To use this version:"
echo "  ${PROJECT_ROOT}/bin/Release/kero-tools --help"
echo ""
