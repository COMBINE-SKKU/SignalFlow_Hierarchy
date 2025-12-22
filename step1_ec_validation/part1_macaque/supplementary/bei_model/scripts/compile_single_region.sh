#!/bin/bash

# ==============================================================================
# Single Region Variability Study - Compilation Script
# ==============================================================================
# Compiles the single region HRF parameter variability study program
# ==============================================================================

# Configuration
SOURCE_DIR="../src"
BIN_DIR="../bin"
SOURCE_FILE="single_region_variability.c"
EXECUTABLE="single_region_variability"

# Create bin directory if it doesn't exist
mkdir -p "${BIN_DIR}"

echo "==============================================================================="
echo "COMPILING SINGLE REGION VARIABILITY STUDY"
echo "==============================================================================="
echo "Source file: ${SOURCE_DIR}/${SOURCE_FILE}"
echo "Output: ${BIN_DIR}/${EXECUTABLE}"
echo ""

# Compilation with optimization and SIMD support
gcc -O3 -march=native -msse -msse2 -msse3 -msse4.1 -msse4.2 -mavx \
    -fopenmp -ffast-math -funroll-loops \
    -o "${BIN_DIR}/${EXECUTABLE}" \
    "${SOURCE_DIR}/${SOURCE_FILE}" \
    -lm

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful!"
    echo "Executable created: ${BIN_DIR}/${EXECUTABLE}"
    echo ""
    echo "Usage: ./${EXECUTABLE} <random_seed> <results_dir>"
    echo ""
else
    echo "✗ Compilation failed!"
    exit 1
fi

echo "==============================================================================="