#!/bin/bash

# ==============================================================================
# Single Region Variability Study - Run Script
# ==============================================================================
# Runs the single region HRF parameter variability analysis
# Generates neural signal (raw & downsampled) and BOLD signals with 
# 6 variability levels from 10% to 100%
# ==============================================================================

# Configuration
EXECUTABLE="../bin/single_region_variability"
RESULTS_DIR="../results/single_region_variability"
RANDOM_SEED=12345

# Create results directory
echo "Setting up results directory..."
mkdir -p "${RESULTS_DIR}"

echo "==============================================================================="
echo "SINGLE REGION HRF PARAMETER VARIABILITY STUDY"
echo "==============================================================================="
echo "Executable: ${EXECUTABLE}"
echo "Results directory: ${RESULTS_DIR}"
echo "Random seed: ${RANDOM_SEED}"
echo ""
echo "This simulation will generate:"
echo "  - Neural signal (raw): ~520,000 timepoints at 1ms resolution"
echo "  - Neural signal (downsampled): 490 timepoints at BOLD resolution (1s)"
echo "  - BOLD signals: 6 files with variability levels 10%, 28%, 46%, 64%, 82%, 100%"
echo ""
echo "All BOLD signals use the same underlying neural activity but different"
echo "HRF parameters to isolate the effect of parameter variability."
echo ""

# Check if executable exists
if [ ! -f "${EXECUTABLE}" ]; then
    echo "ERROR: Executable not found: ${EXECUTABLE}"
    echo "Run ./compile_single_region.sh first to compile the program."
    exit 1
fi

echo "==============================================================================="
echo "STARTING SIMULATION"
echo "==============================================================================="

# Run the simulation
"${EXECUTABLE}" ${RANDOM_SEED} "${RESULTS_DIR}"

if [ $? -eq 0 ]; then
    echo ""
    echo "==============================================================================="
    echo "SIMULATION COMPLETED SUCCESSFULLY"
    echo "==============================================================================="
    echo ""
    echo "Generated files in ${RESULTS_DIR}:"
    echo "  Neural signals:"
    echo "    - neural_signal_raw.txt          → Raw neural signal (1ms resolution)"
    echo "    - neural_signal_downsampled.txt  → Downsampled neural signal (BOLD resolution)"
    echo ""
    echo "  BOLD signals with HRF variability:"
    echo "    - bold_variability_10_percent.txt  → 10% parameter variability"
    echo "    - bold_variability_28_percent.txt  → 28% parameter variability"
    echo "    - bold_variability_46_percent.txt  → 46% parameter variability"
    echo "    - bold_variability_64_percent.txt  → 64% parameter variability"
    echo "    - bold_variability_82_percent.txt  → 82% parameter variability"
    echo "    - bold_variability_100_percent.txt → 100% parameter variability"
    echo ""
    echo "All files include headers with parameter information and metadata."
    echo "Each BOLD signal is generated from the SAME neural activity using"
    echo "different HRF parameters to isolate variability effects."
    echo ""
    echo "You can now:"
    echo "  1. Plot neural vs BOLD signals to see the relationship"
    echo "  2. Compare BOLD signals across variability levels"
    echo "  3. Analyze the impact of HRF parameter heterogeneity"
    echo ""
    echo "==============================================================================="
else
    echo ""
    echo "✗ Simulation failed"
    exit 1
fi