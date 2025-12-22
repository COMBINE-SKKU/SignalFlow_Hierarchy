#!/bin/bash

# ==============================================================================
# BW Parameter Variability Comparison Study  
# ==============================================================================
# This script uses the dual-BOLD functionality to compare BOLD signals with
# and without BW parameter variability using the SAME underlying neural signal.
#
# Single simulation run produces:
# 1. Variable BOLD: Uses heterogeneous BW parameters (loaded from file) 
# 2. Control BOLD: Uses homogeneous BW parameters (hardcoded defaults)
# 
# This ensures identical neural dynamics for proper variability comparison.
# ==============================================================================

# Configuration
EXECUTABLE="../bin/bei_model"
CONFIG_DIR="../config" 
RESULTS_DIR="../results/variability_study"
BW_DIR="${CONFIG_DIR}/bw_params"

# Simulation parameters
SUBJECT_ID="variability_comparison"
RANDOM_SEED=12345
GLOBAL_COUPLING=1.30
PARAM_FILE="${CONFIG_DIR}/variability_params.txt"

# Variability settings
VARIABILITY_LEVEL=1.0  # 30% variability
REGIONS=40

# Create directories
echo "Setting up directories..."
mkdir -p "${RESULTS_DIR}"
mkdir -p "${BW_DIR}"

# Generate parameter file
echo "Creating parameter file..."
PARAM_LINE="${REGIONS} ${GLOBAL_COUPLING} 0.15 0.0 0.00316228 490000 1000"
echo "${PARAM_LINE}" > "${PARAM_FILE}"
echo "Parameters: ${PARAM_LINE}"

# ==============================================================================
# GENERATE HETEROGENEOUS BW PARAMETERS
# ==============================================================================
echo ""
echo "==============================================================================="
echo "BW PARAMETER VARIABILITY COMPARISON"
echo "==============================================================================="

# Generate heterogeneous BW parameters file
echo "Generating heterogeneous BW parameters (${VARIABILITY_LEVEL}% variability)..."
python3 generate_bw_params.py \
    --regions ${REGIONS} \
    --variability ${VARIABILITY_LEVEL} \
    --seed ${RANDOM_SEED} \
    --output "${BW_DIR}/bw_heterogeneous.txt"

echo "BW parameter file created: ${BW_DIR}/bw_heterogeneous.txt"

# ==============================================================================  
# SINGLE SIMULATION WITH DUAL-BOLD OUTPUT
# ==============================================================================
echo ""
echo "Running single simulation with dual-BOLD output..."
echo "  - Variable BOLD: Uses heterogeneous BW parameters (${VARIABILITY_LEVEL}% variability)"
echo "  - Control BOLD:  Uses homogeneous BW parameters (0% variability, defaults)"
echo ""

${EXECUTABLE} "${PARAM_FILE}" "${SUBJECT_ID}" ${RANDOM_SEED} "${RESULTS_DIR}" "${BW_DIR}/bw_heterogeneous.txt"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Simulation completed successfully!"
    echo ""
else
    echo "✗ Simulation failed"
    exit 1
fi

# ==============================================================================
# RESULTS SUMMARY
# ==============================================================================
echo "==============================================================================="
echo "VARIABILITY COMPARISON STUDY RESULTS"
echo "==============================================================================="
echo ""
echo "Generated files in ${RESULTS_DIR}:"
echo "  - ${SUBJECT_ID}_BOLD_timeseries.txt  → Variable BOLD (${VARIABILITY_LEVEL}% BW variability)"
echo "  - ${SUBJECT_ID}_BOLD_control.txt     → Control BOLD (0% BW variability)"
echo ""
echo "BW parameters used:"
echo "  - Variable model: ${BW_DIR}/bw_heterogeneous.txt (${VARIABILITY_LEVEL}% variability)"
echo "  - Control model:  Built-in defaults (homogeneous)"
echo ""
echo "COMPARISON ANALYSIS:"
echo "Both BOLD signals are generated from the SAME neural activity using:"
echo "  ✓ Same random seed: ${RANDOM_SEED}"
echo "  ✓ Same neural parameters and dynamics"
echo "  ✓ Same excitatory synaptic activity input"
echo "  ✗ Different BW parameters (heterogeneous vs homogeneous)"
echo ""
echo "You can now analyze:"
echo "  1. Signal variance: Compare variability between the two BOLD signals"
echo "  2. Cross-correlation: Measure similarity despite BW differences"
echo "  3. Spectral analysis: Compare frequency content"
echo "  4. Regional differences: Examine how variability affects different regions"
echo ""
echo "This isolates the pure effect of BW parameter variability on BOLD characteristics."

# ==============================================================================
# CLEANUP  
# ==============================================================================
echo ""
echo "Cleaning up temporary files..."
rm -f "${PARAM_FILE}"
echo "Cleanup complete."
echo ""
echo "Variability comparison study finished!"
echo "==============================================================================="