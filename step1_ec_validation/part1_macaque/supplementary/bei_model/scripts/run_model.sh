#!/bin/bash

# --- DUAL-BOLD BEI Model Runner ---
# This script runs the BEI model which generates TWO BOLD timeseries:
# 1. Variable BOLD: Uses heterogeneous BW parameters (_BOLD_timeseries.txt)
# 2. Control BOLD: Uses homogeneous BW parameters (_BOLD_control.txt)

# --- Configuration ---
# These paths are relative to THIS script's location (the 'scripts' dir)
EXECUTABLE="../bin/bei_model"
CONFIG_DIR="../config"
RESULTS_DIR="../results"
ITERATIONS=50

# Parameters for the param file
GLOBAL_COUPLING=1.00  # Optimized value from previous analysis
PARAM_FILE="${CONFIG_DIR}/temp_params.txt"

# BW parameter settings (optional)
BW_FILE=""  # Use "" for generated parameters
BW_VARIABILITY=0.5  # 15% variability for heterogeneous parameters

# Simulation parameters (regions, G, J_NMDA, J_i, sigma, timesteps, TR)
# Format: %d %f %f %f %f %d %d
# Updated for 490 BOLD timepoints: 490000 timesteps / 1000 TR = 490 points
PARAM_LINE="40 $GLOBAL_COUPLING 0.15 0.0 0.001 490000 1000"

# --- Setup ---
echo "Ensuring directories exist..."
# Use relative paths to create directories one level up
mkdir -p $CONFIG_DIR
mkdir -p $RESULTS_DIR

echo "Creating temporary parameter file at: $PARAM_FILE"
echo "$PARAM_LINE" > $PARAM_FILE

# --- Main Loop ---
echo "Starting $ITERATIONS dual-BOLD simulation iterations..."
echo "Each iteration will generate:"
echo "  - \$SUBJECT_ID_BOLD_timeseries.txt (variable BW parameters)"
echo "  - \$SUBJECT_ID_BOLD_control.txt (homogeneous BW parameters)"

for i in $(seq 1 $ITERATIONS)
do
    # Use the iteration number to create a unique subject ID and random seed
    SUBJECT_ID="iter_${i}"
    RANDOM_SEED=$i

    echo "--- Running Iteration $i/$ITERATIONS ---"
    echo "  Subject ID: $SUBJECT_ID"
    echo "  Random Seed: $RANDOM_SEED"
    echo "  BW Variability: $BW_VARIABILITY"
    
    # Run the dual-BOLD model with BW parameters
    # Usage: bei_model <param_file> <subject_id> <random_seed> <results_dir> [bw_file] [bw_variability]
    $EXECUTABLE $PARAM_FILE $SUBJECT_ID $RANDOM_SEED $RESULTS_DIR "$BW_FILE" $BW_VARIABILITY
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Generated: ${SUBJECT_ID}_BOLD_timeseries.txt"
        echo "  ✓ Generated: ${SUBJECT_ID}_BOLD_control.txt"
    else
        echo "  ✗ Error in iteration $i"
    fi
    
    echo "Iteration $i complete."
done

# --- Cleanup ---
echo "Cleaning up temporary parameter file..."
rm $PARAM_FILE

echo "Simulation finished."
echo "$ITERATIONS result files are saved in $RESULTS_DIR/"