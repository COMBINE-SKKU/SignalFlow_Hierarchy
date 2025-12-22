#!/usr/bin/env python3
"""
Quick Analysis of BW Parameter Variability Effect

This script provides a basic analysis of the generated BOLD signals
to demonstrate the effect of BW parameter variability.
"""

import numpy as np
import argparse
import os

def load_bold_data(filepath):
    """Load BOLD timeseries data from file."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        header = lines[0].strip().split()
        n_regions = int(header[0])
        n_timepoints = int(header[1])
        
        # Parse data
        bold_data = np.zeros((n_regions, n_timepoints))
        for i in range(1, n_regions + 1):
            values = list(map(float, lines[i].strip().split()))
            bold_data[i-1, :] = values
        
        return bold_data, n_regions, n_timepoints
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, 0, 0

def analyze_variability(variable_bold, control_bold, regions, timepoints):
    """Analyze the effect of BW parameter variability."""
    
    print("\n" + "="*80)
    print("BW PARAMETER VARIABILITY ANALYSIS")
    print("="*80)
    
    # 1. Signal variance across regions
    var_variance_across_regions = np.var(variable_bold, axis=0).mean()  # Mean variance across time
    ctrl_variance_across_regions = np.var(control_bold, axis=0).mean()
    
    # 2. Signal variance across time
    var_variance_across_time = np.var(variable_bold, axis=1).mean()  # Mean variance across regions
    ctrl_variance_across_time = np.var(control_bold, axis=1).mean()
    
    # 3. Cross-correlation between models
    correlations = []
    for region in range(regions):
        corr = np.corrcoef(variable_bold[region, :], control_bold[region, :])[0, 1]
        correlations.append(corr)
    mean_correlation = np.mean(correlations)
    min_correlation = np.min(correlations)
    max_correlation = np.max(correlations)
    
    # 4. Signal amplitude differences
    var_mean_amplitude = np.mean(np.abs(variable_bold))
    ctrl_mean_amplitude = np.mean(np.abs(control_bold))
    
    # 5. Signal range (max - min) per region
    var_signal_ranges = np.max(variable_bold, axis=1) - np.min(variable_bold, axis=1)
    ctrl_signal_ranges = np.max(control_bold, axis=1) - np.min(control_bold, axis=1)
    var_mean_range = np.mean(var_signal_ranges)
    ctrl_mean_range = np.mean(ctrl_signal_ranges)
    
    print(f"Dataset: {regions} regions × {timepoints} timepoints")
    print()
    print("VARIANCE ANALYSIS:")
    print(f"  Variable BOLD - Spatial variance:  {var_variance_across_regions:.6f}")
    print(f"  Control BOLD  - Spatial variance:  {ctrl_variance_across_regions:.6f}")
    print(f"  Spatial variance ratio:             {var_variance_across_regions/ctrl_variance_across_regions:.3f}x")
    print()
    print(f"  Variable BOLD - Temporal variance: {var_variance_across_time:.6f}")
    print(f"  Control BOLD  - Temporal variance: {ctrl_variance_across_time:.6f}")
    print(f"  Temporal variance ratio:            {var_variance_across_time/ctrl_variance_across_time:.3f}x")
    print()
    print("CORRELATION ANALYSIS:")
    print(f"  Mean cross-correlation:   {mean_correlation:.4f}")
    print(f"  Min cross-correlation:    {min_correlation:.4f}")
    print(f"  Max cross-correlation:    {max_correlation:.4f}")
    print(f"  Correlation range:        {max_correlation - min_correlation:.4f}")
    print()
    print("AMPLITUDE ANALYSIS:")
    print(f"  Variable BOLD mean amplitude: {var_mean_amplitude:.4f}")
    print(f"  Control BOLD mean amplitude:  {ctrl_mean_amplitude:.4f}")
    print(f"  Amplitude difference:         {abs(var_mean_amplitude - ctrl_mean_amplitude):.4f}")
    print()
    print("SIGNAL RANGE ANALYSIS:")
    print(f"  Variable BOLD mean range:  {var_mean_range:.4f}")
    print(f"  Control BOLD mean range:   {ctrl_mean_range:.4f}")
    print(f"  Range difference:          {abs(var_mean_range - ctrl_mean_range):.4f}")
    print(f"  Range ratio:               {var_mean_range/ctrl_mean_range:.3f}x")
    print()
    
    # Regional variability
    region_var_diff = np.abs(var_signal_ranges - ctrl_signal_ranges)
    most_affected_region = np.argmax(region_var_diff)
    least_affected_region = np.argmin(region_var_diff)
    
    print("REGIONAL IMPACT:")
    print(f"  Most affected region:  #{most_affected_region} (range diff: {region_var_diff[most_affected_region]:.4f})")
    print(f"  Least affected region: #{least_affected_region} (range diff: {region_var_diff[least_affected_region]:.4f})")
    print()
    
    print("INTERPRETATION:")
    if var_variance_across_regions > ctrl_variance_across_regions * 1.1:
        print("  ✓ BW parameter variability increases spatial signal diversity")
    else:
        print("  - BW parameter variability has minimal effect on spatial diversity")
        
    if mean_correlation > 0.9:
        print(f"  ✓ Strong correlation ({mean_correlation:.3f}) confirms same neural input")
    else:
        print(f"  ! Moderate correlation ({mean_correlation:.3f}) suggests different dynamics")
        
    if var_mean_range > ctrl_mean_range * 1.05:
        print("  ✓ Variable BW parameters increase BOLD signal dynamic range")
    else:
        print("  - Variable BW parameters have minimal effect on signal range")
    
    print("="*80)

def main():
    parser = argparse.ArgumentParser(description='Analyze BW parameter variability effect on BOLD signals')
    parser.add_argument('--variable', type=str, 
                       default='../results/variability_study/variability_comparison_BOLD_timeseries.txt',
                       help='Variable BOLD timeseries file')
    parser.add_argument('--control', type=str,
                       default='../results/variability_study/variability_comparison_BOLD_control.txt',
                       help='Control BOLD timeseries file')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not os.path.exists(args.variable):
        print(f"Error: Variable BOLD file not found: {args.variable}")
        print("Run ./compare_variability.sh first to generate the data.")
        return
        
    if not os.path.exists(args.control):
        print(f"Error: Control BOLD file not found: {args.control}")
        print("Run ./compare_variability.sh first to generate the data.")
        return
    
    # Load data
    print("Loading BOLD data...")
    variable_bold, v_regions, v_timepoints = load_bold_data(args.variable)
    control_bold, c_regions, c_timepoints = load_bold_data(args.control)
    
    if variable_bold is None or control_bold is None:
        print("Error loading BOLD data files.")
        return
    
    if v_regions != c_regions or v_timepoints != c_timepoints:
        print(f"Error: Data dimensions mismatch")
        print(f"Variable: {v_regions}×{v_timepoints}, Control: {c_regions}×{c_timepoints}")
        return
    
    # Analyze
    analyze_variability(variable_bold, control_bold, v_regions, v_timepoints)

if __name__ == "__main__":
    main()