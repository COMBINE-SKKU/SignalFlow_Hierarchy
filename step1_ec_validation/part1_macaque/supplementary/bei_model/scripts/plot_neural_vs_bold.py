#!/usr/bin/env python3
"""
Neural vs BOLD Signal Visualization
Shows the relationship between underlying neural activity and BOLD response
using a stem plot for neural activity and line plot for BOLD signal
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

# Set matplotlib parameters for publication quality
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_signal(filepath):
    """Load signal from text file, skipping header lines."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines
    data_lines = [line for line in lines if not line.startswith('#')]
    return np.array([float(line.strip()) for line in data_lines])

def main():
    # --- 1. Load Data ---
    results_dir = '../results/single_region_variability'
    
    # Load neural signal (raw, 1ms resolution)
    neural_raw = load_signal(os.path.join(results_dir, 'neural_signal_raw.txt'))
    
    # Load BOLD signal (default/control)
    bold_signal = load_signal(os.path.join(results_dir, 'bold_default.txt'))
    
    print(f"Loaded neural signal: {len(neural_raw)} timepoints")
    print(f"Loaded BOLD signal: {len(bold_signal)} timepoints")
    
    # --- 2. Solve the Time Problem (Downsampling) ---
    # We have 1000 neural points for every 1 BOLD point (1ms vs 1000ms TR)
    tr_duration_ms = 1000
    num_trs = len(bold_signal)
    
    # Take only the neural data corresponding to the BOLD recording period
    # Skip burn-in period neural data
    neural_main_phase = neural_raw[-num_trs * tr_duration_ms:]
    
    # Reshape neural data to (num_trs, 1000) and average along axis 1
    neural_tr_averaged = np.mean(
        neural_main_phase.reshape(num_trs, tr_duration_ms), 
        axis=1
    )
    
    print(f"Neural signal after TR averaging: {len(neural_tr_averaged)} timepoints")
    
    # --- 3. Solve the Unit Problem (Z-Scoring) ---
    # Now both signals have the same number of points. Let's z-score them.
    neural_zscored = stats.zscore(neural_tr_averaged)
    bold_zscored = stats.zscore(bold_signal)
    
    # Select a 50-second window to plot (from the middle)
    total_points = len(bold_signal)
    window_size = 50
    start_idx = (total_points - window_size) // 2
    end_idx = start_idx + window_size
    time = np.arange(window_size)
    
    print(f"Plotting window: {start_idx} to {end_idx} (50 timepoints)")
    
    # --- 4. Create the Plot ---
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot the BOLD signal as a smooth, strong line
    ax.plot(
        time, 
        bold_zscored[start_idx:end_idx], 
        color='black', 
        linewidth=2.5, 
        label='BOLD Signal (z-score)',
        zorder=2
    )
    
    # Plot the NEURAL signal as a "lollipop" stem plot
    # We set the bottom of the stems to 0 for clarity
    markerline, stemlines, baseline = ax.stem(
        time, 
        neural_zscored[start_idx:end_idx], 
        linefmt='red', 
        markerfmt='ro', 
        basefmt=' ',  # No baseline
        bottom=0,
        label='Neural Activity (z-score, avg per TR)'
    )
    plt.setp(markerline, markersize=4)
    plt.setp(stemlines, linewidth=1.5, alpha=0.6)
    
    # --- 5. Final Touches ---
    ax.set_title('Underlying Neural Activity and BOLD Response', fontsize=16, fontweight='bold')
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Standard Deviations (z-score)', fontsize=12)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Legend and formatting
    ax.legend(frameon=False, fontsize=11, loc='upper right')
    ax.set_xlim([0, window_size-1])
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.tick_params(axis='both', labelsize=10)
    
    # Calculate and display correlation
    correlation = np.corrcoef(neural_zscored[start_idx:end_idx], bold_zscored[start_idx:end_idx])[0, 1]
    ax.text(0.02, 0.98, f'Correlation: r = {correlation:.3f}', 
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Save figure
    output_path = os.path.join(results_dir, 'neural_vs_bold_stem.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    
    output_png = os.path.join(results_dir, 'neural_vs_bold_stem.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"PNG version saved to: {output_png}")
    
    plt.close()
    
    # --- 6. Print Summary Statistics ---
    print(f"\n=== ANALYSIS SUMMARY ===")
    print(f"Neural activity (z-score): mean={np.mean(neural_zscored):.3f}, std={np.std(neural_zscored):.3f}")
    print(f"BOLD signal (z-score): mean={np.mean(bold_zscored):.3f}, std={np.std(bold_zscored):.3f}")
    print(f"Overall correlation: r = {np.corrcoef(neural_zscored, bold_zscored)[0, 1]:.3f}")
    print(f"Window correlation: r = {correlation:.3f}")

if __name__ == "__main__":
    main()