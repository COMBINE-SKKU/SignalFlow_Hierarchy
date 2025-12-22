#!/usr/bin/env python3
"""
Single Region Variability Study - Publication Quality Figure
Creates a clean figure showing z-scored BOLD responses with varying HRF parameter variability
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import os
from scipy import stats

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

def z_score(signal):
    """Z-score normalize a signal."""
    return stats.zscore(signal)

def main():
    # Data directory
    results_dir = '../results/single_region_variability'
    
    # Time window parameters (50 points in the middle)
    window_size = 50
    total_points = 490
    start_idx = (total_points - window_size) // 2
    end_idx = start_idx + window_size
    time = np.arange(window_size)
    
    # Load BOLD signals
    variability_levels = [0, 30, 60, 90]
    bold_signals = {}
    bold_signals_zscore = {}
    
    for var_level in variability_levels:
        if var_level == 0:
            filepath = os.path.join(results_dir, 'bold_control.txt')
        else:
            filepath = os.path.join(results_dir, f'bold_variability_{var_level:02d}_percent.txt')
        signal = load_signal(filepath)
        bold_signals[var_level] = signal[start_idx:end_idx]
        # Z-score normalize
        bold_signals_zscore[var_level] = z_score(bold_signals[var_level])
    
    # Create figure with single plot
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    
    # Color scheme - viridis colormap
    cmap = cm.get_cmap('viridis')
    colors = [cmap(i / (len(variability_levels) - 1)) for i in range(len(variability_levels))]
    
    # Plot z-scored BOLD signals with different variability levels
    for i, var_level in enumerate(variability_levels):
        label = f'{var_level}%' if var_level in [30, 90] else None  # Only label first and last
        if var_level == 0:
            ax.plot(time, bold_signals_zscore[var_level], color='black', linewidth=1.5, 
                    label=label, alpha=0.9)
        else:
            ax.plot(time, bold_signals_zscore[var_level], color=colors[i], linewidth=1.5, 
                    label=label, alpha=0.9)
    
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('BOLD signal (z-score)', fontsize=11)
    ax.tick_params(axis='both', labelsize=9)
    
    # Add colorbar for variability gradient
    sm = cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(vmin=10, vmax=100))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02, aspect=20)
    cbar.set_label('HRF variability (%)', fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Adjust x-axis
    ax.set_xlim([0, window_size-1])
    ax.set_xticks([0, 10, 20, 30, 40])
    
    # Add horizontal line at y=0 for reference
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    # Set y-axis limits for better visualization
    ax.set_ylim([-3, 3])
    
    # Also save as PNG for quick viewing
    output_png = os.path.join(results_dir, 'single_region_variability_figure.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"PNG version saved to: {output_png}")
    
    plt.close()

if __name__ == "__main__":
    main()