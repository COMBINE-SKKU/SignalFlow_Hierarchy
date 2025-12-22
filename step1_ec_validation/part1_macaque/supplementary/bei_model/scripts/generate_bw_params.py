#!/usr/bin/env python3
"""
BW Parameter File Generator for Variability Study

This script generates Balloon-Windkessel parameter files with different
variability levels to study the impact on BOLD signal characteristics.
"""

import numpy as np
import argparse
import os

def generate_bw_parameters(num_regions=40, variability=0.0, output_file=None, seed=42):
    """
    Generate BW parameters with specified variability.
    
    Parameters:
    -----------
    num_regions : int
        Number of brain regions (default: 40)
    variability : float
        Variability as fraction (0.0 = no variability, 0.15 = 15% variability)
    output_file : str
        Output file path (if None, prints to stdout)
    seed : int
        Random seed for reproducibility
    """
    
    # Set random seed for reproducibility
    np.random.seed(seed)
    
    # Default BW parameters (homogeneous baseline)
    default_rho = 0.34
    default_alpha = 0.32
    default_tau = 0.98
    default_kappa = 1.0 / 0.65
    default_gamma = 1.0 / 0.41
    default_V0 = 0.02
    
    # Generate parameters for each region
    params = []
    for i in range(num_regions):
        if variability == 0.0:
            # No variability - all regions identical
            rho = default_rho
            alpha = default_alpha
            tau = default_tau
            kappa = default_kappa
            gamma = default_gamma
            V0 = default_V0
        else:
            # Add variability: param = baseline * (1 + variability * random_factor)
            # where random_factor is in [-1, 1]
            rho = default_rho * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            alpha = default_alpha * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            tau = default_tau * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            kappa = default_kappa * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            gamma = default_gamma * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            V0 = default_V0 * (1.0 + variability * (np.random.random() - 0.5) * 2.0)
            
            # Ensure all parameters remain positive
            rho = max(0.1, rho)
            alpha = max(0.1, alpha)
            tau = max(0.1, tau)
            kappa = max(0.1, kappa)
            gamma = max(0.1, gamma)
            V0 = max(0.001, V0)
        
        params.append([i, rho, alpha, tau, kappa, gamma, V0])
    
    # Format output
    header = f"{num_regions}\n"
    lines = [header]
    for param_set in params:
        line = f"{param_set[0]:d} {param_set[1]:.6f} {param_set[2]:.6f} {param_set[3]:.6f} {param_set[4]:.6f} {param_set[5]:.6f} {param_set[6]:.6f}\n"
        lines.append(line)
    
    # Output to file or stdout
    if output_file:
        with open(output_file, 'w') as f:
            f.writelines(lines)
        print(f"BW parameters written to {output_file}")
        print(f"Variability: {variability*100:.1f}%, Regions: {num_regions}, Seed: {seed}")
    else:
        print(''.join(lines), end='')
    
    return params

def main():
    parser = argparse.ArgumentParser(description='Generate BW parameter files for variability study')
    parser.add_argument('--regions', type=int, default=40, help='Number of regions (default: 40)')
    parser.add_argument('--variability', type=float, default=0.15, help='Variability fraction (default: 0.15)')
    parser.add_argument('--output', type=str, help='Output file path')
    parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')
    
    args = parser.parse_args()
    
    generate_bw_parameters(
        num_regions=args.regions,
        variability=args.variability,
        output_file=args.output,
        seed=args.seed
    )

if __name__ == "__main__":
    main()