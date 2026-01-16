# Part 3 Empirical Human - iEC Framework Validation with HCP Data

Validates the iEC framework using **human empirical fMRI data** from the Human Connectome Project, demonstrating effectiveness on real-world neuroimaging data.

## Dataset

- **440 HCP subjects** split into training (220) and testing (220) sets
- **Two network resolutions**: Schaefer100 (100 regions) and MMP360 (360 regions)
- **8 EC algorithms**: rDCM, VAR, FASK, CCD, BOSS, LiNGAM, GRASP, Patel

## Pipeline Overview

### Supplementary Analysis
**Directory:** `supplementary`
- **VAR parameter optimization**: Tests lambda values 10-1000 on HCP subset
- **Statistical validation**: t-tests with FDR correction on connectivity matrices
- **Performance evaluation**: Goodness-of-fit using Hopf simulations
- **Deconvolution test**: Investigate the impact of the HRF on EC estimation

### Stage 1: Data Preprocessing
**File:** `stage1_data_preprocessing.m`
- **EC preprocessing**: Group-averages with algorithm-specific thresholding
- **Bayesian algorithm filtering**: Removes edges detected in <10% subjects (excludes rDCM/VAR)
- **Empirical data processing**: Computes FC, FCD, time-delay matrices for train/test sets

### Stage 2: Integration and Weight Optimization  
**File:** `stage2_integration.m`
- **Weight optimization**: Uses Hopf simulations to find optimal integration coefficients
- **Algorithm combinations**: Full (8 algorithms) and reduced (VAR+FASK) variants
- **iEC construction**: Builds integrated connectivity matrices using optimal weights
- **Statistical testing**: Significance analysis with FDR correction

### Stage 3: Simulation Validation
**File:** `stage3_simulation_validation.m`
- **Performance comparison**: Tests all algorithms + iEC using goodness-of-fit metrics  
- **Comprehensive visualization**: GOF boxplots, FC/FCD comparisons, time-delay analysis
- **Stochastic simulation handling**: Option to use pre-stored results for exact replication

## Quick Start

**For figure replication** (uses pre-computed results):
```matlab
cd step1_ec_validation/part3_empirical_human
% Set use_stored_results = true in stage3_simulation_validation.m
run('stage3_simulation_validation.m')
```

**For complete pipeline**:
```matlab
run('supplementary_analysis.m')          % Parameter optimization
run('stage1_data_preprocessing.m')       % Data preprocessing
run('stage2_integration.m')              % Weight optimization + iEC
run('stage3_simulation_validation.m')    % Validation + visualization
```

**Expected Output**: Goodness-of-fit boxplots, FC/FCD comparison matrices, and time-delay analysis showing iEC outperforms individual EC algorithms.

## Important Notes

- **Hopf simulations are stochastic**: New runs produce different results
- **For exact replication**: Use pre-stored results in `[parcellation]_gof_results.mat`
- **Performance metric**: Goodness-of-fit = FC correlation - FCD KS statistic

## Output

Demonstrates that **iEC significantly outperforms individual EC algorithms** in capturing both static and dynamic properties of human brain networks.