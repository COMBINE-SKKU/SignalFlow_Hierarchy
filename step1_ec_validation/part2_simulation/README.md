# Part 2 Simulation - iEC Framework Validation with Hopf Model

Validates the iEC framework using **Hopf model simulations** on directed human structural connectivity at higher resolutions than macaque data.

## Purpose

Tests iEC effectiveness on:
- **Higher resolution networks** (100-360 regions vs 40 macaque regions)  
- **Complex human connectivity patterns** with known ground truth directionality
- **Multiple network variants** (100 directed SC matrices) for robust validation

## Approach

### Directed SC Generation
- Uses human group-level structural connectivity from DTI
- **Single hemisphere analysis** (avoids DTI limitations with corpus callosum)  
- Applies `randmio_dir_connected` to create 100 directed SC variants
- Preserves network properties while inducing known directionality

## Pipeline Overview

### Supplementary Analysis
**Directory:** `supplementary`
- Generates directed SC matrices using `randmio_dir_connected`
- Compares topological properties (clustering, degree) 
- Visualizes connectivity matrices and correlations
- Visualizes Hopf model signal characteristics
- Compared two simulation models (Hopf vs. BEI)

### Stage 1: Run Algorithms with Hopf Model
**File:** `stage1_run_algorithms.m`
- **Lambda optimization**: Tests VAR parameters (1e1 to 1e8) across 10 iterations
- **Directed SC generation**: Creates 100 directed variants for validation
- **Hopf simulations**: Generates synthetic BOLD from directed SC matrices
- **EC estimation**: Applies 8 algorithms to synthetic data

### Stage 2: Integration and Testing  
**File:** `stage2_integrate_and_test.m`
- **Integration**: Combines EC results using training set (50 iterations)
- **Validation**: Tests iEC against ground truth using test set (50 iterations)  
- **Statistical analysis**: Compares performance with significance testing
- **Visualization**: Creates correlation boxplots and beta distributions

## Quick Start

**For validation results** (loads pre-computed):
```matlab
cd step1_ec_validation/part2_simulation  
run('stage2_integrate_and_test.m')
```

**For complete pipeline**:
```matlab
run('supplementary_analysis.m')           % Generate directed SC
run('stage1_run_algorithms.m')            % Hopf sims + EC estimation
run('stage2_integrate_and_test.m')        % Integration + validation
```

**Expected Output**: Correlation boxplots comparing iEC vs individual algorithms, beta weight distributions, and statistical validation results across 100 directed SC variants.

## Output

Validates iEC framework on controlled synthetic data with known ground truth, demonstrating effectiveness across multiple network complexities and scales.