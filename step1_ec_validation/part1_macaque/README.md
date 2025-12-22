# Part 1 Macaque - EC Validation Analysis

Validates the iEC framework against macaque tract-tracing ground truth data to replicate **Figure 2**.

## Pipeline Overview

### Supplementary Analysis
**Directory:** `supplementary`
- Signal characteristics analysis (non-Gaussian structure, skewness comparison)
- FLN bidirectionality measurement 
- Algorithm parameter sensitivity (VAR lambda)
- Leave-One-Out analysis for iEC model performance
- Comparison of simulation model (Hopf vs. BEI)

### Stage 1: Run EC Algorithms
**File:** `stage1_run_algorithms.m`
- Runs 9 EC algorithms on 19 macaque subjects (40 regions)
- Processes original and deconvolved BOLD signals (rsHRF)
- Algorithms: rDCM, VAR, GC, FASK, CCD, BOSS, LiNGAM, GRASP, Patel

### Stage 2: Integration
**File:** `stage2_integrate.m`
- Estimates integration weights using training subjects (9 subjects)
- Creates integrated EC (iEC) and VFL subset combination
- Statistical testing with FDR correction

### Stage 3: Main Analysis & Figure Generation
**File:** `stage3_main_analysis.m`
- Tests iEC vs ground truth FLN using test subjects (10 subjects)
- Statistical comparisons, F1 scores, SLN validation
- Generates Figure 2 components and supplementary analyses

## Quick Start

**For Figure 2 replication** (loads pre-computed results):
```matlab
cd step1_ec_validation/part1_macaque
run('stage3_main_analysis.m')
```

**For complete pipeline**:
```matlab
run('stage1_run_algorithms.m')       % EC estimation  
run('stage2_integrate.m')            % Integration
run('stage3_main_analysis.m')        # Validation & figures
```