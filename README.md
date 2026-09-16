# State-dependent signal flow hierarchy in the human cerebral cortex

Code and data for:

Oh Y, Ann Y, Lee J-J, Ito T, Froudist-Walsh S, Paquola C, Milham M, Spreng RN, Margulies D, Bernhardt B, Woo C-W, Hong S-J. **State-dependent signal flow hierarchy in the human cerebral cortex.** *Nature Neuroscience* (2026). https://doi.org/10.1038/s41593-026-02389-8

Published online 1 September 2026. Article page: https://www.nature.com/articles/s41593-026-02389-8

If you use this code, please cite the paper:

```bibtex
@article{Oh2026SignalFlow,
  title   = {State-dependent signal flow hierarchy in the human cerebral cortex},
  author  = {Oh, Younghyun and Ann, Yejin and Lee, Jae-Joong and Ito, Takuya and Froudist-Walsh, Sean and Paquola, Casey and Milham, Michael and Spreng, R. Nathan and Margulies, Daniel and Bernhardt, Boris and Woo, Choong-Wan and Hong, Seok-Jun},
  journal = {Nature Neuroscience},
  year    = {2026},
  doi     = {10.1038/s41593-026-02389-8},
  url     = {https://www.nature.com/articles/s41593-026-02389-8}
}
```

----

## Abstract

Understanding how information flows across distributed brain networks is central to linking brain structure, dynamics and function. Here we present a neuroimaging framework that combines integrated effective connectivity (iEC) and unconstrained signal flow mapping for data-driven identification of human cerebral functional hierarchies. Simulations and empirical validation show that iEC recovers connectome directionality and aligns with histologically defined feedforward and feedback pathways. The iEC-derived hierarchy exhibits a monotonically increasing level along the axis where the sensorimotor, association and paralimbic areas are sequentially ordered, consistent with predictions from the structural model of laminar connectivity. This hierarchy is not fixed but flexibly reorganizes across brain states; it becomes flatter during externally oriented processing and steeper during internally focused conditions, reflecting increased engagement of interoceptive regions. Our study indicates that macroscale directed functional connectivity can reveal biologically grounded, state-dependent principles of signal flow in the human brain.

## System Requirements

**Operating System:**
- Linux (tested on Ubuntu 20.04, 22.04)
- Windows and macOS may work but are not officially supported

**Software Dependencies:**
- **MATLAB R2020a or later** (tested on R2023a, R2024a)
- **Java Runtime Environment 8+** (tested on OpenJDK 11; required for Tetrad-based causal discovery algorithms)
- **[cifti-matlab toolbox](https://github.com/Washington-University/cifti-matlab)** v2.1.0 (for HCP data processing)
- **XCode** (macOS users only, for rDCM C compilation)

**Hardware Requirements:**
- **RAM:** ≥16GB (32GB+ recommended for MMP360 analyses)
- **CPU:** Multi-core processor (parallel processing utilized)
- **Storage:** ≥50GB free space for datasets and results
- No specialized hardware required

All specialized toolboxes and utilities are included in the `utils/` directory.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/COMBINE-SKKU/SignalFlow_Hierarchy.git
   ```

2. Download and install [cifti-matlab](https://github.com/Washington-University/cifti-matlab)

3. Add paths in MATLAB:
   ```matlab
   addpath(genpath('/path/to/SignalFlow_Hierarchy'));
   addpath('/path/to/cifti-matlab');
   ```

4. Verify Java is accessible from MATLAB:
   ```matlab
   version -java  % Should return Java version
   ```

**Typical install time:** <5 minutes (excluding MATLAB installation)

## Project Structure

This repository contains the complete computational pipeline for iEC framework validation and cortical hierarchy analysis. The analyses span three validation approaches plus hierarchical mapping:

- **`step1_ec_validation/`** - Three-part validation of iEC framework effectiveness
- **`step2_ec_hierarchy/`** - Network profiling and cortical hierarchy analysis  
- **`step3_state_hierarchies/`** - State-dependent hierarchy reorganization analysis
- **`data/`** - Neuroimaging datasets and atlases
- **`utils/`** - Core algorithms, signal flow mapping, and visualization tools
- **`ec_algorithms/`** - Effective connectivity estimation methods

Each directory contains comprehensive README files with specific usage instructions and analysis details.

## Core Framework

**iEC Integration:** Bayesian optimization combines multiple effective connectivity algorithms (rDCM, VAR, FASK, CCD, BOSS, LiNGAM, GRASP, Patel) with algorithm-specific weights optimized via Hopf model simulations.

**Signal Flow Mapping:** Unconstrained, data-driven discovery of directed information flow patterns across brain networks without requiring prior structural connectivity constraints.

## Demo

The demo replicates all main figures using pre-computed results.

**Run the demo:**
```matlab
cd step1_ec_validation/part1_macaque
run('stage3_main_analysis.m')           % Figure 2

cd ../../step2_ec_hierarchy
run('stage1_network_profile_main.m')    % Figure 4
run('stage2_cortical_hierarchy_main.m') % Figure 5

cd ../step3_state_hierarchies
run('stage1_main.m')                    % Figure 6
```

**Expected output:**
- Figure 2: iEC validation against macaque tract-tracing ground truth (correlation plots, F1 scores)
- Figure 4: Network structure analysis (matrix visualization, edge distributions, signal flow)
- Figure 5: Cortical hierarchy maps with cytoarchitectonic validation
- Figure 6: State-dependent hierarchy reorganization across brain states

**Expected run time:** ~10 minutes total on a standard desktop (Intel i7, 16GB RAM)

Results are saved in each directory's `results/` folder.

## Instructions for Use

### Replicating Paper Results

Run the demo commands above to generate all main figures from pre-computed iEC matrices.

### Running Complete Pipeline

To recompute EC estimates from raw data (requires significant computation time):

```matlab
% Step 1: EC validation
cd step1_ec_validation/part1_macaque
run('stage1_run_algorithms.m')          % ~2-4 hours
run('stage2_integrate.m')
run('stage3_main_analysis.m')

% Step 2: Hierarchy analysis (uses pre-computed iEC)
cd ../../step2_ec_hierarchy
run('stage1_network_profile_main.m')
run('stage2_cortical_hierarchy_main.m')
```

### Using Your Own Data

To apply the iEC framework to custom fMRI data:

1. **Prepare time series:** Extract ROI time series as a `[time x regions]` matrix
2. **Run EC algorithms:**
   ```matlab
   addpath(genpath('ec_algorithms'));
   ec_results = run_ec_algorithms(timeseries, TR);
   ```
3. **Integrate results:**
   ```matlab
   addpath(genpath('utils'));
   iEC = integrate_ec(ec_results, weights);
   ```
4. **Compute hierarchy:**
   ```matlab
   hierarchy = computeHierarchyLevels(iEC);
   ```

See individual README files in each step directory for detailed parameter options.
