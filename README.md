## References

Oh YH, Ann YJ, Lee JJ, Ito T, Froudist-Walsh S, Paquola C, Milham M, Spreng RN, Margulies D, Bernhardt B, Woo CW, Hong SK. [**In vivo cartography of state-dependent signal flow hierarchy in the human cerebral cortex **], Arxiv 2024.

----

## Background (Abstract)

Understanding the principle of information flow across distributed brain networks is of paramount importance in neuroscience. Here, we introduce a novel neuroimaging framework leveraging integrated effective connectivity (iEC) and unconstrained signal flow mapping for data-driven discovery of the human cerebral functional hierarchy, one of the principal brain scaffolds for macroscale neural dynamics. Simulation and empirical validation demonstrated the high fidelity of iEC in recovering connectome directionality and its potential relationship with histologically defined feedforward and feedback pathways. Notably, the iEC-based hierarchy revealed a functional axis composed of sequentially ordered sensorimotor-association-paralimbic areas, a pattern supported by the Structural Model of laminar connectivity. This hierarchy was further demonstrated to flexibly reorganize according to brain state, elevating exteroceptive regions during external focus and interoceptive regions during an internally oriented state. Our study highlights the unique role of macroscale directed functional connections in uncovering a neurobiologically grounded, state-dependent signal flow hierarchy.

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
   git clone https://github.com/cocoanlab/Signalflow_Hierarchy.git
   ```

2. Download and install [cifti-matlab](https://github.com/Washington-University/cifti-matlab)

3. Add paths in MATLAB:
   ```matlab
   addpath(genpath('/path/to/Signalflow_Hierarchy'));
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
