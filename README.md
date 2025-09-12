## References

Oh YH, Ann YJ, Lee JJ, Ito T, Froudist-Walsh S, Paquola C, Milham M, Spreng RN, Margulies D, Bernhardt B, Woo CW, Hong SK. [**In vivo cartography of state-dependent signal flow hierarchy in the human cerebral cortex **], Arxiv 2024.

----

## Background (Abstract)

Understanding the principle of information flow across distributed brain networks is of paramount importance in neuroscience. Here, we introduce a novel neuroimaging framework leveraging integrated effective connectivity (iEC) and unconstrained signal flow mapping for data-driven discovery of the human cerebral functional hierarchy, one of the principal brain scaffolds for macroscale neural dynamics. Simulation and empirical validation demonstrated the high fidelity of iEC in recovering connectome directionality and its potential relationship with histologically defined feedforward and feedback pathways. Notably, the iEC-based hierarchy revealed a functional axis composed of sequentially ordered sensorimotor-association-paralimbic areas, a pattern supported by the Structural Model of laminar connectivity. This hierarchy was further demonstrated to flexibly reorganize according to brain state, elevating exteroceptive regions during external focus and interoceptive regions during an internally oriented state. Our study highlights the unique role of macroscale directed functional connections in uncovering a neurobiologically grounded, state-dependent signal flow hierarchy.

## System Requirements

**Operating System:** Linux (tested and validated)
- Windows and macOS may work but are not officially supported

**Software Dependencies:**
- **MATLAB R2020a or later** (core computational environment)
- **Java Runtime Environment** (required for Tetrad-based causal discovery algorithms)
- **[cifti-matlab toolbox](https://github.com/Washington-University/cifti-matlab)** (for HCP data processing)
- **XCode** (macOS users only, for rDCM C compilation)

**Hardware Recommendations:**
- **RAM:** ≥16GB (32GB+ recommended for large-scale analyses)
- **CPU:** Multi-core processor (parallel processing utilized)
- **Storage:** ≥50GB free space for datasets and results

All specialized toolboxes and utilities are included in the `utils/` directory.

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

## Quick Start

1. **Clone repository and add paths:**
   ```matlab
   addpath(genpath('/path/to/Signalflow_Hierarchy'));
   addpath('/path/to/cifti-matlab');  % Add cifti-matlab toolbox
   ```

2. **Run validation analyses:**
   ```matlab
   % iEC framework validation (three approaches)
   cd step1_ec_validation/part1_macaque; run('stage3_main_analysis.m');
   cd ../part2_simulation; run('stage2_integrate_and_test.m');
   cd ../part3_empirical_human; run('stage3_simulation_validation.m');
   
   % Cortical hierarchy analysis
   cd ../../step2_ec_hierarchy; run('stage1_network_profile_main.m'); run('stage2_cortical_hierarchy_main.m');
   
   % State-dependent hierarchy analysis
   cd ../step3_state_hierarchies; run('stage1_main.m');
   ```

3. **Results are saved in each directory's `results/` folder**
