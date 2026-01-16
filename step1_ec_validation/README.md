# Step 1: EC Validation

This directory validates the integrated effective connectivity (iEC) framework using three complementary approaches. Each part demonstrates the framework's effectiveness on different data types and validation strategies.

## Directory Structure

### Part 1 Macaque (`part1_macaque/`)
**Validation against ground truth connectivity**
- Uses macaque fMRI data (19 subjects, 40 regions) with tract-tracing ground truth (FLN)
- Validates iEC against known anatomical connectivity patterns
- **Output**: Figure 2 - iEC validation with macaque ground truth

### Part 2 Simulation (`part2_simulation/`)
**Validation with controlled synthetic data**
- Uses Hopf model simulations with directed human structural connectivity
- Tests framework on higher resolution networks (100-360 regions) with known ground truth
- Generates 100 random directed SC variants for robust validation
- **Output**: Simulation validation results and performance metrics

### Part 3 Empirical Human (`part3_empirical_human/`)
**Validation with human neuroimaging data**
- Uses HCP dataset (440 subjects) with multiple network resolutions
- Validates iEC performance against individual algorithms using goodness-of-fit metrics
- **Output**: Human empirical validation results and comparative analysis

## Quick Start

```matlab
cd step1_ec_validation/part1_macaque
run('stage3_main_analysis.m')    % Figure 2 - macaque validation
```

**Expected Output**: Figure 2 panels showing iEC validation against macaque tract-tracing ground truth, including correlation analysis and F1 score comparisons.

## Analysis Flow

1. **Ground Truth Validation** (Part 1): Tests against known anatomical connectivity
2. **Synthetic Data Validation** (Part 2): Tests on controlled simulations with known directionality
3. **Empirical Data Validation** (Part 3): Tests on real human data with performance metrics

Each part contains its own README.md with detailed pipeline descriptions and usage instructions.