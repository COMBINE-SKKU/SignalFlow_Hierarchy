# Step 2: EC Hierarchy Analysis

Analyzes hierarchical organization and signal flow patterns of the iEC matrix to generate Figures 4 and 5.

## Pipeline Overview

### Stage 1: Network Profiling 
**File:** `stage1_network_profile_main.m`

Network structure analysis for **Figure 4**:
- Matrix visualization and edge distribution analysis
- Network degree computations and connection ratios
- Modular patterns across 7 functional networks (Yeo)
- Signal flow between 22 cortical modules
- Cross-module analysis (unimodal vs heteromodal)

### Stage 2: Cortical Hierarchy
**File:** `stage2_cortical_hierarchy_main.m` 

Hierarchy validation for **Figure 5**:
- Derives hierarchy levels from iEC matrix
- PC1 gradient correlation analysis with divergent regions
- Cytoarchitectonic validation across 6 cortical types
- 27-module signal flow patterns and joy plot visualization

## Quick Start

```matlab
cd step2_ec_hierarchy
run('stage1_network_profile_main.m')    % Figure 4
run('stage2_cortical_hierarchy_main.m') # Figure 5  
```

## Key Analyses

**Network Structure**: Heavy-tail distributions, modular connectivity, signal flow dynamics
**Hierarchy Validation**: Functional gradients, cytoarchitectonic correspondence, multi-scale flow

## Output

- **Figure 4**: Matrix visualization, edge distributions, network boxplots, signal flow
- **Figure 5**: Hierarchy maps, PC1 correlation, cytoarchitecture validation, joy plots
- **Supplementary**: Surface plots, statistical analysis results

Demonstrates that iEC reveals neurobiologically meaningful cortical hierarchy aligned with anatomical and functional organization principles.