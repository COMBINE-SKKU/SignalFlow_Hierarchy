# Step 3: State Hierarchies Analysis

Analyzes state-dependent changes in cortical hierarchy and signal flow patterns across different brain states to generate Figure 6.

## Pipeline Overview

### Stage 1: State-Dependent Analysis
**File:** `stage1_main.m`

Multi-state cortical hierarchy analysis for **Figure 6**:
- State-specific hierarchy computation across resting, movie-watching, and pain states
- Cortical zone comparison with statistical validation (zones 1-4)
- 27-module signal flow analysis for each brain state
- Delta analysis with FDR-corrected significance testing
- Joy plot visualization of state-dependent signal flow changes

## Quick Start

```matlab
cd step3_state_hierarchies
run('stage1_main.m')    % Figure 6 components
```

**Expected Output**: Figure 6 panels showing state-specific hierarchy maps, cortical zone comparisons, 27-module signal flow visualizations, and joy plots with significance markers.

## Key Analyses

**State Comparison**: Hierarchy levels across cortical zones with error bars and significance testing
**Signal Flow**: 27-module directed connectivity patterns for each brain state
**Statistical Testing**: FDR correction with multiple significance thresholds (p<0.1, p<0.05, FDR<0.05)

## Output

- **Figure 6a**: State-specific hierarchy surface plots (movie, pain)
- **Figure 6b**: Cortical zone hierarchy comparisons with error bars
- **Figure 6c**: 27-module signal flow visualizations for pain and movie states
- **Joy plots**: State difference visualization with significance markers

Demonstrates state-dependent reorganization of cortical hierarchy and signal flow patterns, revealing how brain states modulate hierarchical organization.