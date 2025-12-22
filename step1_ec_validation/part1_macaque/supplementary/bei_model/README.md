# BEI Model - Brain E-to-E Inhibition Model

Simplified brain network model with excitatory-to-excitatory connections only, FIC control, and variable Balloon-Windkessel parameters.

## Quick Start

```bash
# Activate conda environment
conda activate hbnm

# Compile
make

# Test run
make test

# Run simulation
./bin/bei_model config/param_default.txt FLN40 12345 results/
```

## Usage

```bash
./bin/bei_model <param_file> <subject_id> <random_seed> <results_dir> [bw_file] [bw_variability]
```

**Parameters:**
- `param_file`: Parameter configuration (see `config/param_default.txt`)
- `subject_id`: Subject identifier (e.g., "FLN40") 
- `random_seed`: Random seed for reproducibility
- `results_dir`: Output directory for results
- `bw_file`: Optional BW parameter file (use "" for default)
- `bw_variability`: BW parameter variability (0.0-1.0, default 0.15)

## Output Files

- `{subject_id}_BOLD_timeseries.txt`: 40×490 BOLD signal matrix
- `{subject_id}_FIC_parameters.txt`: Final FIC inhibition strengths

## Configuration

Edit `config/param_default.txt`:
```
40 1.30 0.15 0.0 0.00316228 490000 1000
```
Format: `regions global_coupling exc_coupling init_inhibition noise timesteps TR`

## Directory Structure

```
bei_model/
├── src/           # Source code
├── input/         # Input connectivity data
├── config/        # Parameter files  
├── bin/           # Compiled executable
├── results/       # Output files
└── Makefile       # Build system
```