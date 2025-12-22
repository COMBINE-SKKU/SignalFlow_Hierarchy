# ===============================================================================
# BEI MODEL - Brain E-to-E Inhibition Model MAKEFILE
# ===============================================================================
#
# Simplified brain network model with E-to-E connections only
# Features: FIC control, variable BW parameters, burn-in phase
#
# Requirements:
#   - System OpenBLAS/LAPACKE libraries
#   - GCC compiler
#   - Conda environment 'hbnm'
#
# Usage:
#   make            # Build the executable
#   make clean      # Clean build artifacts
#   make test       # Test compilation and run example
#   make help       # Show help
#
# ===============================================================================

# Project configuration
TARGET = bin/bei_model
SRCDIR = src
SOURCES = $(SRCDIR)/bei_model.c
BINDIR = bin

# Compiler configuration (avoid FSL conflicts)
override PATH := /usr/bin:/bin:/usr/sbin:/sbin:$(PATH)
CC := /usr/bin/gcc
CFLAGS = -O3 -ffast-math -march=native -fopenmp -msse2 -Wall -Wextra
LDFLAGS = -lm -fopenmp

# Default target
all: $(TARGET)

# Create bin directory if it doesn't exist
$(BINDIR):
	@echo "Creating bin directory..."
	@mkdir -p $(BINDIR)

# Build target
$(TARGET): $(SOURCES) | $(BINDIR)
	@echo "==============================================================================="
	@echo "COMPILING BEI MODEL"
	@echo "==============================================================================="
	@echo "Environment: $(CONDA_DEFAULT_ENV)"
	@echo "Compiler: $(shell which $(CC))"
	@echo "Sources: $(SOURCES)"
	@echo "Target: $(TARGET)"
	@echo ""
	$(CC) $(CFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS)
	@echo ""
	@echo "✓ Compilation successful: $(TARGET)"
	@echo "  File size: $(shell ls -lh $(TARGET) | awk '{print $$5}')"
	@echo "  To run: make test"

# Test compilation and run example
test: $(TARGET)
	@echo ""
	@echo "==============================================================================="
	@echo "TESTING BEI MODEL"
	@echo "==============================================================================="
	@echo "Running example simulation with default parameters..."
	@mkdir -p results/test
	@echo "Command: $(TARGET) config/param_default.txt FLN40 12345 results/test"
	@$(TARGET) config/param_default.txt FLN40 12345 results/test "" 0.0
	@if [ -f "results/test/FLN40_BOLD_timeseries.txt" ]; then \
		echo ""; \
		echo "✓ Test simulation successful!"; \
		HEADER=$$(head -n1 results/test/FLN40_BOLD_timeseries.txt); \
		echo "  Output: $$HEADER"; \
		echo "  Files: results/test/FLN40_BOLD_timeseries.txt"; \
		echo "         results/test/FLN40_FIC_parameters.txt"; \
	else \
		echo ""; \
		echo "✗ Test simulation failed!"; \
		exit 1; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	@rm -f $(TARGET)
	@rm -rf results/test
	@if [ -d "$(BINDIR)" ] && [ -z "$$(ls -A $(BINDIR) 2>/dev/null)" ]; then \
		rmdir $(BINDIR); \
		echo "  Removed empty $(BINDIR) directory"; \
	fi
	@echo "Clean complete"

# Show build configuration
config:
	@echo "==============================================================================="
	@echo "BEI MODEL BUILD CONFIGURATION"
	@echo "==============================================================================="
	@echo "Environment:"
	@echo "  Conda environment: $(CONDA_DEFAULT_ENV)"
	@echo "  Conda prefix: $(CONDA_PREFIX)"
	@echo "  Working directory: $(shell pwd)"
	@echo ""
	@echo "Compiler:"
	@echo "  CC: $(shell which $(CC))"
	@echo "  Version: $(shell $(CC) --version | head -n1)"
	@echo "  CFLAGS: $(CFLAGS)"
	@echo "  LDFLAGS: $(LDFLAGS)"
	@echo ""
	@echo "Project:"
	@echo "  Sources: $(SOURCES)"
	@echo "  Target: $(TARGET)"
	@echo "  Source files:"
	@for src in $(SOURCES); do \
		if [ -f "$$src" ]; then \
			echo "    ✓ $$src"; \
		else \
			echo "    ✗ $$src (MISSING)"; \
		fi \
	done
	@echo ""
	@echo "Input files:"
	@ls -1 input/ | head -5 | while read file; do echo "    ✓ input/$$file"; done

# Development build (with debug symbols)
dev: CFLAGS = -O1 -g -Wall -Wextra
dev: $(TARGET)
	@echo "Development build complete (with debug symbols)"

# Release build (maximum optimization)
release: CFLAGS = -O3 -ffast-math -march=native -fopenmp -msse2 -DNDEBUG
release: clean $(TARGET)
	@echo "Release build complete (maximum optimization)"

# Check environment
check-env:
	@echo "Checking environment..."
ifndef CONDA_PREFIX
	@echo "ERROR: No conda environment active"
	@echo "   Run: conda activate hbnm"
	@exit 1
endif
	@echo "✓ Conda environment active: $(CONDA_DEFAULT_ENV)"
	@echo "✓ Conda prefix: $(CONDA_PREFIX)"
	@echo "✓ GCC available: $(shell which gcc)"

# Run multiple simulations script
simulate:
	@echo "Running multiple BEI simulations..."
	@if [ ! -f "scripts/run_simulations.sh" ]; then \
		echo "ERROR: scripts/run_simulations.sh not found"; \
		echo "Create this script for batch simulations"; \
		exit 1; \
	fi
	@chmod +x scripts/run_simulations.sh
	@scripts/run_simulations.sh

# Help target
help:
	@echo "==============================================================================="
	@echo "BEI MODEL - Build System"
	@echo "==============================================================================="
	@echo ""
	@echo "Prerequisites:"
	@echo "  1. System libraries: Standard GCC and math library"
	@echo "  2. Conda environment: conda activate hbnm"
	@echo ""
	@echo "Available targets:"
	@echo "  all          - Build the executable (default)"
	@echo "  test         - Build and run test simulation"
	@echo "  clean        - Remove build artifacts"
	@echo "  config       - Show detailed build configuration"
	@echo "  check-env    - Verify conda environment"
	@echo "  dev          - Development build (with debug symbols)"
	@echo "  release      - Release build (maximum optimization)"
	@echo "  simulate     - Run batch simulations (requires scripts/run_simulations.sh)"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Example workflow:"
	@echo "  conda activate hbnm          # Activate environment"
	@echo "  make                         # Build executable"
	@echo "  make test                    # Test with example"
	@echo "  ./bin/bei_model config/param_default.txt FLN40 12345 results/"
	@echo ""
	@echo "Model parameters (config/param_default.txt):"
	@echo "  Format: regions G exc_coupling init_inh noise timesteps TR"
	@echo "  Example: 40 1.30 0.15 0.0 0.00316228 490000 1000"

.PHONY: all test clean config dev release check-env simulate help