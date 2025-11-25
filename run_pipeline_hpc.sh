#!/bin/bash
# BTEH Pipeline Runner for HPC
# Handles folder names with spaces correctly

set -e  # Exit on error

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  BTEH Pipeline - HPC Execution                                 ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Define script directory (handles spaces)
SCRIPT_DIR="scripts_RemoteSensing for E and C_"

# Verify directory exists
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "ERROR: Script directory not found: $SCRIPT_DIR"
    echo "Current contents:"
    ls -la
    exit 1
fi

echo "✓ Found script directory: $SCRIPT_DIR"
echo ""

# Function to run script with proper quoting
run_script() {
    local script_name="$1"
    local args="${@:2}"
    
    echo "═══════════════════════════════════════════════════════════════"
    echo "Running: $script_name $args"
    echo "═══════════════════════════════════════════════════════════════"
    
    Rscript "$SCRIPT_DIR/$script_name" $args
    
    if [ $? -eq 0 ]; then
        echo "✓ Completed: $script_name"
    else
        echo "✗ FAILED: $script_name"
        exit 1
    fi
    echo ""
}

# H2O Training
echo "STAGE 1: H2O Model Training"
echo "────────────────────────────────────────────────────────────────"
run_script "03_h2o_train.R" --run A --mode REPRO
run_script "03_h2o_train.R" --run B --mode REPRO

# SSDM Training
echo "STAGE 2: SSDM Model Training"
echo "────────────────────────────────────────────────────────────────"
run_script "04_ssdm_train.R" --run A --mode REPRO
run_script "04_ssdm_train.R" --run B --mode REPRO

# Comparisons
echo "STAGE 3: Model Comparisons"
echo "────────────────────────────────────────────────────────────────"
run_script "05_h2o_vs_ssdm.R"
run_script "05b_temporal_comparison.R"
run_script "05a_h2o_vs_ssdm_panel.R"

# Appendix
echo "STAGE 4: Appendix Figures"
echo "────────────────────────────────────────────────────────────────"
run_script "07_appendix.R"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  Pipeline Complete!                                            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "End time: $(date)"
echo ""
echo "Results saved in: results/"
echo "Logs saved in: logs/"
