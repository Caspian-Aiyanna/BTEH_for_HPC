#!/bin/bash
# HPC Environment Setup Script for BTEH Project
# This script handles problematic package installations on HPC

echo "=========================================="
echo "BTEH HPC Environment Setup"
echo "=========================================="

# Step 1: Load R module (adjust version as needed)
echo "Step 1: Loading R module..."
module load R/4.5.1 || module load R || echo "Warning: Could not load R module. Continuing..."

# Step 2: Set library paths for ICU
echo "Step 2: Setting library paths..."
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

# Step 3: Backup original renv.lock
echo "Step 3: Backing up renv.lock..."
if [ -f "renv.lock" ] && [ ! -f "renv_original.lock" ]; then
    cp renv.lock renv_original.lock
    echo "✓ Backed up to renv_original.lock"
fi

# Step 4: Use HPC-compatible lock file
echo "Step 4: Using HPC-compatible lock file..."
if [ -f "renv_hpc.lock" ]; then
    cp renv_hpc.lock renv.lock
    echo "✓ Using renv_hpc.lock (without amt/ctmm/Bessel/Rmpfr)"
else
    echo "⚠ Warning: renv_hpc.lock not found, using existing renv.lock"
fi

# Step 5: Install stringi separately (to avoid ICU issues)
echo "Step 5: Installing stringi with system ICU..."
Rscript -e "if (!require('stringi', quietly = TRUE)) { install.packages('stringi', configure.args='--disable-pkg-config', repos='https://cloud.r-project.org') }" || {
    echo "⚠ Warning: stringi installation failed. Trying alternative method..."
    Rscript -e "install.packages('stringi', repos='https://cloud.r-project.org')"
}

# Step 6: Restore remaining packages with renv
echo "Step 6: Restoring packages with renv..."
Rscript -e "renv::restore()" && {
    echo "=========================================="
    echo "✓ SUCCESS: Environment setup complete!"
    echo "=========================================="
    exit 0
} || {
    echo "=========================================="
    echo "✗ ERROR: renv::restore() failed"
    echo "=========================================="
    echo ""
    echo "Troubleshooting steps:"
    echo "1. Check the error message above"
    echo "2. Try: export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH"
    echo "3. Try: R -e \"install.packages('stringi', configure.args='--disable-pkg-config')\""
    echo "4. See HPC_STRINGI_FIX.md for more solutions"
    exit 1
}
