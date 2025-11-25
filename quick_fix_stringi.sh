#!/bin/bash
# QUICK FIX for stringi ICU error on HPC
# Copy and paste this entire block into your HPC terminal

echo "🔧 Quick Fix: Installing stringi with system ICU..."

# Set library paths
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH

# Install stringi separately
R --vanilla << 'EOF'
install.packages("stringi", 
                 configure.args="--disable-pkg-config",
                 repos="https://cloud.r-project.org")
quit(save = "no")
EOF

echo "✅ stringi installed. Now run: R -e 'renv::restore()'"
