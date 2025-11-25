#!/bin/bash
# BTEH Project: Quick Deployment Reference
# Choose your study and run the appropriate section

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  BTEH Project: Dual-Study Deployment                          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Choose your study:"
echo "  [1] Remote Sensing (HPC) - 123 packages, no system libs"
echo "  [2] EMS Paper (Local)    - 165 packages, needs MPFR/GSL"
echo ""
read -p "Enter choice (1 or 2): " choice

case $choice in
  1)
    echo ""
    echo "═══ Remote Sensing Deployment (HPC) ═══"
    echo ""
    echo "Step 1: Use HPC lock file"
    cp renv_remotesensing.lock renv.lock
    echo "✓ Copied renv_remotesensing.lock → renv.lock"
    
    echo ""
    echo "Step 2: Fix stringi ICU issue"
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
    echo "✓ Set LD_LIBRARY_PATH"
    
    echo ""
    echo "Step 3: Install stringi separately"
    R --vanilla << 'EOF'
install.packages("stringi", 
                 configure.args="--disable-pkg-config",
                 repos="https://cloud.r-project.org")
quit(save = "no")
EOF
    
    echo ""
    echo "Step 4: Restore all packages"
    R -e "renv::restore()"
    
    echo ""
    echo "✅ Remote Sensing environment ready!"
    echo "   Scripts: scripts_RemoteSensing for E and C_/"
    ;;
    
  2)
    echo ""
    echo "═══ EMS Paper Deployment (Local) ═══"
    echo ""
    echo "⚠️  System libraries required:"
    echo "   - libmpfr-dev"
    echo "   - libgmp-dev"
    echo "   - libgsl-dev"
    echo ""
    read -p "Are system libraries installed? (y/n): " libs_ok
    
    if [[ $libs_ok != "y" ]]; then
      echo ""
      echo "Install system libraries first:"
      echo "  Ubuntu/Debian: sudo apt-get install libmpfr-dev libgmp-dev libgsl-dev"
      echo "  RedHat/CentOS: sudo yum install mpfr-devel gmp-devel gsl-devel"
      echo ""
      exit 1
    fi
    
    echo ""
    echo "Step 1: Use EMS lock file"
    cp renv_ems.lock renv.lock
    echo "✓ Copied renv_ems.lock → renv.lock"
    
    echo ""
    echo "Step 2: Restore all packages"
    R -e "renv::restore()"
    
    echo ""
    echo "✅ EMS Paper environment ready!"
    echo "   Scripts: scripts_ems_paper/"
    ;;
    
  *)
    echo "Invalid choice. Exiting."
    exit 1
    ;;
esac

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  Deployment Complete!                                          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
