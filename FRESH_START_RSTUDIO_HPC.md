# 🚀 Fresh Start: RStudio Server on HPC

## Quick Start Guide for Clean GitHub → HPC Deployment

---

## 📋 Prerequisites Checklist

Before you start, make sure you have:
- [ ] Pushed all changes to GitHub
- [ ] RStudio Server URL from HPC admin
- [ ] SSH access to HPC
- [ ] Your HPC username and password

---

## 🎯 STEP-BY-STEP GUIDE

### **Step 1: Clean Up Old Trial Directory (SSH)**

```bash
# SSH into HPC
ssh harin@eslab-hpc03

# Remove old trial directory
rm -rf ~/trial_BTEH

# Verify it's gone
ls ~ | grep trial
```

---

### **Step 2: Clone Fresh from GitHub (SSH)**

```bash
# Clone your repository
cd ~
git clone https://github.com/Caspian-Aiyanna/BTEH_for_HPC.git BTEH

# Navigate into project
cd BTEH

# Verify structure
ls -la
# Should see: BTEH/, renv.lock, renv_remotesensing.lock, etc.

# Check the BTEH subfolder
ls -la BTEH/
# Should see: scripts_RemoteSensing for E and C_/, data/, R/, etc.
```

---

### **Step 3: Set Up R Environment (SSH)**

```bash
# Load R module (adjust version as needed)
module load R/4.5.1

# Set library path for stringi fix
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Install stringi separately (ICU fix)
R --vanilla << 'EOF'
install.packages("stringi", 
                 configure.args="--disable-pkg-config",
                 repos="https://cloud.r-project.org")
quit(save = "no")
EOF
```

**Expected output**: stringi should install successfully

---

### **Step 4: Restore R Packages (SSH)**

```bash
# Still in ~/BTEH directory
R -e "renv::restore()"
```

**This will take 15-30 minutes**. You should see:
- Installing 123 packages
- No errors about MPFR, GSL, or Rmpfr
- All packages install successfully

---

### **Step 5: Open RStudio Server**

1. **Open your browser** and navigate to RStudio Server URL:
   - Usually: `https://rstudio.hpc.university.edu`
   - Or: `http://eslab-hpc03:8787`

2. **Login** with your HPC credentials:
   - Username: `harin`
   - Password: [your HPC password]

3. **You should see the RStudio interface**

---

### **Step 6: Set Working Directory in RStudio**

In the **RStudio Console**, run:

```r
# Navigate to project
setwd("~/BTEH/BTEH")

# Verify you're in the right place
getwd()
# Should show: /home/harin/BTEH/BTEH

# List files
list.files()
# Should see: config.yml, R/, scripts_RemoteSensing for E and C_/, data/, etc.
```

---

### **Step 7: Verify Environment in RStudio**

```r
# Check renv status
renv::status()
# Should show: "No issues found" or "The library is already synchronized"

# Test key packages
library(terra)
library(sf)
library(h2o)
library(SSDM)
library(dplyr)
library(ggplot2)

# All should load without errors
cat("\n✅ All packages loaded successfully!\n")
```

---

### **Step 8: Test H2O Initialization**

```r
# Initialize H2O
library(h2o)
h2o.init(nthreads = 1, max_mem_size = "4G")

# Should see H2O cluster info
# Then shutdown
h2o.shutdown(prompt = FALSE)
```

---

### **Step 9: Verify Data Files**

```r
# Check data directories
list.files("data/clean/OG")
list.files("data/envi/A")
list.files("data/envi/B")
list.files("data/shp")

# All should show your data files
```

---

### **Step 10: Run Your First Test**

#### **Option A: Quick Test (Recommended First)**

```r
# Navigate to scripts directory
setwd("scripts_RemoteSensing for E and C_")

# Source the script with test parameters
source("03_h2o_train.R")

# Or run from terminal in RStudio
system('Rscript 03_h2o_train.R --run B --mode FAST --species E3B')
```

#### **Option B: Full Pipeline**

```r
# Go back to project root
setwd("~/BTEH/BTEH")

# Navigate to scripts
setwd("scripts_RemoteSensing for E and C_")

# Run all scripts in order
source("03_h2o_train.R")  # H2O training
source("04_ssdm_train.R")  # SSDM training
source("05_h2o_vs_ssdm.R")  # Comparison
source("05b_temporal_comparison.R")  # Temporal analysis
source("05a_h2o_vs_ssdm_panel.R")  # Panels
source("07_appendix.R")  # Appendix
```

---

## 🎯 RECOMMENDED WORKFLOW

### **For Interactive Work in RStudio:**

1. **Set working directory to scripts folder:**
```r
setwd("~/BTEH/BTEH/scripts_RemoteSensing for E and C_")
```

2. **Run scripts one at a time:**
```r
source("03_h2o_train.R")
```

3. **Monitor progress in Console**

4. **Check results:**
```r
setwd("~/BTEH/BTEH")
list.files("results/H2O/B/", recursive = TRUE)
```

---

### **For Long-Running Jobs (Terminal in RStudio):**

1. **Open Terminal** in RStudio (Tools → Terminal → New Terminal)

2. **Navigate to scripts:**
```bash
cd ~/BTEH/BTEH/"scripts_RemoteSensing for E and C_"
```

3. **Run script:**
```bash
Rscript 03_h2o_train.R --run B --mode REPRO
```

4. **Monitor in real-time** - output will show in terminal

---

## 📊 Expected Timeline

| Task | Time |
|------|------|
| Clone from GitHub | 1-2 min |
| Install stringi | 2-5 min |
| renv::restore() | 15-30 min |
| Test run (FAST mode, 1 species) | 5-10 min |
| Full pipeline (REPRO mode) | 8-12 hours |

---

## ✅ Success Indicators

You'll know everything is working when:

1. ✅ `renv::status()` shows no issues
2. ✅ All packages load without errors
3. ✅ H2O initializes successfully
4. ✅ Test script completes without errors
5. ✅ Results appear in `results/` folder

---

## 🐛 Quick Troubleshooting

### **Issue: renv::restore() fails**
```r
# Update renv first
install.packages("renv")
renv::restore()
```

### **Issue: H2O won't start**
```r
# Check Java
system("java -version")

# Try with less memory
h2o.init(nthreads = 1, max_mem_size = "2G")
```

### **Issue: Can't find data files**
```r
# Check you're in right directory
getwd()
# Should be: /home/harin/BTEH/BTEH

# List data
list.files("data", recursive = TRUE)
```

### **Issue: Script path errors**
```r
# Always run scripts from within their directory
setwd("~/BTEH/BTEH/scripts_RemoteSensing for E and C_")
source("03_h2o_train.R")
```

---

## 💡 Pro Tips

1. **Use RStudio Projects**: File → Open Project → Navigate to BTEH.Rproj
2. **Save your workspace**: Session → Save Workspace As...
3. **Use Terminal for long runs**: Keeps RStudio responsive
4. **Monitor disk space**: `system("df -h ~")`
5. **Check logs**: `readLines("../../logs/03_h2o_train_B.log")`

---

## 📁 Directory Structure Reference

```
~/BTEH/
├── BTEH/                                    # ← Main project folder
│   ├── scripts_RemoteSensing for E and C_/ # ← Your scripts HERE
│   │   ├── 03_h2o_train.R
│   │   ├── 04_ssdm_train.R
│   │   ├── 05_h2o_vs_ssdm.R
│   │   ├── 05a_h2o_vs_ssdm_panel.R
│   │   ├── 05b_temporal_comparison.R
│   │   └── 07_appendix.R
│   ├── data/                                # ← Your data
│   ├── R/                                   # ← Helper functions
│   ├── results/                             # ← Output goes here
│   ├── logs/                                # ← Log files
│   └── config.yml                           # ← Configuration
├── renv.lock                                # ← Package dependencies
├── renv_remotesensing.lock                  # ← HPC version
└── .Rprofile                                # ← Activates renv
```

---

## 🎬 YOUR FIRST COMMAND IN RSTUDIO

```r
# Copy-paste this entire block into RStudio Console:

# Set working directory
setwd("~/BTEH/BTEH")

# Verify environment
cat("Working directory:", getwd(), "\n")
cat("renv status:\n")
renv::status()

# Load key packages
library(terra)
library(h2o)
library(dplyr)

cat("\n✅ Environment ready!\n")
cat("Next: Navigate to scripts and run your analysis\n")
cat("Command: setwd('scripts_RemoteSensing for E and C_')\n")
```

---

**You're all set! Start with the test run, then move to full pipeline.** 🚀

**Good luck with your HPC runs!** 🐘📊
