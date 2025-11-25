# 🎯 QUICK START: Your First RStudio Session on HPC

## Copy-Paste Commands (In Order)

### 1️⃣ In RStudio Console (First Time Setup)
```r
# Set working directory
setwd("~/BTEH/BTEH")

# Check environment
renv::status()

# Load packages to verify
library(terra)
library(h2o)
library(dplyr)
library(ggplot2)

cat("\n✅ Ready to run!\n")
```

### 2️⃣ Navigate to Scripts
```r
setwd("scripts_RemoteSensing for E and C_")
getwd()  # Should show: /home/harin/BTEH/BTEH/scripts_RemoteSensing for E and C_
```

### 3️⃣ Run Your First Test
```r
# Quick test (5-10 minutes)
source("03_h2o_train.R")
```

### 4️⃣ Check Results
```r
setwd("~/BTEH/BTEH")
list.files("results/H2O/", recursive = TRUE)
```

---

## 📂 Key Directories

```r
# Project root
setwd("~/BTEH/BTEH")

# Scripts
setwd("~/BTEH/BTEH/scripts_RemoteSensing for E and C_")

# Results
list.files("~/BTEH/BTEH/results/")

# Logs
list.files("~/BTEH/BTEH/logs/")
```

---

## 🚀 Full Pipeline (After Test Works)

```r
# From scripts directory
setwd("~/BTEH/BTEH/scripts_RemoteSensing for E and C_")

# Run all scripts in order
source("03_h2o_train.R")
source("04_ssdm_train.R")
source("05_h2o_vs_ssdm.R")
source("05b_temporal_comparison.R")
source("05a_h2o_vs_ssdm_panel.R")
source("07_appendix.R")
```

---

## 🆘 If Something Goes Wrong

```r
# Check where you are
getwd()

# Go back to project root
setwd("~/BTEH/BTEH")

# Check renv
renv::status()

# Restart R session
# Session → Restart R

# Try again
```

---

**That's it! Start with the test, then run the full pipeline.** 🎉
