#!/usr/bin/env Rscript
# Pre-flight check before overnight run

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  BTEH: Pre-Flight Check for Overnight Run                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Check 1: Working directory
cat("Check 1: Working directory\n")
wd <- getwd()
cat("  Current:", wd, "\n")
if (grepl("BTEH$", wd)) {
  cat("  ✓ Correct directory\n\n")
} else {
  cat("  ✗ Wrong directory! Should end with 'BTEH'\n")
  cat("  Run: setwd('~/BTEH/BTEH')\n\n")
}

# Check 2: renv status
cat("Check 2: renv environment\n")
status <- try(renv::status(), silent = TRUE)
if (!inherits(status, "try-error")) {
  cat("  ✓ renv is active\n\n")
} else {
  cat("  ✗ renv issue detected\n")
  cat("  Run: renv::restore()\n\n")
}

# Check 3: Key packages
cat("Check 3: Key packages\n")
pkgs <- c("terra", "sf", "h2o", "SSDM", "dplyr", "ggplot2")
all_ok <- TRUE
for (pkg in pkgs) {
  if (require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "NOT FOUND\n")
    all_ok <- FALSE
  }
}
if (all_ok) cat("\n  ✓ All packages available\n\n")

# Check 4: H2O initialization
cat("Check 4: H2O initialization\n")
h2o_ok <- try({
  library(h2o)
  h2o.init(nthreads = 1, max_mem_size = "2G", startH2O = TRUE)
  h2o.shutdown(prompt = FALSE)
  TRUE
}, silent = TRUE)

if (!inherits(h2o_ok, "try-error")) {
  cat("  ✓ H2O can initialize\n\n")
} else {
  cat("  ✗ H2O initialization failed\n")
  cat("  Check Java installation\n\n")
}

# Check 5: Data files
cat("Check 5: Data files\n")
data_dirs <- c("data/clean/OG", "data/envi/A", "data/envi/B", "data/shp")
for (dir in data_dirs) {
  if (dir.exists(dir)) {
    n_files <- length(list.files(dir))
    cat("  ✓", dir, "-", n_files, "files\n")
  } else {
    cat("  ✗", dir, "NOT FOUND\n")
  }
}
cat("\n")

# Check 6: Scripts exist
cat("Check 6: Scripts\n")
scripts <- c(
  "scripts_RemoteSensing for E and C_/03_h2o_train.R",
  "scripts_RemoteSensing for E and C_/04_ssdm_train.R"
)
for (script in scripts) {
  if (file.exists(script)) {
    cat("  ✓", basename(script), "\n")
  } else {
    cat("  ✗", basename(script), "NOT FOUND\n")
  }
}
cat("\n")

# Check 7: Disk space
cat("Check 7: Disk space\n")
if (.Platform$OS.type == "windows") {
  disk_info <- system("wmic logicaldisk get size,freespace,caption", intern = TRUE)
  cat("  ", disk_info[3], "\n")
  cat("  (Need at least 5 GB free)\n\n")
} else {
  system("df -h .")
  cat("\n")
}

# Summary
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  Pre-Flight Check Complete                                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("If all checks passed, you're ready to run:\n")
cat("  ./run_overnight.bat\n\n")

cat("Estimated completion time: 10-14 hours\n")
cat("Start now to have results by morning!\n\n")
