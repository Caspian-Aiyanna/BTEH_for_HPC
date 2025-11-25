#!/usr/bin/env Rscript
# Comprehensive dependency analysis for both studies
# Creates two separate renv.lock files optimized for each study

library(jsonlite)

# Load original lock file
lock_data <- fromJSON("renv.lock", simplifyVector = FALSE)

cat("=== BTEH Project: Dual-Study Dependency Analysis ===\n\n")

# ============================================================================
# STUDY 1: Remote Sensing (E and C) - For HPC
# ============================================================================
cat("STUDY 1: Remote Sensing for E and C\n")
cat("Location: scripts_RemoteSensing for E and C_/\n")
cat("Purpose: H2O and SSDM model comparison\n\n")

# Core packages used in RemoteSensing scripts (from grep analysis)
remotesensing_core <- c(
    # Data manipulation
    "dplyr", "tidyr", "readr", "tibble", "purrr", "stringr",
    # Spatial
    "terra", "sf", "raster",
    # Modeling
    "h2o", "SSDM",
    # Visualization
    "ggplot2", "scales", "grid",
    # Utilities
    "tools", "lubridate"
)

# System library issues to remove for HPC
hpc_problematic <- c("amt", "ctmm", "Bessel", "Rmpfr", "gsl")

cat("Core packages needed:", length(remotesensing_core), "\n")
cat("Removing HPC-problematic packages:", length(hpc_problematic), "\n\n")

# ============================================================================
# STUDY 2: EMS Paper - Movement Ecology
# ============================================================================
cat("STUDY 2: EMS Paper (Movement Ecology)\n")
cat("Location: scripts_ems_paper/\n")
cat("Purpose: SSF analysis, DBSCAN thinning, replicates\n\n")

# Core packages used in EMS paper scripts
ems_core <- c(
    # Movement ecology (CRITICAL - needs system libraries)
    "amt", # Animal movement tools
    "survival", # Cox models for SSF
    # Spatial
    "terra", "sf", "raster",
    # Data manipulation
    "dplyr", "tidyr", "readr", "tibble", "purrr", "stringr",
    # Modeling
    "h2o", "SSDM",
    # Visualization
    "ggplot2", "scales", "grid", "cowplot", "viridisLite",
    # Mapping
    "rnaturalearth",
    # Utilities
    "tools", "lubridate", "optparse"
)

cat("Core packages needed:", length(ems_core), "\n")
cat("Includes amt (requires system libs: MPFR, GSL)\n\n")

# ============================================================================
# Create Lock File 1: Remote Sensing (HPC-compatible)
# ============================================================================
cat("Creating: renv_remotesensing.lock (HPC-compatible)\n")

# Get all dependencies recursively
get_all_deps <- function(pkg_names, lock_data) {
    all_pkgs <- character(0)
    to_check <- pkg_names

    while (length(to_check) > 0) {
        pkg <- to_check[1]
        to_check <- to_check[-1]

        if (pkg %in% all_pkgs || !(pkg %in% names(lock_data$Packages))) {
            next
        }

        all_pkgs <- c(all_pkgs, pkg)
        pkg_info <- lock_data$Packages[[pkg]]

        # Get dependencies
        deps <- character(0)
        if (!is.null(pkg_info$Imports)) deps <- c(deps, pkg_info$Imports)
        if (!is.null(pkg_info$Depends)) deps <- c(deps, pkg_info$Depends)
        if (!is.null(pkg_info$LinkingTo)) deps <- c(deps, pkg_info$LinkingTo)

        # Clean dependency names (remove version specs)
        deps <- gsub("\\s*\\(.*\\)", "", deps)
        deps <- deps[!deps %in% c("R", "methods", "stats", "utils", "graphics", "grDevices", "grid")]

        to_check <- c(to_check, setdiff(deps, all_pkgs))
    }

    all_pkgs
}

rs_all_deps <- get_all_deps(remotesensing_core, lock_data)
rs_all_deps <- setdiff(rs_all_deps, hpc_problematic)

# Create RemoteSensing lock file
rs_lock <- lock_data
rs_lock$Packages <- lock_data$Packages[rs_all_deps]

write_json(rs_lock, "renv_remotesensing.lock", pretty = TRUE, auto_unbox = TRUE)
cat("✓ Saved:", length(rs_lock$Packages), "packages\n")
cat("✓ Removed:", paste(hpc_problematic, collapse = ", "), "\n\n")

# ============================================================================
# Create Lock File 2: EMS Paper (Full dependencies)
# ============================================================================
cat("Creating: renv_ems.lock (Full dependencies including amt)\n")

ems_all_deps <- get_all_deps(ems_core, lock_data)

# Create EMS lock file
ems_lock <- lock_data
ems_lock$Packages <- lock_data$Packages[ems_all_deps]

write_json(ems_lock, "renv_ems.lock", pretty = TRUE, auto_unbox = TRUE)
cat("✓ Saved:", length(ems_lock$Packages), "packages\n")
cat("✓ Includes amt and dependencies (requires MPFR, GSL on system)\n\n")

# ============================================================================
# Summary Report
# ============================================================================
cat("=== SUMMARY ===\n\n")

cat("Original renv.lock:          ", length(lock_data$Packages), "packages\n")
cat("renv_remotesensing.lock:     ", length(rs_lock$Packages), "packages (HPC-ready)\n")
cat("renv_ems.lock:               ", length(ems_lock$Packages), "packages (needs system libs)\n\n")

# Packages unique to each study
rs_only <- setdiff(rs_all_deps, ems_all_deps)
ems_only <- setdiff(ems_all_deps, rs_all_deps)
shared <- intersect(rs_all_deps, ems_all_deps)

cat("Shared packages:             ", length(shared), "\n")
cat("RemoteSensing-only:          ", length(rs_only), "\n")
if (length(rs_only) > 0) cat("  -", paste(head(rs_only, 10), collapse = ", "), "\n")
cat("EMS-only:                    ", length(ems_only), "\n")
if (length(ems_only) > 0) cat("  -", paste(head(ems_only, 10), collapse = ", "), "\n\n")

# System library requirements
cat("=== SYSTEM LIBRARY REQUIREMENTS ===\n\n")
cat("renv_remotesensing.lock:\n")
cat("  ✓ No system libraries required\n")
cat("  ✓ Ready for HPC deployment\n\n")

cat("renv_ems.lock:\n")
cat("  ⚠ Requires system libraries:\n")
cat("    - libmpfr-dev (for Rmpfr)\n")
cat("    - libgmp-dev (for gmp)\n")
cat("    - libgsl-dev (for gsl)\n")
cat("  ⚠ Contact HPC admin or install locally\n\n")

cat("=== USAGE ===\n\n")
cat("For Remote Sensing (HPC):\n")
cat("  cp renv_remotesensing.lock renv.lock\n")
cat("  renv::restore()\n\n")

cat("For EMS Paper (Local or HPC with system libs):\n")
cat("  cp renv_ems.lock renv.lock\n")
cat("  renv::restore()\n\n")

cat("✓ Analysis complete!\n")
