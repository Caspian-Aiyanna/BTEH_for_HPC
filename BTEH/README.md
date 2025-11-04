# 🐘 BTEH – Elephant Habitat Modelling & Uncertainty Analysis

This repository contains the full, reproducible modelling pipeline for **Main Paper** of the BTEH project — quantifying **habitat suitability, temporal change, and algorithmic uncertainty** for African elephants before and after fence removal in Kariega Game Reserve, South Africa.

All workflows are implemented in **R** (with optional Google Earth Engine extraction scripts) and are designed to run both locally and on **HPC clusters** via SLURM or PBS job schedulers.

---

## 🌿 Project Overview

The pipeline integrates multiple modelling frameworks:

| Framework | Description |
|------------|--------------|
| **H2O AutoML** | Machine-learning ensemble (GBM, XGBoost, DNN, DRF, etc - upto 100 models) for habitat suitability |
| **SSDM** | Classical stacked species distribution modelling (GLM, GAM, RF, ANN, GBM) |
| **SSF/RSF** | Step-selection functions derived from GPS telemetry (used for triangulation) |
| **Comparison & Uncertainty** | Between-method, temporal, and replicate analyses |

All outputs are automatically written under `results/` and logged under `logs/`.

---

## 📁 Folder Structure

```
.
├── config.yml                # Global config (mode, cores, etc.)
├── Makefile                  # Optional automation targets
├── renv.lock                 # Reproducible R environment
├── R/                        # Utility functions
│   ├── utils_h2o.R
│   ├── utils_io.R
│   ├── utils_kendall.R
│   ├── utils_plot.R
│   └── utils_repro.R
├── scripts/                  # Main pipeline scripts
│   ├── 02_dbscan_thin_degrees.R     # Occurrence thinning
│   ├── 03_h2o_train.R               # AutoML training
│   ├── 04_ssdm_train.R              # SSDM ensemble training
│   ├── 05_h2o_vs_ssdm.R             # Between-method comparison
│   ├── 06a_uncertainty.R            # Uncertainty decomposition
│   └── BTEH_GEE_Extract.js          # GEE environmental extraction
├── data/
│   ├── clean/       # Pre-processed telemetry CSVs (E1B–E6B, A/B)
│   ├── envi/        # Environmental stacks (A/B)
│   ├── occ/         # Thinned & replicate occurrences
│   └── shp/         # AOI shapefiles (HV20233.*)
├── results/          # Model outputs (H2O, SSDM, compare, uncertainty, etc.)
├── logs/             # Run-time logs
├── plans/            # Variable-selection plans (Kendall results)
└── hpc/              # SLURM job scripts
```

---

## ⚙️ Environment Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/Caspian-Aiyanna/ems_paper.git
   cd ems_paper
   ```

2. **Restore the R environment**
   ```bash
   Rscript -e "renv::restore(prompt = FALSE)"
   ```

3. **(Optional)** Check directory structure
   ```r
   fs::dir_tree('.', recurse = 2)
   ```

---

## 🚀 Running the Pipeline

### 🧩 Option 1 — Reproducible (“REPRO”) Run

Use **single-core deterministic** execution to eliminate randomness.  
Ideal for publication and final validation.

```bash
sbatch hpc/BTEH_REPRO.slurm
```

Or locally:

```bash
Rscript scripts/03_h2o_train.R --run A --mode REPRO
Rscript scripts/04_ssdm_train.R --run A --mode REPRO
Rscript scripts/05_h2o_vs_ssdm.R --run A --mode REPRO

Rscript scripts/03_h2o_train.R --run B --mode REPRO
Rscript scripts/04_ssdm_train.R --run B --mode REPRO
Rscript scripts/05_h2o_vs_ssdm.R --run B --mode REPRO
```

**Outputs:**  
- `results/H2O/<RUN>/<SP>/prediction_<SP>.tif`  
- `results/SSDM/<RUN>/<SP>/ESDM_<SP>.tif`  
- `results/compare/h2o_vs_ssdm/<RUN>/metrics/…`  
- Logs under `logs/03_*.log`, `logs/04_*.log`, `logs/05_*.log`

---

### ⚡ Option 2 — Fast (“FAST”) Run

Parallelized execution for development and testing.  
Results are near-identical but may vary slightly due to parallel randomness.

```bash
sbatch hpc/BTEH_FAST.slurm
```

Or locally (multi-core machine):

```bash
Rscript scripts/03_h2o_train.R --run A --mode FAST
Rscript scripts/04_ssdm_train.R --run A --mode FAST
Rscript scripts/05_h2o_vs_ssdm.R --run A --mode FAST
```

---

### 🧬 Option 3 — Species-Array Mode (HPC only)

Runs each elephant (E1B–E6B) as a separate array job:

```bash
sbatch hpc/BTEH_FAST_array.slurm
```

This distributes species across nodes and merges results automatically.

---

## 🧾 Log Files

All scripts log progress and warnings to the `logs/` folder:
- `03_h2o_train_A.log`, `04_ssdm_train_B.log`, etc.
- `05_compare_methods_A.log` and `05_compare_methods_B.log`
- Environment information and timing details are captured for reproducibility.

---

## 📊 Expected Outputs

- `results/H2O/` – AutoML rasters, models, leaderboards  
- `results/SSDM/` – Ensemble rasters, algorithm summaries  
- `results/compare/` – Metrics, hotspot overlaps, maps  
- `results/uncertainty/` – Variance, stability, gain/loss tables  
- `results/figures/` – Final patchwork panels for publication

Each result folder includes `.csv` summaries and `.tif` rasters ready for visualization.

---

## 🧠 Notes for HPC Users

- **Modules:** Load `R/4.3.2`, `GDAL`, `GEOS`, and `PROJ` if required.  
- **Logs:** Check progress with `tail -f logs/BTEH_REPRO_*.out`.  
- **Monitoring:** Use `squeue -u $USER` (SLURM) or `qstat` (PBS).  
- **Storage:** Prefer `/scratch` for heavy intermediate files.

---

## 📚 Citation

If you use this workflow or data structure, please cite:

> Aiyanna C., et al. (2025). *Automated ensemble modelling and uncertainty quantification for African elephant habitat connectivity*.  
> Environmental Modelling & Software. (In preparation)

---

## 🧰 Contact

**Author:** Caspian Aiyanna  
**Institution:** [Your University / Research Group]  
**Email:** [harinaiyanna.cheriyandaraveendra@phd.unipd.it]  
**GitHub:** [https://github.com/Caspian-Aiyanna](https://github.com/Caspian-Aiyanna)

---

### 🏁 Quick Summary

| Mode | Purpose | HPC script | CPU | Typical runtime |
|------|----------|-------------|-----|-----------------|
| **REPRO** | Final, deterministic results | `hpc/BTEH_REPRO.slurm` | 1 | 24–48 h |
| **FAST** | Development, parallel | `hpc/BTEH_FAST.slurm` | 8–16 | 3–8 h |
| **FAST-ARRAY** | Multi-species parallel | `hpc/BTEH_FAST_array.slurm` | 8 each | 1–4 h |

---

*This repository embodies the principles of open, reproducible ecological modelling and provides a modular foundation for cross-framework SDM benchmarking.*
