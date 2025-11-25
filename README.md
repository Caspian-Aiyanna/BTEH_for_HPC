# BTEH Project: Elephant Habitat Modeling

Species Distribution Models (SDMs) for Asian elephants using H2O AutoML and SSDM ensemble methods.

## Requirements

- R 4.4.2 or higher
- RStudio (recommended)
- Java 8-17 (for H2O)
- 8+ GB RAM
- 5+ GB disk space

## Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/Caspian-Aiyanna/BTEH_for_HPC.git
cd BTEH_for_HPC
```

### 2. Install Packages

Open RStudio and open the project file `BTEH.Rproj`. Then run:

```r
renv::restore()
```

This installs all required packages. Takes 15-30 minutes.

### 3. Run Analysis

In RStudio console:

```r
# Set working directory to project root
setwd("C:/path/to/BTEH/BTEH")

# Run H2O models for both runs
source("scripts_RemoteSensing for E and C_/03_h2o_train.R")

# Run SSDM models for both runs
source("scripts_RemoteSensing for E and C_/04_ssdm_train.R")

# Generate comparison plots
source("scripts_RemoteSensing for E and C_/05_h2o_vs_ssdm.R")
source("scripts_RemoteSensing for E and C_/05b_temporal_comparison.R")
source("scripts_RemoteSensing for E and C_/05a_h2o_vs_ssdm_panel.R")

# Generate appendix figures
source("scripts_RemoteSensing for E and C_/07_appendix.R")
```

## Project Structure

```
BTEH/
├── data/                    # Input data
│   ├── clean/              # GPS telemetry data
│   ├── envi/               # Environmental layers
│   └── shp/                # Shapefiles
├── scripts_RemoteSensing for E and C_/  # Main analysis scripts
├── R/                      # Helper functions
├── results/                # Output files
├── logs/                   # Log files
└── config.yml             # Configuration
```

## Scripts

### Main Analysis
- `03_h2o_train.R` - Train H2O AutoML models
- `04_ssdm_train.R` - Train SSDM ensemble models
- `05_h2o_vs_ssdm.R` - Compare methods
- `05a_h2o_vs_ssdm_panel.R` - Before/after panels
- `05b_temporal_comparison.R` - Temporal analysis
- `07_appendix.R` - Supplementary figures

### Configuration
Edit `config.yml` to change:
- Run mode (FAST or REPRO)
- Number of threads
- Memory allocation
- Species to analyze

## Output

After running all scripts:

```
results/
├── H2O/                    # H2O model predictions
│   ├── A/                 # Run A (3 species)
│   └── B/                 # Run B (6 species)
├── SSDM/                   # SSDM ensemble predictions
│   ├── A/
│   └── B/
└── compare/                # Comparison plots and metrics
```

## Troubleshooting

### Package installation fails
```r
# Update renv
install.packages("renv")
renv::restore()
```

### H2O won't start
```r
# Check Java
system("java -version")

# Initialize with less memory
library(h2o)
h2o.init(nthreads = 1, max_mem_size = "4G")
```

### Script can't find config.yml
Make sure you're running from project root:
```r
setwd("C:/path/to/BTEH/BTEH")
file.exists("config.yml")  # Should return TRUE
```

## Citation

If you use this code, please cite:

[Add your citation here]

## License

[Add license information]

## Contact

[Add contact information]
