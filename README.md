# Mass Cytometry Explorations

Replication analysis of the Meyer et al. (2025) study:

> Meyer, Jackson et al. (2025). *A stratification system for breast cancer based on basoluminal tumor cells and spatial tumor architecture.* Cancer Cell 43(9):1637–1655.e9. https://doi.org/10.1016/j.ccell.2025.06.019

Data are loaded from the `imcdatasets` Bioconductor package using the `Meyer_2025_TripleNegativeBreastCancer` dataset.

## Hardware requirements

These scripts are configured for a **16 GB RAM** environment (e.g. WSL2). Images are stored on disk via HDF5 to avoid loading ~21 GB into memory. The **subset dataset** (`full_dataset = FALSE`, 125 images / 60 patients) is used for all image-level and spatial analyses; the full SCE (`full_dataset = TRUE`, 450 images / 215 patients) is used for single-cell-level steps only.

Recommended `.wslconfig` settings:
```ini
[wsl2]
memory=14GB
swap=8GB
processors=8
```

## Project structure

```
notebooks/
  01_data_acquisition.qmd   # Step 1: load data, inspect structure, validate linkage (Quarto)
  02_quality_control.qmd    # Step 2: cell/image/marker QC (Quarto)
  03_preprocessing.qmd      # Step 3: transformation, normalisation, EDA (Quarto)
R/
  01_data_acquisition.R     # Step 1: plain R script equivalent
data/
  h5_cache/                 # HDF5-backed images (auto-created, git-ignored)
results/                    # Saved .rds / .csv outputs
figures/                    # Saved plots
```

## Dependencies

Install via `BiocManager`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "imcdatasets",
    "SingleCellExperiment",
    "SpatialExperiment",
    "cytomapper",
    "HDF5Array",
    "ExperimentHub",
    "scater",
    "scuttle",
    "ggplot2",
    "patchwork",
    "dplyr",
    "reshape2",
    "tidyr",
    "ggridges"
))
```

## Running the pipeline

### Option A — Quarto notebooks (recommended, human-readable)

Open any `.qmd` file in RStudio or VS Code and click **Render**, or from the terminal:

```bash
quarto render notebooks/01_data_acquisition.qmd
```

This produces a self-contained `01_data_acquisition.html` file with narrative text,
code, and all figures inline — readable in any browser without R installed.

### Option B — Plain R scripts

Run each script in order:

```r
source("R/01_data_acquisition.R")
source("R/02_quality_control.R")
source("R/03_preprocessing.R")
```

Intermediate objects are saved to `results/` so each step can be re-run independently.

## Citation

Please cite the original study and the `imcdatasets` package:

- Meyer, Jackson et al. (2025). Cancer Cell 43(9):1637–1655.e9.
- Damond N, Steenbuck N, Eling N, Fischer J, Hoch T, Meyer L (2026). *imcdatasets: Collection of publicly available imaging mass cytometry (IMC) datasets.*
