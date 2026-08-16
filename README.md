# Mass Cytometry Explorations

Replication analysis of the Jackson, Fischer et al. (2020) study:

> Jackson, Fischer et al. (2020). *The single-cell pathology landscape of breast cancer.* Nature 578(7796):615–620. https://doi.org/10.1038/s41586-019-1876-x

Data are loaded from the `imcdatasets` Bioconductor package using the `JacksonFischer_2020_BreastCancer` dataset.

## Hardware requirements

These scripts are configured for a **16 GB RAM** environment (e.g. WSL2). Images are stored on disk via HDF5 to avoid loading ~19 GB into memory. The **Basel subset** (`full_dataset = FALSE`, `cohort = "Basel"`, 100 images / 100 patients) is used for all image-level and spatial analyses; the full SCE (`full_dataset = TRUE`, Basel + Zurich cohorts) is used for single-cell-level steps only.

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

- Jackson, Fischer et al. (2020). Nature 578(7796):615–620.
- Damond N, Steenbuck N, Eling N, Fischer J, Hoch T, Meyer L (2026). *imcdatasets: Collection of publicly available imaging mass cytometry (IMC) datasets.*
