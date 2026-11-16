# Mass Cytometry Explorations

## Environment and validation

- Create the Conda environment from `environment.yaml`, then activate `mass-cytometry` before rendering:
  ```bash
  conda env create -f environment.yaml
  conda activate mass-cytometry
  ```
- The primary executable units are Quarto notebooks. Render an individual stage (the equivalent of a focused check) from the repository root:
  ```bash
  quarto render notebooks/01_data_acquisition.qmd
  quarto render notebooks/02_quality_control.qmd
  quarto render notebooks/03_preprocessing.qmd
  quarto render notebooks/04_dimensionality_reduction.qmd
  quarto render notebooks/05_clustering_and_annotation.qmd
  ```
- There is no configured automated test suite, linter, formatter, or project-wide build command. Validate notebook changes by rendering the changed stage and any downstream stage whose input/output contract changed.

## Pipeline architecture

The repository is a sequential R/Quarto replication analysis of the `JacksonFischer_2020_BreastCancer` `imcdatasets` dataset. The notebooks are the authoritative workflow:

1. `01_data_acquisition.qmd` loads the Basel cohort, stores images and segmentation masks as on-disk HDF5, validates image linkage, and writes `results/sce_raw.rds`.
2. `02_quality_control.qmd` applies image, cell, and marker QC to that object and writes `results/sce_qc_filtered.rds`.
3. `03_preprocessing.qmd` selects an analysis assay, filters markers, adds a marker-wise `scaled` assay, flags variable markers in `rowData(...$use_for_dr)`, and writes `results/sce_preprocessed.rds`.
4. `04_dimensionality_reduction.qmd` uses those features to calculate PCA plus optional UMAP/t-SNE, stores embeddings in `reducedDims`, and writes `results/sce_dimensionality_reduced.rds`.
5. `05_clustering_and_annotation.qmd` runs seeded k-means on leading PCA coordinates, writes `colData(...$cluster_id)`, and saves `results/sce_clustered_annotated.rds`.

Each stage also writes CSV summaries under its own `results/*_summary/` directory. Treat the saved `.rds` object, its assay names, `rowData`, `colData`, metadata, and reduced-dimension names as the interface between stages; update downstream notebooks when changing any of them.

`R/01_data_acquisition.R` is a root-relative script alternative for Step 1 only. Later stages are implemented in Quarto notebooks, so do not assume similarly named R scripts exist.

## Repository-specific conventions

- Notebooks support being rendered from either the repository root or `notebooks/`, but normalize their working directory to `notebooks/`; paths to shared artifacts therefore use `../results`, `../data`, and `../figures`.
- Keep image work constrained to the Basel subset (`full_dataset = FALSE`, `cohort = "Basel"`). On a 16 GB machine, images and masks must remain HDF5-backed in `data/h5_cache/`; the acquisition notebook batch-loads images to avoid out-of-memory failures.
- The ExperimentHub cache is project-local at `data/ExperimentHub/`. Both it and HDF5 data are runtime caches and are ignored, as are generated `results/` and `figures/` outputs.
- Preserve the notebooks' fail-fast checks for missing upstream `.rds` files and required embeddings rather than adding silent fallbacks.
- Step 5 deliberately uses PCA coordinates for clustering and UMAP only for visualization. Its `cluster_annotation_template.csv` is an editable input: preserve existing annotation labels when regenerating it, and keep `cache: false` so edits are reread on every render.
