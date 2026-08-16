#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
suppressPackageStartupMessages({
    library(imcdatasets)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(cytomapper)
    library(HDF5Array)
    library(ExperimentHub)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})
#
#
#
#
#
# Pin ExperimentHub cache to the project directory so data is only
# downloaded once and reused across sessions.
Sys.setenv(EXPERIMENT_HUB_CACHE = file.path(getwd(), "..", "data", "ExperimentHub"))

# HDF5 directory for on-disk image storage
h5_dir <- file.path("..", "data", "h5_cache")

# Output directories
dir.create(file.path("..", "results", "acquisition_summary"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("..", "figures", "01_acquisition"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(h5_dir, recursive = TRUE, showWarnings = FALSE)
#
#
#
#
#
#
#
#
#
#
#
#
#
sce <- Meyer_2025_TripleNegativeBreastCancer(
    data_type    = "sce",
    full_dataset = FALSE   # 125 images / 60 patients
)
sce
#
#
#
#
#
#
#
#
#
masks <- Meyer_2025_TripleNegativeBreastCancer(
    data_type    = "masks",
    full_dataset = FALSE,
    on_disk      = TRUE,
    h5FilesPath  = h5_dir
)
masks
#
#
#
#
#
#
#
#
#
#
#
#
images <- Meyer_2025_TripleNegativeBreastCancer(
    data_type   = "images",
    on_disk     = TRUE,
    h5FilesPath = h5_dir   # REQUIRED on 16 GB — in-memory needs ~21 GB
)
images
#
#
#
#
#
#
#
#
#
#
sce_names    <- unique(sce$image_name)
images_names <- mcols(images)$image_name
masks_names  <- mcols(masks)$image_name

ok_img  <- all(images_names %in% sce_names)
ok_mask <- all(masks_names  %in% sce_names)

cat(sprintf("Images → SCE : %s  (%d / %d matched)\n",
            ifelse(ok_img,  "✓ OK", "✗ FAIL"),
            sum(images_names %in% sce_names), length(images_names)))

cat(sprintf("Masks  → SCE : %s  (%d / %d matched)\n",
            ifelse(ok_mask, "✓ OK", "✗ FAIL"),
            sum(masks_names  %in% sce_names), length(masks_names)))

if (!ok_img || !ok_mask) {
    stop("Linkage validation FAILED — check dataset version or re-download.")
}
#
#
#
#
#
#
#
#
#
#
#
assayNames(sce)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
names(colData(sce))
#
#
#
#
#
#
#
as.data.frame(rowData(sce))
#
#
#
#
#
#
#
#
#
#
#
cells_per_image <- as.data.frame(colData(sce)) |>
    dplyr::count(image_name, name = "n_cells")

summary(cells_per_image$n_cells)
#
#
#
#| label: fig-cells-per-image
#| fig-cap: "Cells per image, sorted ascending. The dashed red line marks a typical QC threshold of 100 cells — images below this may represent failed or poor-quality acquisitions."

ggplot(cells_per_image,
       aes(x = reorder(image_name, n_cells), y = n_cells)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_hline(yintercept = 100, colour = "red",
               linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = 3, y = 130,
             label = "QC threshold (100 cells)",
             hjust = 0, colour = "red", size = 3) +
    coord_flip() +
    labs(
        title    = "Cells per image",
        subtitle = sprintf("n = %d images", nrow(cells_per_image)),
        x = NULL, y = "Cell count"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 6))
#
#
#
#
#
#
#
cat("Cells per patient group:\n")
print(table(sce$patient_patientgroup))

cat(sprintf("\nUnique patients: %d\n", length(unique(sce$patient_id))))
#
#
#
#| label: fig-patient-groups
#| fig-cap: "Total cell count per patient, broken down by patient group. Each point represents one patient. Violin shape shows the distribution across patients within each group."

cells_per_patient <- as.data.frame(colData(sce)) |>
    dplyr::count(patient_id, patient_patientgroup, name = "n_cells")

ggplot(cells_per_patient,
       aes(x = patient_patientgroup, y = n_cells,
           fill = patient_patientgroup)) +
    geom_violin(alpha = 0.6, colour = "grey40") +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, colour = "grey20") +
    labs(
        title    = "Total cells per patient by group",
        subtitle = "Each point = one patient",
        x = "Patient group", y = "Total cell count"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 30, hjust = 1))
#
#
#
#
#
#| label: fig-panel
#| fig-cap: "39-channel IMC antibody panel. Each point is one marker; the metal tag (channel name) is shown to the right."

panel <- as.data.frame(rowData(sce))

if (all(c("channel_name", "marker_name") %in% colnames(panel))) {
    panel$marker_index <- seq_len(nrow(panel))

    ggplot(panel,
           aes(x = 1,
               y = reorder(marker_name, marker_index),
               label = channel_name)) +
        geom_point(colour = "steelblue", size = 3) +
        geom_text(hjust = -0.25, size = 3) +
        xlim(0.9, 1.5) +
        labs(
            title    = "Antibody panel",
            subtitle = sprintf("%d channels", nrow(panel)),
            x = NULL, y = "Marker"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.x = element_blank())
}
#
#
#
#
#
#
#
# Marker panel → CSV
write.csv(panel,
          file = file.path("..", "results", "acquisition_summary",
                           "marker_panel.csv"),
          row.names = FALSE)

# Cells per image → CSV
write.csv(cells_per_image,
          file = file.path("..", "results", "acquisition_summary",
                           "cells_per_image.csv"),
          row.names = FALSE)

# Raw SCE → RDS (used as input for Step 2: QC)
saveRDS(sce,
        file = file.path("..", "results", "sce_raw.rds"))

cat("Saved:\n")
cat("  results/acquisition_summary/marker_panel.csv\n")
cat("  results/acquisition_summary/cells_per_image.csv\n")
cat(sprintf("  results/sce_raw.rds  (%d cells × %d markers)\n",
            ncol(sce), nrow(sce)))
#
#
#
#
#
#
#
sessionInfo()
#
#
#
#
#
#
#
#
#
#
#
#
#
