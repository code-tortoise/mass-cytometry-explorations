# =============================================================================
# Step 1: Data Acquisition & Inspection
# =============================================================================
#
# Study: Jackson, Fischer et al. (2020). The single-cell pathology landscape
#        of breast cancer.
#        Nature 578(7796):615-620.
#        https://doi.org/10.1038/s41586-019-1876-x
#
# Dataset: JacksonFischer_2020_BreastCancer (imcdatasets)
#   - Basel cohort subset (100 images / 100 patients) used for image/mask
#     loading (16 GB RAM); in-memory images need ~19 GB
#   - 42-channel IMC panel, FFPE TMA cores
#
# Outputs (written to results/):
#   - sce_raw.rds          : raw SingleCellExperiment (Basel subset)
#   - acquisition_summary/ : plain-text and CSV summaries
# =============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────

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

# Reproducibility: pin ExperimentHub cache location explicitly
# (avoids repeated downloads across sessions)
Sys.setenv(EXPERIMENT_HUB_CACHE = file.path(getwd(), "data", "ExperimentHub"))

# Output directories
dir.create(file.path("results", "acquisition_summary"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path("figures", "01_acquisition"),      recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path("data",    "h5_cache"),            recursive = TRUE,
           showWarnings = FALSE)

cat(rep("=", 70), "\n", sep = "")
cat("STEP 1: Data Acquisition & Inspection\n")
cat(rep("=", 70), "\n\n", sep = "")


# ── 1a. Load single-cell data (Basel subset: 100 images / 100 patients) ───────
#
# Memory: ~513 MB in RAM — safe on 16 GB
# full_dataset = FALSE + cohort = "Basel"  →  100-image Basel subset
# full_dataset = TRUE                      →  full Basel + Zurich cohorts

cat("Loading SingleCellExperiment (Basel subset) ...\n")

sce <- JacksonFischer_2020_BreastCancer(
    data_type    = "sce",
    full_dataset = FALSE,
    cohort       = "Basel"
)

cat("  Done.\n\n")


# ── 1b. Load cell segmentation masks (on-disk HDF5) ───────────────────────────
#
# Memory without on_disk: ~269 MB — could load in memory, but HDF5 is safer
# HDF5 files are written once and reused across sessions

h5_dir <- file.path("data", "h5_cache")

cat("Loading segmentation masks (on-disk HDF5) ...\n")

masks <- JacksonFischer_2020_BreastCancer(
    data_type    = "masks",
    full_dataset = FALSE,
    cohort       = "Basel",
    on_disk      = TRUE,
    h5FilesPath  = h5_dir
)

cat("  Done.\n\n")


# ── 1c. Load multichannel images (on-disk HDF5 — required on 16 GB) ───────────
#
# WARNING: Images require ~19 GB in RAM (Basel, in-memory).
# on_disk = TRUE writes HDF5 files to h5_dir so only small chunks are read
# during downstream operations.  This is mandatory on 16 GB hardware.
#
# Note: full_dataset = TRUE is not supported for data_type = "images".

cat("Loading multichannel images (on-disk HDF5) ...\n")
cat("  [This may take 10-30 min on first run — HDF5 files are written to",
    h5_dir, "]\n")

images <- JacksonFischer_2020_BreastCancer(
    data_type   = "images",
    cohort      = "Basel",
    on_disk     = TRUE,
    h5FilesPath = h5_dir
)

cat("  Done.\n\n")


# ── 1d. Validate object linkage ───────────────────────────────────────────────
#
# All three objects are linked by image_name. Mismatches indicate a version
# mismatch or download error.

cat("Validating linkage between SCE, images, and masks ...\n")

sce_image_names    <- unique(sce$image_name)
images_image_names <- mcols(images)$image_name
masks_image_names  <- mcols(masks)$image_name

ok_img  <- all(images_image_names %in% sce_image_names)
ok_mask <- all(masks_image_names  %in% sce_image_names)

cat(sprintf("  Images linked to SCE : %s (%d / %d)\n",
            ifelse(ok_img,  "OK", "FAIL"),
            sum(images_image_names %in% sce_image_names),
            length(images_image_names)))

cat(sprintf("  Masks linked to SCE  : %s (%d / %d)\n\n",
            ifelse(ok_mask, "OK", "FAIL"),
            sum(masks_image_names  %in% sce_image_names),
            length(masks_image_names)))

if (!ok_img || !ok_mask) {
    warning("Linkage validation FAILED — check dataset version or re-download.")
}


# ── 1e. Structural inspection ─────────────────────────────────────────────────

cat(rep("-", 60), "\n", sep = "")
cat("SCE object overview\n")
cat(rep("-", 60), "\n", sep = "")
print(sce)
cat("\n")

cat("Assay names:\n")
cat(" ", paste(assayNames(sce), collapse = ", "), "\n\n")

cat("colData columns (cell-level metadata):\n")
cat(" ", paste(names(colData(sce)), collapse = "\n  "), "\n\n")

cat("rowData columns (marker-level metadata):\n")
print(as.data.frame(rowData(sce)))
cat("\n")

cat("Images object:\n")
print(images)
cat("\n")

cat("Masks object:\n")
print(masks)
cat("\n")


# ── 1f. EDA: Dataset-level summaries ──────────────────────────────────────────

cat(rep("-", 60), "\n", sep = "")
cat("EDA: Dataset-level summaries\n")
cat(rep("-", 60), "\n", sep = "")

# --- Cells per image ---
cells_per_image <- as.data.frame(colData(sce)) |>
    dplyr::count(image_name, name = "n_cells")

cat("\nCells per image (summary across", nrow(cells_per_image), "images):\n")
print(summary(cells_per_image$n_cells))

cat(sprintf("Unique images in subset: %d\n", length(unique(sce$image_name))))

# --- Tumour grades ---
cat("\nCells per tumour grade:\n")
print(table(sce$tumor_grade))

# --- Cells per image ---
cells_per_image_grade <- as.data.frame(colData(sce)) |>
    dplyr::count(image_name, tumor_grade, name = "n_cells")
cat("\nCells per image (summary):\n")
print(summary(cells_per_image_grade$n_cells))

# --- Marker panel ---
panel <- as.data.frame(rowData(sce))
cat(sprintf("\nPanel: %d markers\n", nrow(panel)))

# Save panel to CSV
write.csv(panel,
          file = file.path("results", "acquisition_summary", "marker_panel.csv"),
          row.names = FALSE)
cat("  Marker panel saved to results/acquisition_summary/marker_panel.csv\n")

# Save cells-per-image table
write.csv(cells_per_image,
          file = file.path("results", "acquisition_summary", "cells_per_image.csv"),
          row.names = FALSE)
cat("  Cells-per-image saved to results/acquisition_summary/cells_per_image.csv\n\n")


# ── 1g. EDA: Figures ──────────────────────────────────────────────────────────

cat(rep("-", 60), "\n", sep = "")
cat("EDA: Generating figures\n")
cat(rep("-", 60), "\n", sep = "")

# --- Figure 1: Cells per image (sorted bar chart) ---
p_cells_per_image <- ggplot(cells_per_image,
                             aes(x = reorder(image_name, n_cells),
                                 y = n_cells)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_hline(yintercept = 100, colour = "red", linetype = "dashed",
               linewidth = 0.6) +
    annotate("text", x = 1, y = 130, label = "QC threshold (100 cells)",
             hjust = 0, colour = "red", size = 2.8) +
    coord_flip() +
    labs(
        title    = "Step 1 — Cells per image",
        subtitle = sprintf("JacksonFischer_2020_BreastCancer (Basel)  |  n = %d images",
                           nrow(cells_per_image)),
        x        = "Image",
        y        = "Cell count"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 5))

ggsave(
    filename = file.path("figures", "01_acquisition", "cells_per_image.pdf"),
    plot     = p_cells_per_image,
    width    = 8, height = max(6, nrow(cells_per_image) * 0.12)
)
cat("  Saved: figures/01_acquisition/cells_per_image.pdf\n")

# --- Figure 2: Cells per image by tumour grade (violin + jitter) ---
p_group <- ggplot(cells_per_image_grade,
                  aes(x = tumor_grade, y = n_cells,
                      fill = tumor_grade)) +
    geom_violin(alpha = 0.6, colour = "grey40") +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7, colour = "grey20") +
    labs(
        title    = "Step 1 — Total cells per image by tumour grade",
        subtitle = "Each point = one image",
        x        = "Tumour grade",
        y        = "Total cell count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 30, hjust = 1))

ggsave(
    filename = file.path("figures", "01_acquisition", "cells_per_patient_group.pdf"),
    plot     = p_group,
    width    = 7, height = 5
)
cat("  Saved: figures/01_acquisition/cells_per_patient_group.pdf\n")

# --- Figure 3: Panel dot plot (metal tag × target) ---
if (all(c("channel_name", "marker_name") %in% colnames(panel))) {
    panel_plot <- panel
    panel_plot$marker_index <- seq_len(nrow(panel_plot))

    p_panel <- ggplot(panel_plot,
                      aes(x = 1, y = reorder(marker_name, marker_index),
                          label = channel_name)) +
        geom_point(colour = "steelblue", size = 3) +
        geom_text(hjust = -0.3, size = 2.8) +
        labs(
            title    = "Step 1 — Antibody panel",
            subtitle = sprintf("%d channels", nrow(panel_plot)),
            x        = NULL,
            y        = "Marker"
        ) +
        xlim(0.9, 1.5) +
        theme_bw(base_size = 10) +
        theme(axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.x = element_blank())

    ggsave(
        filename = file.path("figures", "01_acquisition", "antibody_panel.pdf"),
        plot     = p_panel,
        width    = 5, height = max(5, nrow(panel_plot) * 0.3)
    )
    cat("  Saved: figures/01_acquisition/antibody_panel.pdf\n")
}

cat("\n")


# ── 1h. Save raw SCE object ───────────────────────────────────────────────────

cat(rep("-", 60), "\n", sep = "")
cat("Saving raw SCE object\n")
cat(rep("-", 60), "\n", sep = "")

out_sce <- file.path("results", "sce_raw.rds")
saveRDS(sce, file = out_sce)
cat(sprintf("  Saved: %s  (%d cells x %d markers)\n\n", out_sce,
            ncol(sce), nrow(sce)))


# ── 1i. Session info ──────────────────────────────────────────────────────────

sink(file.path("results", "acquisition_summary", "session_info.txt"))
sessionInfo()
sink()

cat(rep("=", 70), "\n", sep = "")
cat("STEP 1 COMPLETE\n")
cat("  Output SCE  : results/sce_raw.rds\n")
cat("  Figures     : figures/01_acquisition/\n")
cat("  Summaries   : results/acquisition_summary/\n")
cat("  HDF5 images : data/h5_cache/\n")
cat(rep("=", 70), "\n", sep = "")
