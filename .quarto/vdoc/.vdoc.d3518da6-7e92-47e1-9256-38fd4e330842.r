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
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Resolve project/notebook paths robustly, without absolute paths.
if (basename(getwd()) == "notebooks") { # Check if running from notebooks/ directory.
  notebooks_dir <- normalizePath(getwd())
} else if (dir.exists(file.path(getwd(), "notebooks"))) { # Check if running from project root.
  notebooks_dir <- normalizePath(file.path(getwd(), "notebooks"))
} else { # If neither then stop with an error message.
  stop("Run from project root or notebooks/ directory.")
}

if (!identical(normalizePath(getwd()), notebooks_dir)) { # Set working directory to notebooks/ if not already there.
  setwd(notebooks_dir)
} 

cat("Working directory:", getwd(), "\n") # Tell the user where the notebook is running from.
#
#
#
dir.create(file.path("..", "results", "preprocessing_summary"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("..", "figures", "03_preprocessing"),
           recursive = TRUE, showWarnings = FALSE)
#
#
#
#
#
#
#
sce_path <- file.path("..", "results", "sce_qc_filtered.rds")
if (!file.exists(sce_path)) {
  stop("Missing results/sce_qc_filtered.rds. Please render notebooks/02_quality_control.qmd first.")
}

sce_qc <- readRDS(sce_path)
sce_qc
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
available_assays <- assayNames(sce_qc)
preferred_assays <- c("exprs", "quant_norm", "counts")
analysis_assay <- preferred_assays[preferred_assays %in% available_assays][1]

if (is.na(analysis_assay)) {
  stop("No supported assay found. Expected one of: exprs, quant_norm, counts.")
}

analysis_mat <- assay(sce_qc, analysis_assay)
cat("Selected assay:", analysis_assay, "\n")
#
#
#
#
#
#
#
#
min_pct_nonzero <- 1
n_top_variable_markers <- 30L

marker_metrics <- data.frame(
  marker_name = rownames(sce_qc),
  mean_expr = rowMeans(analysis_mat, na.rm = TRUE),
  median_expr = apply(analysis_mat, 1, median, na.rm = TRUE),
  variance_expr = apply(analysis_mat, 1, var, na.rm = TRUE),
  pct_cells_nonzero = rowMeans(analysis_mat > 0, na.rm = TRUE) * 100
) |>
  dplyr::mutate(
    keep_marker = pct_cells_nonzero >= min_pct_nonzero & is.finite(variance_expr) & variance_expr > 0
  ) |>
  dplyr::arrange(dplyr::desc(variance_expr))

n_keep <- sum(marker_metrics$keep_marker)
cat(sprintf("Markers before preprocessing: %d\n", nrow(marker_metrics)))
cat(sprintf("Markers retained: %d\n", n_keep))
#
#
#
#| fig-cap: "Marker variance versus prevalence (% cells with non-zero signal). Red points fail marker preprocessing."
base_plot <- ggplot(marker_metrics,
       aes(x = pct_cells_nonzero, y = variance_expr, colour = keep_marker)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick")) +
  labs(
    title = "Marker preprocessing diagnostics",
    subtitle = sprintf("Keep threshold: >= %.1f%% non-zero cells", min_pct_nonzero),
    x = "% cells with non-zero signal",
    y = "Variance",
    colour = "Keep marker"
  ) +
  theme_bw()
base_plot


#
#
#
#
#
#
#
#
keep_markers <- marker_metrics$marker_name[marker_metrics$keep_marker] # Keep only markers that are TRUE for keep_marker obj.
sce_preprocessed <- sce_qc[keep_markers, ] # In sce_qc, keep only the selected markers in keep_markers, retaining all cells.

# Select assay for downstream analysis, as chosen in chunk r assay-selection.
analysis_mat_keep <- assay(sce_preprocessed, analysis_assay)

# Add scaled assay for downstream methods that assume marker-wise centering/scaling, using chosen assay.
scaled_mat <- t(scale(t(analysis_mat_keep))) # Transpose, scale, and transpose back to maintain marker-wise scaling.
scaled_mat[!is.finite(scaled_mat)] <- 0 # Scale introduces NaN/Inf for constant markers, replace with 0.
assay(sce_preprocessed, "scaled") <- scaled_mat # Add the scaled matrix as a new assay in the preprocessed SingleCellExperiment object.

# Flag top variable markers (for DR/clustering feature selection).
marker_metrics_keep <- marker_metrics |>
  dplyr::filter(keep_marker)

top_n <- min(n_top_variable_markers, nrow(marker_metrics_keep))
top_markers <- marker_metrics_keep$marker_name[seq_len(top_n)]

rowData(sce_preprocessed)$marker_variance <- marker_metrics_keep$variance_expr[
  match(rownames(sce_preprocessed), marker_metrics_keep$marker_name)
]
rowData(sce_preprocessed)$pct_cells_nonzero <- marker_metrics_keep$pct_cells_nonzero[
  match(rownames(sce_preprocessed), marker_metrics_keep$marker_name)
]
rowData(sce_preprocessed)$use_for_dr <- rownames(sce_preprocessed) %in% top_markers

metadata(sce_preprocessed)$preprocessing <- list(
  analysis_assay = analysis_assay,
  min_pct_nonzero = min_pct_nonzero,
  n_top_variable_markers = n_top_variable_markers,
  n_markers_retained = nrow(sce_preprocessed),
  n_cells_retained = ncol(sce_preprocessed)
)

cat(sprintf("Cells retained: %d\n", ncol(sce_preprocessed)))
cat(sprintf("Markers retained: %d\n", nrow(sce_preprocessed)))
cat(sprintf("Top variable markers flagged for DR: %d\n", sum(rowData(sce_preprocessed)$use_for_dr)))
#
#
#
#| fig-cap: "Top variable markers retained for dimensionality reduction."
top_markers_tbl <- marker_metrics_keep |>
  dplyr::slice_head(n = top_n) |>
  dplyr::mutate(
    marker_name = factor(marker_name, levels = rev(marker_name))
  )

ggplot(top_markers_tbl, aes(x = marker_name, y = variance_expr)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top variable markers",
    subtitle = sprintf("Top %d by variance in '%s' assay", top_n, analysis_assay),
    x = NULL,
    y = "Variance"
  ) +
  theme_bw()
#
#
#
#
#
#
#
params_tbl <- data.frame(
  parameter = c("analysis_assay", "min_pct_nonzero", "n_top_variable_markers",
                "n_cells_retained", "n_markers_retained"),
  value = c(
    analysis_assay,
    as.character(min_pct_nonzero),
    as.character(n_top_variable_markers),
    as.character(ncol(sce_preprocessed)),
    as.character(nrow(sce_preprocessed))
  )
)

write.csv(marker_metrics,
          file = file.path("..", "results", "preprocessing_summary", "marker_preprocessing_metrics.csv"),
          row.names = FALSE)

write.csv(params_tbl,
          file = file.path("..", "results", "preprocessing_summary", "preprocessing_parameters.csv"),
          row.names = FALSE)

saveRDS(sce_preprocessed, file = file.path("..", "results", "sce_preprocessed.rds"))

params_tbl 

cat("Saved:\n")
cat("  results/preprocessing_summary/marker_preprocessing_metrics.csv\n")
cat("  results/preprocessing_summary/preprocessing_parameters.csv\n")
cat("  results/sce_preprocessed.rds\n")
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
