# Figure 2A: region annotation on UMAP and spatial coordinates.
# The same script can be used for Figure S3A and Figure S4B by changing sample_id.
#
# Per-sample clustering parameters are recorded in:
#   config/sample_specific_parameters.csv
# Region marker genes used to name the clusters are recorded in:
#   config/region_marker_sets.csv

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

plot_region_annotation <- function(sample_id, reduction = "umap", region_col = NULL) {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  region_col <- region_col %||% cfg$region_column
  if (!region_col %in% colnames(obj@meta.data)) {
    stop("Region column not found in Seurat metadata: ", region_col, call. = FALSE)
  }

  obj[[region_col]][, 1] <- factor(obj[[region_col]][, 1])
  present_regions <- levels(obj[[region_col]][, 1])
  colors <- region_colors[present_regions]
  colors[is.na(colors)] <- "#BDBDBD"

  # UMAP shows transcriptome-defined region separation; SpatialDimPlot checks
  # that the same labels map to anatomically coherent tissue areas.
  p_umap <- DimPlot(obj, reduction = reduction, group.by = region_col, cols = colors) +
    ggtitle(display_sample(cfg$sample_id)) +
    theme_void() +
    theme(legend.position = "right")
  p_spatial <- SpatialDimPlot(obj, group.by = region_col, cols = colors, pt.size.factor = 1.6) +
    ggtitle(display_sample(cfg$sample_id)) +
    theme(legend.position = "right")
  plot <- p_umap + p_spatial + plot_layout(widths = c(1, 1.25))
  save_panel(plot, paste0("figure2A_", display_sample(cfg$sample_id), "_region_annotation.pdf"), 10, 4)
  invisible(plot)
}

# Example:
# plot_region_annotation("Tonsil2")
