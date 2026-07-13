# Figure 2C: marker-gene violin plots supporting region annotation.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

region_marker_genes <- list(
  "Follicular DZ" = c("AICDA", "MKI67"),
  "Inter LZ-DZ" = c("GKV3-11", "PDCD1"),
  "Follicular LZ" = c("LYAR", "FCRL1"),
  "Mantle zone" = c("CD7", "CD5", "LAPTM5"),
  "Inter-Follicular CD5+" = c("LAG3", "HBA2"),
  "Vascular/Mesenchymal" = c("LUM", "DCN"),
  "Vascular" = c("DST", "SPARCL1", "CXCL1", "TNFRSF6B")
)

plot_region_markers <- function(sample_id, assay = "Spatial", region_col = NULL, marker_genes = region_marker_genes) {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  region_col <- region_col %||% cfg$region_column
  DefaultAssay(obj) <- assay
  Idents(obj) <- obj[[region_col]][, 1]

  features <- unique(unlist(marker_genes, use.names = FALSE))
  features <- features[features %in% rownames(obj)]
  if (length(features) == 0) {
    stop("None of the marker genes were found in the Seurat object.", call. = FALSE)
  }

  # The Word code uses VlnPlot with region identities; keeping the same display
  # makes it clear which marker genes drove each anatomical annotation.
  p <- VlnPlot(obj, features = features, group.by = region_col, pt.size = 0, stack = TRUE, flip = TRUE) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  save_panel(p, paste0("figure2C_", display_sample(cfg$sample_id), "_region_marker_violin.pdf"), 9, 4.5)
  invisible(p)
}

# Example:
# plot_region_markers("Tonsil2")
