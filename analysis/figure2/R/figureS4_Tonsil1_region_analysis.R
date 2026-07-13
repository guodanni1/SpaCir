# Figure S4: supplementary region analysis for Tonsil1.
#
# Panel map:
#   S4A: nCount/nFeature QC before and after filtering.
#   S4B: UMAP/spatial region annotation.
#   S4C: cell2location cell-type composition by region.
#   S4D: marker-gene dotplot for region annotation.
#   S4E-G: IGH/TRB Shannon diversity comparisons.
#   S4F-H: IGH/TRB abundance-versus-clonotype correlation heatmaps.
#
# Shared code:
#   S4B/C/E-G/F-H reuse Figure 2 helper functions. S4A and S4D are specific to
#   the Tonsil1 supplement and keep their own QC/marker plotting parameters here.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))
source(file.path(figure2_dir, "R", "figure2A_region_annotation_umap_spatial.R"))
source(file.path(figure2_dir, "R", "figure2B_celltype_composition_by_region.R"))
source(file.path(figure2_dir, "R", "figure2G_J_shannon_diversity_by_region.R"))
source(file.path(figure2_dir, "R", "figure2I_K_correlation_heatmaps.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

figureS4_sample_id <- "Tonsil1"

figureS4_regions <- c(
  "Follicular zone",
  "Mantle zone",
  "Subepithelial/Mesenchymal",
  "Epithelium-like structures",
  "Inter-Follicular CD6+",
  "Inter-Follicular XBP1+"
)

figureS4_region_marker_genes <- list(
  "Epithelium-like structures" = c("KRT3", "KRT5", "KRT8", "KRT18", "SPIB", "SYK", "JCHAIN", "SDC1"),
  "Inter-Follicular XBP1+" = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "XBP1"),
  "Subepithelial/Mesenchymal" = c("CCL21", "CCL19", "CCL16", "JCHAIN", "LCN1"),
  "Follicular zone" = c("SERPINA9", "RGS13", "CXCR4", "SMOC2"),
  "Mantle zone" = c("HAPLN1", "MADCAM1", "PLAC8", "CD9"),
  "Inter-Follicular CD6+" = c("CXCR2", "FCER2", "MME", "CD6", "ARHGAP15")
)

plot_figureS4_qc <- function(sample_id = figureS4_sample_id, filtered_rds = NULL) {
  cfg <- sample_row(sample_id)
  obj_raw <- readRDS(file.path(repo_root, cfg$seurat_rds))
  obj_filtered <- if (is.null(filtered_rds)) obj_raw else readRDS(file.path(repo_root, filtered_rds))

  # S4A compares common Visium QC metrics before and after filtering. The
  # original code used nCount_Spatial and nFeature_Spatial.
  p1 <- VlnPlot(obj_raw, features = "nCount_Spatial", pt.size = 0.1) + ggtitle("Tonsil1")
  p2 <- VlnPlot(obj_filtered, features = "nCount_Spatial", pt.size = 0.1) + ggtitle("Tonsil1 filtered")
  p3 <- VlnPlot(obj_raw, features = "nFeature_Spatial", pt.size = 0.1) + ggtitle("Tonsil1")
  p4 <- VlnPlot(obj_filtered, features = "nFeature_Spatial", pt.size = 0.1) + ggtitle("Tonsil1 filtered")
  plot <- (p1 + p2) / (p3 + p4)
  save_panel(plot, "figureS4A_Tonsil1_QC_violin.pdf", 6, 5)
  invisible(plot)
}

plot_figureS4_marker_dotplot <- function(sample_id = figureS4_sample_id, marker_genes = figureS4_region_marker_genes) {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  Idents(obj) <- obj[[cfg$region_column]][, 1]
  features <- unique(unlist(marker_genes, use.names = FALSE))
  features <- features[features %in% rownames(obj)]

  p <- DotPlot(obj, features = features, group.by = cfg$region_column) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_panel(p, "figureS4D_Tonsil1_marker_dotplot.pdf", 9, 3)
  invisible(p)
}

plot_figureS4_all <- function(sample_id = figureS4_sample_id) {
  plot_figureS4_qc(sample_id)
  plot_region_annotation(sample_id)
  plot_celltype_composition(sample_id)
  plot_figureS4_marker_dotplot(sample_id)
  plot_shannon_by_region(sample_id, "IGH", regions = figureS4_regions, output_tag = "figureS4E_G")
  plot_shannon_by_region(sample_id, "TRB", regions = figureS4_regions, output_tag = "figureS4E_G")
  plot_correlation_heatmap(sample_id, "IGH", "Shannon", "figureS4F_H")
  plot_correlation_heatmap(sample_id, "TRB", "Shannon", "figureS4F_H")
}

# Example:
# plot_figureS4_all()
