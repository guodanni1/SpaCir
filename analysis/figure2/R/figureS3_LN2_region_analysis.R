# Figure S3: supplementary region analysis for LN2.
#
# Panel map:
#   S3A: UMAP/spatial region annotation.
#   S3B: cell2location cell-type composition by region.
#   S3C: marker-gene violin plot for region annotation.
#   S3D-F: region-specific cell-type abundance violin plots.
#   S3G-H: IGH/TRB Shannon diversity and IGH SHM comparisons.
#   S3I: IGH/TRB abundance-versus-clonotype correlation heatmaps.
#   S3J-K: random B/T cell-type clonotype diversity and B-T correlation.
#
# Shared code:
#   Region annotation, composition, marker violin, diversity, SHM and correlation
#   panels reuse the Figure 2 functions. S3-specific parameters are declared
#   below so the supplement remains explicit.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))
source(file.path(figure2_dir, "R", "figure2A_region_annotation_umap_spatial.R"))
source(file.path(figure2_dir, "R", "figure2B_celltype_composition_by_region.R"))
source(file.path(figure2_dir, "R", "figure2C_marker_violin.R"))
source(file.path(figure2_dir, "R", "figure2D_E_F_celltype_abundance_violin.R"))
source(file.path(figure2_dir, "R", "figure2G_J_shannon_diversity_by_region.R"))
source(file.path(figure2_dir, "R", "figure2H_shm_rate_by_region.R"))
source(file.path(figure2_dir, "R", "figure2I_K_correlation_heatmaps.R"))
source(file.path(figure2_dir, "R", "random_celltype_diversity_helpers.R"))

figureS3_sample_id <- "LN2"

figureS3_regions <- c(
  "Follicular DZ",
  "Follicular LZ",
  "Mantle zone",
  "Inter-Follicular SDC1+",
  "Inter-Follicular low SDC1+",
  "Inter-Follicular IL7R+",
  "Inter-Follicular CXCL10+",
  "Epithelial/Muscle",
  "Subepithelial/Mesenchymal/Vascular"
)

figureS3_region_marker_genes <- list(
  "Follicular DZ" = c("AICDA", "SUGCT"),
  "Follicular LZ" = c("PDCD1", "COL15A1"),
  "Epithelial/Muscle" = c("TNFRSF6B", "SPARCL1", "SELD", "ACKR1"),
  "Inter-Follicular CXCL10+" = c("CXCL10", "CHI3L1", "HBA1"),
  "Mantle zone" = c("CXCL13", "CLEC4G"),
  "Inter-Follicular IL7R+" = c("IL7R"),
  "Inter-Follicular SDC1+" = c("SDC1")
)

figureS3_abundance_panels <- list(
  FigureS3D = list(
    regions = c("Follicular DZ", "Follicular LZ", "Mantle zone"),
    cell_types = c("B Cycling", "B GC DZ", "B GC LZ", "B GC prePB", "FDC", "T CD4+ TfH GC", "B naive", "B mem")
  ),
  FigureS3E = list(
    regions = c("Inter-Follicular SDC1+", "Inter-Follicular low SDC1+", "Inter-Follicular IL7R+", "Inter-Follicular CXCL10+"),
    cell_types = c("B plasma", "DC pDC", "DC cDC1", "DC cDC2", "DC CCR7+", "T CD8+ cytotoxic", "Macrophages M2", "Macrophages M1", "T CD4+", "T CD8+ naive", "B activated", "Mast", "B GC prePB")
  ),
  FigureS3F = list(
    regions = c("Epithelial/Muscle", "Subepithelial/Mesenchymal/Vascular"),
    cell_types = c("Endo", "DC CCR7+", "T CD8+ naive", "Macrophages M1", "Mast", "B plasma", "B naive", "Macrophages M2")
  )
)

plot_figureS3_abundance_panels <- function(sample_id = figureS3_sample_id) {
  for (panel in names(figureS3_abundance_panels)) {
    params <- figureS3_abundance_panels[[panel]]
    plot_abundance_violin(sample_id, params$regions, params$cell_types, output_tag = panel)
  }
}

plot_figureS3_all <- function(sample_id = figureS3_sample_id) {
  plot_region_annotation(sample_id)
  plot_celltype_composition(sample_id)
  plot_region_markers(sample_id, marker_genes = figureS3_region_marker_genes)
  plot_figureS3_abundance_panels(sample_id)
  plot_shannon_by_region(sample_id, "IGH", regions = figureS3_regions, output_tag = "figureS3G")
  plot_shannon_by_region(sample_id, "TRB", regions = figureS3_regions, output_tag = "figureS3G")
  plot_shm_by_region(sample_id, regions = figureS3_regions)
  plot_correlation_heatmap(sample_id, "IGH", "Shannon", "figureS3I")
  plot_correlation_heatmap(sample_id, "TRB", "Shannon", "figureS3I")
}

# Example:
# plot_figureS3_all()
