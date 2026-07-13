# Figure 4E/F and Figure S8D/E WGCNA and GO panels.
#
# WGCNA and GO settings are kept explicit in config/hdwgcna_go_parameters.csv.
# The final public panels use EarlyGC and LateGC only and the curated GO terms
# selected for the manuscript figures.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "00_setup.R"))

run_figure4_wgcna_go_panels <- function(sample_id = "Tonsil2") {
  obj_info <- load_sample_objects(sample_id)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$filtered_object %||% obj_info$object

  if (!is.null(seurat_obj) && !"GC_Maturity_Label" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- assign_gc_maturity(
      seurat_obj,
      gc_col = "GC_part",
      assay = "SCT",
      early_markers = early_gc_markers,
      late_markers = late_gc_markers,
      score_nbin = gc_maturity_score_nbin
    )$object
  }

  if (!is.null(seurat_obj)) {
    p_umap <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap",
      group.by = "GC_Maturity_Label",
      cols = maturity_colors
    ) +
      ggtitle(NULL) +
      theme_classic()
    save_panel(p_umap, paste0("figure4F_", cfg$output_prefix, "_maturity_umap.pdf"), width = 4.2, height = 3.8)
  }

  go_table <- read_csv_if_exists(cfg$go_terms_csv, stringsAsFactors = FALSE)
  p_go <- plot_go_dotplot(go_table, group_col = "group")
  save_panel(p_go, paste0("figure4E_", cfg$output_prefix, "_go_maturity.pdf"), width = 6.2, height = 4.6)

  wgcna_go <- read_csv_if_exists(cfg$wgcna_go_terms_csv, stringsAsFactors = FALSE)
  p_wgcna <- plot_wgcna_hub_go_dotplot(wgcna_go)
  save_panel(p_wgcna, paste0("figureS8E_", cfg$output_prefix, "_wgcna_hub_go.pdf"), width = 8, height = 7)

  invisible(list(go = go_table, wgcna_go = wgcna_go))
}
