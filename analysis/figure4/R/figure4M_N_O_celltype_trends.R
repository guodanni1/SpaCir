# Figure 4M/N/O cell-type trend panels with EarlyGC versus LateGC significance marks.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "00_setup.R"))

run_figure4_celltype_trend_panels <- function(sample_id = "Tonsil2") {
  obj_info <- load_sample_objects(sample_id)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$filtered_object %||% obj_info$object
  if (is.null(seurat_obj)) {
    stop("A Seurat object is required for cell-type trend panels.", call. = FALSE)
  }
  if (!"GC_Maturity_Label" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- assign_gc_maturity(
      seurat_obj,
      gc_col = "GC_part",
      assay = "SCT",
      early_markers = early_gc_markers,
      late_markers = late_gc_markers,
      score_nbin = gc_maturity_score_nbin
    )$object
  }

  celltype_data <- prepare_celltype_trend(seurat_obj)
  p_b <- plot_celltype_category_trend(celltype_data, "B cell")
  p_t <- plot_celltype_category_trend(celltype_data, "T cell")
  p_o <- plot_celltype_category_trend(celltype_data, "Other cell")

  save_panel(p_b, paste0("figure4M_", cfg$output_prefix, "_b_cell_trend.pdf"), width = 4.8, height = 3.2)
  save_panel(p_t, paste0("figure4N_", cfg$output_prefix, "_t_cell_trend.pdf"), width = 4.8, height = 3.2)
  save_panel(p_o, paste0("figure4O_", cfg$output_prefix, "_other_cell_trend.pdf"), width = 5.4, height = 3.2)

  invisible(celltype_data)
}
