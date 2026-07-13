# Figure 4A/D and GC-part composition panels.
#
# GC maturity definition used in these panels:
#   EarlyGC markers: CD83, CD86, CXCR5, GPR183, SLAMF1, BCL6
#   LateGC markers: IRF4, PRDM1, XBP1, ZBTB20, FOXP3, DUSP2, IRF8,
#                   GADD45B, JCHAIN, TOX2, SDC1
#   Module scoring: Seurat::AddModuleScore(..., nbin = 10) on the SCT assay.
#   GC label rule: each GC part is EarlyGC if its mean EarlyGC score is
#                  greater than or equal to its mean LateGC score; otherwise
#                  it is LateGC.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "00_setup.R"))

run_figure4_gc_annotation_panels <- function(sample_id = "Tonsil2") {
  obj_info <- load_sample_objects(sample_id)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  if (is.null(seurat_obj)) {
    stop("Seurat object is required for GC annotation panels.", call. = FALSE)
  }

  scored <- assign_gc_maturity(
    seurat_obj,
    gc_col = "GC_part",
    assay = "SCT",
    early_markers = early_gc_markers,
    late_markers = late_gc_markers,
    score_nbin = gc_maturity_score_nbin
  )
  filtered <- scored$object
  metadata_df <- filtered@meta.data %>%
    mutate(spatial_barcode_ori = if ("spatial_bc" %in% colnames(.)) .data$spatial_bc else rownames(filtered@meta.data)) %>%
    select(.data$spatial_barcode_ori, .data$GC_part, .data$GC_Maturity_Label)

  p_spatial <- plot_spatial_gc_parts(filtered, group_col = "GC_part", label = TRUE)
  save_panel(p_spatial, paste0("figure4A_", cfg$output_prefix, "_gc_spatial.pdf"), width = 7, height = 6)

  p_dot <- plot_gc_maturity_dotplot(filtered, group_col = "GC_part")
  save_panel(p_dot, paste0("figure4D_", cfg$output_prefix, "_maturity_dotplot.pdf"), width = 4.2, height = 3.5)

  p_cgene <- plot_cgene_bar_by_gc(cfg$igh_alignment_csv, metadata_df)
  save_panel(p_cgene, paste0("figureS8A_", cfg$output_prefix, "_cgene_by_gc.pdf"), width = 4.8, height = 3)

  class_obj <- classify_igh_clones(
    filtered,
    cfg$igh_assay,
    gc_col = "GC_part",
    umi_quantile = clone_expanded_umi_quantile
  )
  p_class <- plot_clone_class_bar_by_gc(class_obj$long, class_obj$summary)
  save_panel(p_class, paste0("figureS8B_", cfg$output_prefix, "_class_by_gc.pdf"), width = 4.8, height = 3)

  save_table(scored$summary, paste0("figure4_", cfg$output_prefix, "_gc_maturity_summary.csv"))
  save_table(scored$parameters, paste0("figure4_", cfg$output_prefix, "_gc_maturity_parameters.csv"))
  save_table(class_obj$parameters, paste0("figure4_", cfg$output_prefix, "_clone_class_parameters.csv"))
  invisible(list(object = filtered, metadata = metadata_df, clone_class = class_obj))
}
