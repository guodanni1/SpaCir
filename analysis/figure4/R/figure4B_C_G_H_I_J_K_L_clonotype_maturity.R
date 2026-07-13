# Figure 4 clonotype, SHM, normalized diversity, and trend panels.
#
# Figure 4B clone-class definition:
#   A (Private Expanded): IGH clonotype detected in exactly one GC site and
#                         total UMI above the sample-specific 95th percentile.
#   B (Shared Expanded): IGH clonotype detected in more than one GC site and
#                        total UMI above the same 95th percentile.
#   C (Private Unexpanded): IGH clonotype detected in exactly one GC site and
#                           total UMI at or below the 95th percentile.
#   D (Shared Unexpanded): IGH clonotype detected in more than one GC site and
#                          total UMI at or below the 95th percentile.
# Diversity panels use abundance-normalized values only.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "00_setup.R"))

run_figure4_clonotype_panels <- function(sample_id = "Tonsil2") {
  obj_info <- load_sample_objects(sample_id)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$filtered_object %||% obj_info$object
  if (is.null(seurat_obj)) {
    stop("A Seurat object is required for clonotype maturity panels.", call. = FALSE)
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

  metadata_df <- seurat_obj@meta.data %>%
    mutate(spatial_barcode_ori = if ("spatial_bc" %in% colnames(.)) .data$spatial_bc else rownames(seurat_obj@meta.data)) %>%
    select(.data$spatial_barcode_ori, .data$GC_part, .data$GC_Maturity_Label)

  class_obj <- classify_igh_clones(
    seurat_obj,
    cfg$igh_assay,
    gc_col = "GC_part",
    umi_quantile = clone_expanded_umi_quantile
  )
  save_table(class_obj$summary, paste0("figure4_", cfg$output_prefix, "_igh_clone_class_summary.csv"))
  save_table(class_obj$parameters, paste0("figure4_", cfg$output_prefix, "_igh_clone_class_parameters.csv"))

  p_detection <- plot_clone_detection(class_obj$summary, title = paste(cfg$display_sample, "IGH Clone Detection in GC Sites"))
  save_panel(p_detection, paste0("figure4B_", cfg$output_prefix, "_clone_detection.pdf"), width = 5.3, height = 4)

  shm_by_class <- read_csv_if_exists(cfg$igh_alignment_csv, stringsAsFactors = FALSE) %>%
    distinct() %>%
    left_join(
      class_obj$summary %>% transmute(aaSeqCDR3_most = .data$IGH_clonotype, class = .data$class),
      by = "aaSeqCDR3_most"
    ) %>%
    filter(!is.na(.data$class)) %>%
    mutate(SHM = if ("VDJ_mutation_rate" %in% colnames(.)) as.numeric(.data$VDJ_mutation_rate) else NA_real_)
  p_shm_class <- plot_violin_with_p(shm_by_class, "class", "SHM", clone_class_colors, y_label = "SHM Frequency")
  save_panel(p_shm_class, paste0("figure4C_", cfg$output_prefix, "_shm_by_class.pdf"), width = 4.6, height = 3.2)

  meta_norm <- calc_abundance_totals(seurat_obj@meta.data)
  igh_metrics <- calc_diversity_norm(seurat_obj, cfg$igh_assay, meta_norm, "B_total")
  trb_metrics <- calc_diversity_norm(seurat_obj, cfg$trb_assay, meta_norm, "T_total")
  p_div <- plot_diversity_norm_pair(igh_metrics, trb_metrics)
  save_panel(p_div, paste0("figure4G_", cfg$output_prefix, "_diversity_norm.pdf"), width = 4.8, height = 3.4)

  p_cdf <- plot_cumulative_clone_frequency(seurat_obj, cfg$igh_assay)
  save_panel(p_cdf, paste0("figure4H_", cfg$output_prefix, "_igh_cumulative_frequency.pdf"), width = 3.8, height = 3)

  p_clone_comp <- plot_late_clone_composition(seurat_obj, cfg$igh_assay, top_n = 100)
  save_panel(p_clone_comp, paste0("figure4I_", cfg$output_prefix, "_igh_clone_composition.pdf"), width = 2.8, height = 3)

  shm_by_maturity <- prepare_shm_table(cfg$igh_alignment_csv, metadata_df, group_col = "GC_Maturity_Label")
  p_shm_maturity <- plot_violin_with_p(shm_by_maturity, "GC_Maturity_Label", "SHM", maturity_colors, y_label = "SHM")
  save_panel(p_shm_maturity, paste0("figure4J_", cfg$output_prefix, "_shm_by_maturity.pdf"), width = 3.4, height = 3)

  class_long <- matrix_to_long(get_assay_matrix(seurat_obj, cfg$igh_assay), feature_col = "IGH_clonotype", value_col = "UMI") %>%
    left_join(seurat_obj@meta.data %>% mutate(spot = rownames(seurat_obj@meta.data)) %>% select(.data$spot, .data$GC_Maturity_Label), by = "spot") %>%
    left_join(class_obj$summary %>% select(.data$IGH_clonotype, .data$class), by = "IGH_clonotype") %>%
    filter(!is.na(.data$class), .data$GC_Maturity_Label %in% maturity_levels)
  class_trend <- class_long %>%
    group_by(.data$GC_Maturity_Label, .data$class) %>%
    summarise(class_UMI = sum(.data$UMI), .groups = "drop") %>%
    group_by(.data$GC_Maturity_Label) %>%
    mutate(Percentage = 100 * .data$class_UMI / sum(.data$class_UMI)) %>%
    ungroup()
  class_raw <- class_long %>%
    group_by(.data$spot, .data$GC_Maturity_Label, .data$class) %>%
    summarise(Percentage = sum(.data$UMI), .groups = "drop")
  p_class_trend <- plot_two_state_line_with_stars(
    point_df = class_trend,
    raw_df = class_raw,
    group_col = GC_Maturity_Label,
    value_col = Percentage,
    line_group_col = class,
    y_label = "Class Percentage (%)",
    color_values = clone_class_colors,
    title = "IGH Class"
  )
  save_panel(p_class_trend, paste0("figure4K_", cfg$output_prefix, "_igh_class_trend.pdf"), width = 4.4, height = 3)

  cgene_data <- prepare_cgene_trend(cfg$igh_alignment_csv, metadata_df)
  cgene_colors <- setNames(rep(cell_type_palette, length.out = length(unique(cgene_data$trend$bestCGene))),
                           sort(unique(cgene_data$trend$bestCGene)))
  p_cgene <- plot_two_state_line_with_stars(
    point_df = cgene_data$trend,
    raw_df = cgene_data$raw,
    group_col = GC_Maturity_Label,
    value_col = UMI_percent,
    line_group_col = bestCGene,
    y_label = "Relative Abundance (%)",
    color_values = cgene_colors,
    title = NULL
  )
  save_panel(p_cgene, paste0("figure4L_", cfg$output_prefix, "_cgene_trend.pdf"), width = 4.3, height = 3)

  invisible(list(clone_class = class_obj, igh_metrics = igh_metrics, trb_metrics = trb_metrics))
}
