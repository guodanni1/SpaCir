# Figure 5O-R and Figure S12J-L: clonotype diversity across SHM groups.
#
# Only normalized Shannon diversity is exported for Shannon panels. Mean
# Levenshtein distance is not abundance-normalized because it is a sequence
# distance metric rather than a count or diversity count.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

run_clonotype_diversity_panels <- function(sample_id = "Tonsil2",
                                           output_tag = "figure5",
                                           group_levels = shm_levels) {
  obj_info <- load_sample_object(sample_id, filtered = TRUE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  shm_meta <- prepare_spot_mutation_metadata(seurat_obj, cfg)
  abundance_meta <- calc_abundance_totals(seurat_obj@meta.data)

  igh_metrics <- calc_diversity_norm(seurat_obj, cfg$igh_assay, abundance_meta, "B_total") %>%
    left_join(shm_meta %>% select(.data$spatial_barcode_ori, .data$SHM_VDJ),
              by = c("Spot" = "spatial_barcode_ori")) %>%
    filter(.data$SHM_VDJ %in% group_levels) %>%
    mutate(SHM_VDJ = factor(.data$SHM_VDJ, levels = group_levels))

  trb_metrics <- calc_diversity_norm(seurat_obj, cfg$trb_assay, abundance_meta, "T_total") %>%
    left_join(shm_meta %>% select(.data$spatial_barcode_ori, .data$SHM_VDJ),
              by = c("Spot" = "spatial_barcode_ori")) %>%
    filter(.data$SHM_VDJ %in% group_levels) %>%
    mutate(SHM_VDJ = factor(.data$SHM_VDJ, levels = group_levels))

  p_igh_shannon <- plot_violin_metric(igh_metrics, "SHM_VDJ", "Shannon_norm", "IGH Shannon Diversity Index")
  save_panel(p_igh_shannon, paste0(output_tag, "O_", cfg$output_prefix, "_IGH_shannon_norm.pdf"), width = 3.1, height = 3)

  p_igh_lv <- plot_violin_metric(igh_metrics, "SHM_VDJ", "MeanLevenshtein", "IGH Mean Levenshtein Distance")
  save_panel(p_igh_lv, paste0(output_tag, "P_", cfg$output_prefix, "_IGH_levenshtein.pdf"), width = 3.1, height = 3)

  p_trb_shannon <- plot_violin_metric(trb_metrics, "SHM_VDJ", "Shannon_norm", "TRB Shannon Diversity Index")
  save_panel(p_trb_shannon, paste0(output_tag, "Q_", cfg$output_prefix, "_TRB_shannon_norm.pdf"), width = 3.1, height = 3)

  cgene_data <- prepare_cgene_trend(cfg$igh_alignment_csv, shm_meta, group_levels = group_levels)
  p_cgene <- plot_cgene_trend(cgene_data)
  save_panel(p_cgene, paste0(output_tag, "R_", cfg$output_prefix, "_cgene_trend.pdf"), width = 3.4, height = 3)

  save_table(igh_metrics, paste0(output_tag, "_", cfg$output_prefix, "_IGH_diversity_norm_metrics.csv"))
  save_table(trb_metrics, paste0(output_tag, "_", cfg$output_prefix, "_TRB_diversity_norm_metrics.csv"))
  save_table(cgene_data$trend, paste0(output_tag, "_", cfg$output_prefix, "_cgene_trend.csv"))

  invisible(list(igh = igh_metrics, trb = trb_metrics, cgene = cgene_data))
}

# Examples:
# run_clonotype_diversity_panels("Tonsil2", output_tag = "figure5", group_levels = shm_levels)
# run_clonotype_diversity_panels("LN2", output_tag = "figureS12", group_levels = c("Low_Mutation", "High_Mutation"))
