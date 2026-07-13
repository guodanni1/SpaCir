# Figure 5A/B and Figure S12A/B: SHM group definition and spatial maps.
#
# SHM groups are defined from avg_15G_VDJ_mutation_rate:
#   No_Mutation: <= 0 or missing
#   Low_Mutation: > 0 and <= 0.05
#   High_Mutation: > 0.05

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

run_shm_group_panels <- function(sample_id = "Tonsil2", output_tag = "figure5") {
  obj_info <- load_sample_object(sample_id, filtered = TRUE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  shm_meta <- prepare_spot_mutation_metadata(seurat_obj, cfg)
  seurat_obj <- attach_shm_group(seurat_obj, shm_meta)

  p_rank <- plot_ranked_shm(shm_meta, cfg$display_sample)
  save_panel(p_rank, paste0(output_tag, "A_", cfg$output_prefix, "_shm_ranked.pdf"), width = 3.4, height = 3)

  p_spatial <- SpatialDimPlot(
    seurat_obj,
    group.by = "SHM_VDJ",
    cols = c(shm_colors, "Other" = "#D9D9D9"),
    pt.size.factor = 1.6
  ) +
    ggtitle(paste(cfg$display_sample, "SHM Group")) +
    theme(legend.position = "right")
  save_panel(p_spatial, paste0(output_tag, "B_", cfg$output_prefix, "_shm_spatial.pdf"), width = 4.2, height = 3.8)

  save_table(shm_meta, paste0(output_tag, "_", cfg$output_prefix, "_shm_group_metadata.csv"))
  invisible(list(object = seurat_obj, shm_meta = shm_meta))
}

# Examples:
# run_shm_group_panels("Tonsil2", output_tag = "figure5")
# run_shm_group_panels("LN2", output_tag = "figureS12")
