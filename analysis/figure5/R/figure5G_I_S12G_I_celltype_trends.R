# Figure 5G-I and Figure S12G-I: cell-type abundance trends across SHM groups.
#
# Curves show mean per-spot cell-type percentages. Significance is tested on
# the original cell2location abundance values.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

run_celltype_trend_panels <- function(sample_id = "Tonsil2",
                                      output_tag = "figure5",
                                      group_levels = shm_levels) {
  obj_info <- load_sample_object(sample_id, filtered = TRUE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  shm_meta <- prepare_spot_mutation_metadata(seurat_obj, cfg)
  cell_frac <- prepare_celltype_fraction(seurat_obj@meta.data, shm_meta)

  p_b <- plot_celltype_trend(cell_frac, "B cell", group_levels = group_levels)
  save_panel(p_b, paste0(output_tag, "G_", cfg$output_prefix, "_b_cell_trend.pdf"), width = 4.8, height = 3.2)

  p_t <- plot_celltype_trend(cell_frac, "T cell", group_levels = group_levels)
  save_panel(p_t, paste0(output_tag, "H_", cfg$output_prefix, "_t_cell_trend.pdf"), width = 4.8, height = 3.2)

  p_o <- plot_celltype_trend(cell_frac, "Other cell", group_levels = group_levels)
  save_panel(p_o, paste0(output_tag, "I_", cfg$output_prefix, "_other_cell_trend.pdf"), width = 5.5, height = 3.2)

  save_table(cell_frac, paste0(output_tag, "_", cfg$output_prefix, "_celltype_fraction_by_shm.csv"))
  invisible(cell_frac)
}

# Examples:
# run_celltype_trend_panels("Tonsil2", output_tag = "figure5", group_levels = shm_levels)
# run_celltype_trend_panels("LN2", output_tag = "figureS12", group_levels = c("Low_Mutation", "High_Mutation"))
