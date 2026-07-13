# Figure 5E/F and Figure S12E/F: module scores versus SHM rate.
#
# High and Low SHM module genes are declared in config/shm_gene_sets.csv.
# Correlations use Spearman tests and the plotted trend line is linear.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

run_module_correlation_panels <- function(sample_id = "Tonsil2", output_tag = "figure5") {
  obj_info <- load_sample_object(sample_id, filtered = TRUE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  shm_meta <- prepare_spot_mutation_metadata(seurat_obj, cfg)
  seurat_obj <- add_shm_module_scores(seurat_obj, assay = "SCT", nbin = 2)

  p_high <- plot_module_correlation(
    seurat_obj,
    shm_meta,
    score_col = "Score_High1",
    color_value = "#E3170D",
    x_label = "High SHM Gene Expression Score"
  )
  save_panel(p_high, paste0(output_tag, "E_", cfg$output_prefix, "_high_module_correlation.pdf"), width = 3.2, height = 3)

  p_low <- plot_module_correlation(
    seurat_obj,
    shm_meta,
    score_col = "Score_Low1",
    color_value = "#8B58A4",
    x_label = "Low SHM Gene Expression Score"
  )
  save_panel(p_low, paste0(output_tag, "F_", cfg$output_prefix, "_low_module_correlation.pdf"), width = 3.2, height = 3)

  invisible(list(object = seurat_obj, shm_meta = shm_meta))
}

# Examples:
# run_module_correlation_panels("Tonsil2", output_tag = "figure5")
# run_module_correlation_panels("LN2", output_tag = "figureS12")
