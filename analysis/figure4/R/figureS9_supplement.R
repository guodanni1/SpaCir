# Figure S9 entry points.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "figure4A_D_gc_maturity_annotation.R"))
source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "figure4B_C_G_H_I_J_K_L_clonotype_maturity.R"))
source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "figure4M_N_O_celltype_trends.R"))

run_figureS9_panels <- function(sample_id = "LN2") {
  annotation <- run_figure4_gc_annotation_panels(sample_id)
  clonotype <- run_figure4_clonotype_panels(sample_id)
  celltype <- run_figure4_celltype_trend_panels(sample_id)
  invisible(list(annotation = annotation, clonotype = clonotype, celltype = celltype))
}
