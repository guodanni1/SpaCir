# Figure S8 entry points.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "figure4A_D_gc_maturity_annotation.R"))
source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure4", "R", "figure4E_F_S8D_E_wgcna_go.R"))

run_figureS8_panels <- function(sample_id = "Tonsil2") {
  annotation <- run_figure4_gc_annotation_panels(sample_id)
  wgcna <- run_figure4_wgcna_go_panels(sample_id)
  invisible(list(annotation = annotation, wgcna = wgcna))
}
