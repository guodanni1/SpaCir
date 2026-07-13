# Figure S1C: IGH and TRB sample coverage curves.

source("analysis/figure1/R/00_setup.R")
suppressPackageStartupMessages({
  library(iNEXT)
})

get_umi_abundance <- function(obj, assay_name) {
  mat <- as.matrix(obj@assays[[assay_name]]@data)
  umi <- rowSums(mat)
  umi[umi > 0]
}

plot_coverage <- function(chain) {
  assay_name <- paste0("aaSeqCDR3_Lv1_UMI3_", chain, "_20G")
  abundance_list <- lapply(seq_len(nrow(sample_config)), function(i) {
    object_name <- sample_config$object_name[i]
    obj <- get(object_name)
    get_umi_abundance(obj, assay_name)
  })
  names(abundance_list) <- paste0(sample_config$display_sample, "_", chain)

  out <- iNEXT(abundance_list, q = 0, datatype = "abundance")
  p <- ggiNEXT(out, type = 2) +
    labs(title = paste(chain, "Sample Coverage Curve"),
         x = "Cumulative UMI Count",
         y = "Sample Coverage") +
    theme_minimal()
  save_panel(p, paste0("figureS1C_", chain, "_sample_coverage_curve.pdf"), width = 5, height = 4)
}

plot_coverage("IGH")
plot_coverage("TRB")
