# Figure 2M: pathway-specific CellChat network and ligand-receptor contribution.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(CellChat)
})

plot_cellchat_pathway <- function(cellchat_rds, sample_id, pathway = "ICOS", output_tag = "figure2M") {
  cellchat <- readRDS(cellchat_rds)
  label <- display_sample(sample_id)

  # The manuscript Figure 2M highlights ICOS; the supplementary panels use the
  # same logic for BAFF, MHC-II, CXCL and other pathways.
  pdf(file.path(output_dir, paste0(output_tag, "_", label, "_", pathway, "_aggregate_network.pdf")), width = 5, height = 4)
  netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  dev.off()

  contrib <- netAnalysis_contribution(cellchat, signaling = pathway)
  ggsave(
    filename = file.path(output_dir, paste0(output_tag, "_", label, "_", pathway, "_lr_contribution.pdf")),
    plot = contrib,
    width = 3,
    height = 2.2,
    limitsize = FALSE
  )

  invisible(contrib)
}

# Example:
# plot_cellchat_pathway("analysis/figure2/outputs/figure2L_LN2_cellchat.rds", "Tonsil2", "ICOS")
