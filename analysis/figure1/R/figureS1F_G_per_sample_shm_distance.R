# Figure S1F/G: per-sample SHM and nearest-neighbor distance panels.

source("analysis/figure1/R/figure1H_shm_rate.R")
source("analysis/figure1/R/figure1J_K_M_N_spatial_expansion.R")

# Per-sample SHM violin plots.
for (s in unique(data$Sample)) {
  p <- plot_violin(data |> filter(Sample == s), "VDJ_mutation_rate", "SHM rate")
  save_panel(p, paste0("figureS1F_", s, "_IGH_SHM_rate.pdf"), width = 4, height = 4)
}

# Per-sample nearest-neighbor distance plots for IGH and TRB.
for (chain in c("IGH", "TRB")) {
  for (s in unique(all_dist$Sample)) {
    df <- all_dist |> filter(Chain == chain, Sample == s, Classification != "Not Expanded")
    p <- ggplot(df, aes(x = Classification, y = MeanNearestNeighborDistance, fill = Classification)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21) +
      geom_signif(comparisons = list(c("Small", "Medium"), c("Medium", "Hyperexpanded")),
                  map_signif_level = TRUE) +
      scale_fill_manual(values = c("Small" = "orange", "Medium" = "blue", "Hyperexpanded" = "red")) +
      labs(title = paste(s, chain), x = "", y = "Mean Nearest Neighbor Distance") +
      theme_classic() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    save_panel(p, paste0("figureS1G_", s, "_", chain, "_nearest_neighbor_distance.pdf"), width = 4, height = 4)
  }
}
