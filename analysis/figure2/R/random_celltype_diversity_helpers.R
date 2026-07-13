# Helper functions for the random cell-type assignment analyses used in
# Figure 2J/K and Figure S2-S4.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

assign_clonotypes_to_celltypes <- function(vdj_df, meta_df, celltype_pattern,
                                           spot_col = "spatial_barcode_ori",
                                           clonotype_col = "aaSeqCDR3_most",
                                           umi_col = "UMI_count",
                                           seed = 1) {
  set.seed(seed)
  abundance <- cell_abundance_long(meta_df) %>%
    filter(grepl(celltype_pattern, CellType)) %>%
    group_by(Spot) %>%
    mutate(Probability = Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()

  vdj_df %>%
    transmute(
      Spot = .data[[spot_col]],
      Clonotype = .data[[clonotype_col]],
      UMI = suppressWarnings(as.numeric(.data[[umi_col]]))
    ) %>%
    filter(Spot %in% abundance$Spot) %>%
    group_by(Spot) %>%
    group_modify(function(.x, .y) {
      probs <- abundance %>% filter(Spot == .y$Spot[[1]], Probability > 0)
      if (nrow(probs) == 0) {
        return(tibble())
      }
      .x$CellType <- sample(probs$CellType, nrow(.x), replace = TRUE, prob = probs$Probability)
      .x
    }) %>%
    ungroup()
}

celltype_shannon <- function(assigned_df) {
  assigned_df %>%
    group_by(Spot, CellType, Clonotype) %>%
    summarise(UMI = sum(replace_na(UMI, 1)), .groups = "drop") %>%
    group_by(Spot, CellType) %>%
    summarise(Shannon = weighted_shannon(UMI), .groups = "drop")
}

plot_celltype_shannon_violin <- function(shannon_df, output_tag) {
  shannon_df$CellType <- factor(shannon_df$CellType, levels = intersect(names(cell_type_colors), unique(shannon_df$CellType)))
  p <- ggplot(shannon_df, aes(x = CellType, y = Shannon, fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
    scale_fill_manual(values = cell_type_colors, na.value = "#BDBDBD") +
    labs(x = NULL, y = "Shannon diversity") +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_panel(p, paste0(output_tag, "_random_celltype_shannon.pdf"), 7, 4)
  invisible(p)
}
