# Figure 2D-F: compare cell2location cell-type abundance across region groups.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

pairwise_wilcox_bh <- function(df, group_col, value_col) {
  groups <- sort(unique(df[[group_col]]))
  pairs <- combn(groups, 2, simplify = FALSE)
  stats <- lapply(pairs, function(pair) {
    x <- df[df[[group_col]] == pair[1], value_col, drop = TRUE]
    y <- df[df[[group_col]] == pair[2], value_col, drop = TRUE]
    p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    data.frame(group1 = pair[1], group2 = pair[2], p = p)
  })
  bind_rows(stats) %>% mutate(p_adj = p.adjust(p, method = "BH"))
}

plot_abundance_violin <- function(sample_id, regions, cell_types, output_tag) {
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  region_col <- cfg$region_column

  dat <- cell_abundance_long(meta) %>%
    mutate(Region = .data[[region_col]]) %>%
    filter(Region %in% regions, CellType %in% cell_types) %>%
    mutate(
      Region = factor(Region, levels = regions),
      CellType = factor(CellType, levels = cell_types)
    )

  # Wilcoxon + BH follows the original plotting logic for region comparisons.
  stats <- dat %>%
    group_by(CellType) %>%
    group_modify(~ pairwise_wilcox_bh(.x, "Region", "Abundance")) %>%
    ungroup()

  p <- ggplot(dat, aes(x = Region, y = Abundance, fill = Region)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
    facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = region_colors, na.value = "#BDBDBD") +
    labs(x = NULL, y = "Cell abundance") +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

  save_panel(p, paste0(output_tag, "_", display_sample(cfg$sample_id), "_abundance_violin.pdf"), 10, 3.2)
  write.csv(stats, file.path(output_dir, paste0(output_tag, "_", display_sample(cfg$sample_id), "_wilcox_BH.csv")), row.names = FALSE)
  invisible(list(plot = p, stats = stats))
}

# Figure 2D follicular-region examples.
# plot_abundance_violin("Tonsil2", c("Follicular LZ", "Inter LZ-DZ", "Follicular DZ", "Mantle zone"),
#                       c("B Cycling", "B GC DZ", "B GC LZ", "B GC prePB", "B plasma", "T CD4+ TfH GC"),
#                       "figure2D")
#
# Figure 2E inter-follicular examples.
# plot_abundance_violin("Tonsil2", c("Inter-Follicular CD5+", "Inter-Follicular CD5-"),
#                       c("T CD8+ naive", "Macrophages M1", "B plasma", "B Cycling"),
#                       "figure2E")
#
# Figure 2F stromal/parenchymal examples.
# plot_abundance_violin("Tonsil2", c("Muscle/Mesenchymal", "Vascular/Mesenchymal", "Vascular"),
#                       c("VSMC", "B mem", "B plasma", "Endo", "DC cDC2", "Monocytes", "Mast"),
#                       "figure2F")
