# Figure 2I/K: Spearman correlations between cell-type abundance and VDJ metrics.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

cor_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

compute_abundance_metric_correlation <- function(sample_id, receptor = c("IGH", "TRB")) {
  receptor <- match.arg(receptor)
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  vdj_path <- if (receptor == "IGH") cfg$igh_vdj_csv else cfg$trb_vdj_csv
  vdj <- read_csv_if_exists(file.path(repo_root, vdj_path), stringsAsFactors = FALSE)

  abundance <- cell_abundance_long(meta) %>% select(Spot, CellType, Abundance)
  metrics <- vdj_spot_metrics(vdj, spot_col = cfg$spot_column)

  abundance %>%
    inner_join(metrics, by = "Spot") %>%
    pivot_longer(c(Shannon, MeanLevenshtein, UniqueClonotypes), names_to = "Metric", values_to = "MetricValue") %>%
    group_by(CellType, Metric) %>%
    summarise(
      n = sum(complete.cases(Abundance, MetricValue)),
      rho = ifelse(n >= 3, cor(Abundance, MetricValue, method = "spearman", use = "complete.obs"), NA_real_),
      p = ifelse(n >= 3, suppressWarnings(cor.test(Abundance, MetricValue, method = "spearman")$p.value), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p, method = "BH"), stars = cor_stars(p_adj), Receptor = receptor)
}

plot_correlation_heatmap <- function(sample_id, receptor = c("IGH", "TRB"), metric = "Shannon",
                                     output_tag = "figure2I_K") {
  receptor <- match.arg(receptor)
  cor_df <- compute_abundance_metric_correlation(sample_id, receptor) %>%
    filter(Metric == metric)

  cor_df$CellType <- factor(cor_df$CellType, levels = rev(intersect(names(cell_type_colors), unique(cor_df$CellType))))
  p <- ggplot(cor_df, aes(x = Receptor, y = CellType, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = stars), size = 2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Spearman r") +
    labs(x = NULL, y = NULL, title = paste(display_sample(sample_id), receptor, metric)) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank())

  save_panel(p, paste0(output_tag, "_", display_sample(sample_id), "_", receptor, "_", metric, "_heatmap.pdf"), 2.6, 6)
  write.csv(cor_df, file.path(output_dir, paste0(output_tag, "_", display_sample(sample_id), "_", receptor, "_", metric, "_correlations.csv")), row.names = FALSE)
  invisible(cor_df)
}

# Examples:
# plot_correlation_heatmap("Tonsil2", "IGH", "Shannon")
# plot_correlation_heatmap("Tonsil2", "TRB", "Shannon")
