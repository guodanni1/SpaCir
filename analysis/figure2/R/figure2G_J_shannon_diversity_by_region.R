# Figure 2G/J: clonotype Shannon diversity by tissue region.
# The same UMI-weighted Shannon calculation is used for IGH and TRB.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

plot_shannon_by_region <- function(sample_id, receptor = c("IGH", "TRB"), regions = NULL, output_tag = "figure2G") {
  receptor <- match.arg(receptor)
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  vdj_path <- if (receptor == "IGH") cfg$igh_vdj_csv else cfg$trb_vdj_csv
  vdj <- read_csv_if_exists(file.path(repo_root, vdj_path), stringsAsFactors = FALSE)

  metrics <- vdj_spot_metrics(vdj, spot_col = cfg$spot_column)
  region_df <- meta %>%
    mutate(Spot = rownames(meta), Region = .data[[cfg$region_column]]) %>%
    select(Spot, Region)

  dat <- metrics %>%
    inner_join(region_df, by = "Spot") %>%
    filter(!is.na(Region))
  if (!is.null(regions)) {
    dat <- dat %>% filter(Region %in% regions)
  }
  dat <- dat %>% mutate(Region = factor(Region, levels = intersect(names(region_colors), unique(Region))))

  p <- ggplot(dat, aes(x = Region, y = Shannon, fill = Region)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
    scale_fill_manual(values = region_colors, na.value = "#BDBDBD") +
    labs(title = paste(display_sample(cfg$sample_id), receptor, "Shannon Diversity"), x = NULL, y = "Shannon Diversity") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  save_panel(p, paste0(output_tag, "_", display_sample(cfg$sample_id), "_", receptor, "_shannon.pdf"), 5, 3)
  write.csv(dat, file.path(output_dir, paste0(output_tag, "_", display_sample(cfg$sample_id), "_", receptor, "_spot_metrics.csv")), row.names = FALSE)
  invisible(dat)
}

# Examples:
# plot_shannon_by_region("Tonsil2", "IGH")
# plot_shannon_by_region("Tonsil2", "TRB")
