# Figure 2B: stacked cell-type composition by annotated tissue region.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

plot_celltype_composition <- function(sample_id) {
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  region_col <- cfg$region_column
  if (!region_col %in% names(meta)) {
    stop("Region column not found in cell2location metadata: ", region_col, call. = FALSE)
  }

  comp <- cell_abundance_long(meta) %>%
    mutate(Region = .data[[region_col]]) %>%
    group_by(Region, CellType) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Region) %>%
    mutate(Percentage = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()

  p <- ggplot(comp, aes(x = Region, y = Percentage, fill = CellType)) +
    geom_col(width = 0.9, color = NA) +
    scale_fill_manual(values = cell_type_colors, na.value = "#BDBDBD") +
    labs(x = NULL, y = "Percentage (%)", fill = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

  save_panel(p, paste0("figure2B_", display_sample(cfg$sample_id), "_celltype_composition.pdf"), 8, 4)
  invisible(comp)
}

# Example:
# plot_celltype_composition("Tonsil2")
