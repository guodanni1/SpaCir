# Figure 2H: IGH somatic hypermutation (SHM) rate by tissue region.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

find_shm_column <- function(df) {
  candidates <- c("SHM_rate", "VDJ_mutation_rate", "VDJC_mutation_rate", "allVAlignments_mutation_rate")
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) {
    stop("No SHM/mutation-rate column found. Expected one of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit
}

plot_shm_by_region <- function(sample_id, regions = NULL) {
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  vdj <- read_csv_if_exists(file.path(repo_root, cfg$igh_vdj_csv), stringsAsFactors = FALSE)
  shm_col <- find_shm_column(vdj)

  region_df <- meta %>%
    mutate(Spot = rownames(meta), Region = .data[[cfg$region_column]]) %>%
    select(Spot, Region)

  dat <- vdj %>%
    transmute(Spot = .data[[cfg$spot_column]], SHM = suppressWarnings(as.numeric(.data[[shm_col]]))) %>%
    inner_join(region_df, by = "Spot") %>%
    filter(!is.na(SHM), !is.na(Region))
  if (!is.null(regions)) {
    dat <- dat %>% filter(Region %in% regions)
  }

  p <- ggplot(dat, aes(x = Region, y = SHM, fill = Region)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
    scale_fill_manual(values = region_colors, na.value = "#BDBDBD") +
    labs(title = "SHM", x = NULL, y = "SHM rate") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  save_panel(p, paste0("figure2H_", display_sample(cfg$sample_id), "_IGH_SHM_region.pdf"), 4, 3)
  invisible(dat)
}

# Example:
# plot_shm_by_region("Tonsil2", c("Follicular LZ", "Inter LZ-DZ", "Follicular DZ"))
