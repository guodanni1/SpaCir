# Figure 2J/K: random cell-type assignment and B/T diversity correlation.
#
# This script wraps the shared helper functions in
# random_celltype_diversity_helpers.R into a panel-level workflow:
#   1. Use cell2location abundance as spot-level assignment probabilities.
#   2. Randomly assign each VDJ clonotype record to a B or T cell subtype.
#   3. Compute clonotype Shannon diversity for each Spot x CellType.
#   4. Plot per-cell-type Shannon violin plots and B-by-T Spearman heatmaps.
#
# The receptor-to-cell-family mapping and statistical settings are recorded in:
#   config/random_celltype_parameters.csv

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))
source(file.path(figure2_dir, "R", "random_celltype_diversity_helpers.R"))

bt_correlation_heatmap <- function(b_shannon, t_shannon, output_tag) {
  b_wide <- b_shannon %>%
    rename(BCellType = CellType, BShannon = Shannon) %>%
    select(Spot, BCellType, BShannon)
  t_wide <- t_shannon %>%
    rename(TCellType = CellType, TShannon = Shannon) %>%
    select(Spot, TCellType, TShannon)

  cor_tbl <- inner_join(b_wide, t_wide, by = "Spot") %>%
    group_by(BCellType, TCellType) %>%
    summarise(
      n = sum(complete.cases(BShannon, TShannon)),
      rho = ifelse(n >= 3, cor(BShannon, TShannon, method = "spearman", use = "complete.obs"), NA_real_),
      p = ifelse(n >= 3, suppressWarnings(cor.test(BShannon, TShannon, method = "spearman")$p.value), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      stars = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  p <- ggplot(cor_tbl, aes(x = TCellType, y = BCellType, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = stars), size = 2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Spearman r") +
    labs(x = "T cell type", y = "B cell type") +
    theme_minimal(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  save_panel(p, paste0(output_tag, "_BT_random_celltype_shannon_correlation.pdf"), 5, 4)
  write.csv(cor_tbl, file.path(output_dir, paste0(output_tag, "_BT_random_celltype_shannon_correlation.csv")), row.names = FALSE)
  invisible(cor_tbl)
}

run_random_celltype_assignment <- function(sample_id,
                                           b_receptor = "IGH",
                                           t_receptor = "TRB",
                                           b_pattern = "^B ",
                                           t_pattern = "^T ",
                                           seed = 1,
                                           output_tag = "figure2J_K") {
  cfg <- sample_row(sample_id)
  meta <- read_csv_if_exists(file.path(repo_root, cfg$cell2location_metadata), stringsAsFactors = FALSE)
  b_vdj_path <- if (b_receptor == "IGH") cfg$igh_vdj_csv else cfg$trb_vdj_csv
  t_vdj_path <- if (t_receptor == "IGH") cfg$igh_vdj_csv else cfg$trb_vdj_csv
  b_vdj <- read_csv_if_exists(file.path(repo_root, b_vdj_path), stringsAsFactors = FALSE)
  t_vdj <- read_csv_if_exists(file.path(repo_root, t_vdj_path), stringsAsFactors = FALSE)

  # B and T assignments use separate abundance patterns so the random draw only
  # chooses among the biologically relevant cell-type family.
  b_assigned <- assign_clonotypes_to_celltypes(
    b_vdj, meta, celltype_pattern = b_pattern,
    spot_col = cfg$spot_column, seed = seed
  )
  t_assigned <- assign_clonotypes_to_celltypes(
    t_vdj, meta, celltype_pattern = t_pattern,
    spot_col = cfg$spot_column, seed = seed
  )

  b_shannon <- celltype_shannon(b_assigned)
  t_shannon <- celltype_shannon(t_assigned)

  prefix <- paste(output_tag, display_sample(sample_id), sep = "_")
  write.csv(b_assigned, file.path(output_dir, paste0(prefix, "_B_assigned_clonotypes.csv")), row.names = FALSE)
  write.csv(t_assigned, file.path(output_dir, paste0(prefix, "_T_assigned_clonotypes.csv")), row.names = FALSE)
  write.csv(b_shannon, file.path(output_dir, paste0(prefix, "_B_celltype_shannon.csv")), row.names = FALSE)
  write.csv(t_shannon, file.path(output_dir, paste0(prefix, "_T_celltype_shannon.csv")), row.names = FALSE)

  plot_celltype_shannon_violin(b_shannon, paste0(prefix, "_B"))
  plot_celltype_shannon_violin(t_shannon, paste0(prefix, "_T"))
  bt_correlation_heatmap(b_shannon, t_shannon, prefix)

  invisible(list(
    b_assigned = b_assigned,
    t_assigned = t_assigned,
    b_shannon = b_shannon,
    t_shannon = t_shannon
  ))
}

# Examples:
# run_random_celltype_assignment("LN2", b_receptor = "IGH", t_receptor = "TRB", output_tag = "figure2J_K")
# run_random_celltype_assignment("LN2", b_receptor = "IGH", t_receptor = "TRB", output_tag = "figureS3J_K")
