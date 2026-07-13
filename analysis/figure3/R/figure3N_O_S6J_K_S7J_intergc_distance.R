# Figure 3N/O and supplementary InterGC count/distance panels.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))
source(file.path(figure3_dir, "R", "figure3H_I_spatial_clonefamily_branches.R"))

plot_intergc_origin_counts <- function(sample_id, edge_df = NULL, output_tag = "figure3N") {
  if (is.null(edge_df)) {
    edge_df <- read_tree_edges(sample_id)
  }
  dat <- prepare_branch_edges(edge_df) %>%
    filter(GC_transition_type == "InterGC") %>%
    mutate(OriginZone = collapse_zone(parent_GC_NewClusters)) %>%
    filter(OriginZone != "Other") %>%
    count(OriginZone, name = "Count")

  p <- ggplot(dat, aes(x = OriginZone, y = Count, fill = OriginZone)) +
    geom_boxplot(stat = "identity", width = 0.55) +
    scale_fill_manual(values = region_colors, na.value = "#BDBDBD") +
    labs(x = NULL, y = "Count per CloneFamily Tree") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  save_panel(p, paste0(output_tag, "_", sample_id, "_intergc_count.pdf"), 2.6, 3)
  save_table(dat, paste0(output_tag, "_", sample_id, "_intergc_count.csv"))
  invisible(dat)
}

plot_intergc_distance <- function(sample_id, edge_df = NULL, n_control = 1, seed = 1,
                                  output_tag = "figure3O") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  if (is.null(edge_df)) {
    edge_df <- read_tree_edges(sample_id)
  }

  intergc <- prepare_branch_edges(edge_df) %>%
    filter(GC_transition_type == "InterGC") %>%
    left_join(coords %>% select(spatial_barcode_ori, x_start = imagecol, y_start = imagerow),
              by = c("parent_spatial_barcode_ori" = "spatial_barcode_ori")) %>%
    left_join(coords %>% select(spatial_barcode_ori, x_end = imagecol, y_end = imagerow),
              by = c("node_spatial_barcode_ori" = "spatial_barcode_ori")) %>%
    mutate(Distance = sqrt((x_start - x_end)^2 + (y_start - y_end)^2), Group = "InterGC") %>%
    select(Distance, Group)

  set.seed(seed)
  control <- bind_rows(lapply(seq_len(max(1, nrow(intergc) * n_control)), function(i) {
    pair <- coords[sample(seq_len(nrow(coords)), 2), ]
    data.frame(
      Distance = sqrt(diff(pair$imagecol)^2 + diff(pair$imagerow)^2),
      Group = "Control"
    )
  }))

  dat <- bind_rows(intergc, control)
  p <- ggplot(dat, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.25) +
    scale_fill_manual(values = c("InterGC" = "#903586", "Control" = "#345D82")) +
    labs(x = NULL, y = "Spatial Distance") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")

  save_panel(p, paste0(output_tag, "_", sample_id, "_intergc_distance.pdf"), 2.4, 3)
  save_table(dat, paste0(output_tag, "_", sample_id, "_intergc_distance.csv"))
  invisible(dat)
}

# Examples:
# plot_intergc_origin_counts("Tonsil2")
# plot_intergc_distance("Tonsil2")
