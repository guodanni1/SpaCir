# Figure 3H-I: spatial clonefamily maps and branch-transition arrows.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

plot_clonefamily_spatial_map <- function(sample_id, clonefamily_df = NULL,
                                         clonefamily_col = "Clonefamily",
                                         output_tag = "figure3H") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  if (is.null(clonefamily_df)) {
    clonefamily_df <- read_csv_if_exists(file.path(repo_root, cfg$clonefamily_table), stringsAsFactors = FALSE)
  }

  dat <- clonefamily_df %>%
    filter(!is.na(.data[[clonefamily_col]]), .data[[clonefamily_col]] != "None") %>%
    distinct(spatial_barcode_ori, Clonefamily = .data[[clonefamily_col]]) %>%
    inner_join(coords, by = "spatial_barcode_ori")

  p <- ggplot(coords, aes(x = imagecol, y = imagerow)) +
    geom_point(color = "grey85", size = 0.25) +
    geom_point(data = dat, aes(color = Clonefamily), size = 0.55) +
    scale_y_reverse() +
    coord_fixed() +
    labs(color = "Clonefamily") +
    theme_void(base_size = 7)

  save_panel(p, paste0(output_tag, "_", sample_id, "_clonefamily_spatial.pdf"), 4, 4)
  invisible(dat)
}

read_tree_edges <- function(sample_id) {
  cfg <- sample_row(sample_id)
  read_csv_if_exists(file.path(repo_root, cfg$tree_edges_csv), stringsAsFactors = FALSE)
}

prepare_branch_edges <- function(edge_df, parent_region_col = "parent_GC_NewClusters",
                                 node_region_col = "node_GC_NewClusters") {
  # Each edge represents a parent-node relationship in the clonal tree. The
  # transition class follows the original Figure 3 logic: InnerGC, InterGC,
  # Out GC, To GC, or Non-GC.
  edge_df %>%
    mutate(
      GC_transition_type = classify_gc_transition(.data[[parent_region_col]], .data[[node_region_col]]),
      branch_group = coalesce(.data$branch_group, paste0(.data$Clonefamily, "_Branch_", .data$branch))
    )
}

filter_branch_edges <- function(edge_df, transition = NULL, parent_pattern = NULL, node_pattern = NULL) {
  out <- edge_df
  if (!is.null(transition)) {
    out <- out %>% filter(GC_transition_type %in% transition)
  }
  if (!is.null(parent_pattern)) {
    out <- out %>% filter(grepl(parent_pattern, parent_GC_NewClusters))
  }
  if (!is.null(node_pattern)) {
    out <- out %>% filter(grepl(node_pattern, node_GC_NewClusters))
  }
  out
}

plot_branch_arrows <- function(sample_id, edge_df = NULL, output_tag = "figure3I",
                               transition = NULL, parent_pattern = "Follicular LZ$",
                               node_pattern = "Inter LZ-DZ$") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  if (is.null(edge_df)) {
    edge_df <- read_tree_edges(sample_id)
  }
  edges <- prepare_branch_edges(edge_df) %>%
    filter_branch_edges(transition = transition, parent_pattern = parent_pattern, node_pattern = node_pattern)

  if (!all(c("parent_spatial_barcode_ori", "node_spatial_barcode_ori") %in% names(edges))) {
    stop("Tree edge table must contain parent_spatial_barcode_ori and node_spatial_barcode_ori.", call. = FALSE)
  }

  arrow_df <- edges %>%
    left_join(coords %>% select(spatial_barcode_ori, x_start = imagecol, y_start = imagerow),
              by = c("parent_spatial_barcode_ori" = "spatial_barcode_ori")) %>%
    left_join(coords %>% select(spatial_barcode_ori, x_end = imagecol, y_end = imagerow),
              by = c("node_spatial_barcode_ori" = "spatial_barcode_ori"))

  p <- ggplot(coords, aes(x = imagecol, y = imagerow)) +
    geom_point(color = "grey90", size = 0.2) +
    geom_curve(
      data = arrow_df,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = branch_group),
      curvature = 0.2,
      arrow = arrow(length = unit(1.5, "mm")),
      linewidth = 0.35,
      alpha = 0.9
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(color = "Branch") +
    theme_void(base_size = 7) +
    theme(legend.position = "right")

  save_panel(p, paste0(output_tag, "_", sample_id, "_branch_arrows.pdf"), 6, 4)
  save_table(arrow_df, paste0(output_tag, "_", sample_id, "_branch_arrows.csv"))
  invisible(arrow_df)
}

# Examples:
# plot_clonefamily_spatial_map("Tonsil2")
# plot_branch_arrows("Tonsil2")
