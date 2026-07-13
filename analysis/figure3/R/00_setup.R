# Shared settings and helper functions for SpaCir Figure 3 and Figure S5-S7.
# Run scripts from the repository root, or set SPACIR_ROOT before sourcing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(purrr)
})

repo_root <- Sys.getenv("SPACIR_ROOT", unset = getwd())
figure3_dir <- file.path(repo_root, "analysis", "figure3")
input_config <- file.path(figure3_dir, "config", "figure3_inputs.csv")
output_dir <- file.path(figure3_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_config <- read.csv(input_config, stringsAsFactors = FALSE)

sample_levels <- c("LN1", "LN2", "Tonsil1", "Tonsil2")
sample_colors <- c(
  "LN1" = "#CB5D17",
  "LN2" = "#0C6CA7",
  "Tonsil1" = "#E9CC54",
  "Tonsil2" = "#823379"
)

region_colors <- c(
  "Follicular DZ" = "#E4572E",
  "Follicular LZ" = "#E9C94F",
  "Inter LZ-DZ" = "#9BCB66",
  "Mantle zone" = "#F29B7A",
  "Inter-Follicular CD5+" = "#7C66A2",
  "Inter-Follicular CD5-" = "#6FC7B6",
  "Inter-Follicular SDC1+" = "#68BFAE",
  "Inter-Follicular low SDC1+" = "#B8D989",
  "Inter-Follicular IL7R+" = "#BCD58A",
  "Inter-Follicular CXCL10+" = "#8A6FA6",
  "Inter-Follicular CD6+" = "#7B669A",
  "Inter-Follicular XBP1+" = "#8DD3C7",
  "Muscle/Mesenchymal" = "#F0B5D1",
  "Vascular" = "#F3E999",
  "Vascular/Mesenchymal" = "#C46693",
  "Subepithelial/Mesenchymal/Vascular" = "#6C6E70",
  "Epithelium-like structures" = "#F0B3C7",
  "Unannotated regions" = "#9A9A9A"
)

gc_colors <- c(
  "GC1" = "#8DD3C7", "GC2" = "#FFFFB3", "GC3" = "#BEBADA",
  "GC4" = "#FB8072", "GC5" = "#80B1D3", "GC6" = "#FDB462",
  "GC7" = "#B3DE69", "GC8" = "#FCCDE5", "GC9" = "#D9D9D9",
  "GC10" = "#BC80BD", "GC11" = "#CCEBC5", "GC12" = "#FFED6F",
  "GC13" = "#66C2A5"
)

transition_colors <- c(
  "InnerGC" = "#A6CEE3",
  "InterGC" = "#1086BB",
  "Out GC" = "#FB9A99",
  "To GC" = "#CAB2D6",
  "Non-GC" = "#FDBF6F"
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

sample_row <- function(sample_id) {
  row <- sample_config[sample_config$sample_id == sample_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one row in figure3_inputs.csv for sample: ", sample_id, call. = FALSE)
  }
  row
}

read_csv_if_exists <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  read.csv(path, ...)
}

save_panel <- function(plot, filename, width, height) {
  out <- file.path(output_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved: ", out)
}

save_table <- function(df, filename) {
  out <- file.path(output_dir, filename)
  write.csv(df, out, row.names = FALSE)
  message("Saved: ", out)
}

get_spatial_coordinates <- function(seurat_obj, image_name = NULL) {
  if (is.null(image_name)) {
    image_name <- names(seurat_obj@images)[1]
  }
  coords <- seurat_obj@images[[image_name]]@coordinates
  coords$spatial_barcode_ori <- rownames(coords)
  coords
}

mean_nearest_neighbor <- function(coords) {
  coords <- coords[, c("imagecol", "imagerow")]
  coords <- coords[complete.cases(coords), , drop = FALSE]
  if (nrow(coords) <= 1) {
    return(0)
  }
  d <- as.matrix(dist(coords))
  diag(d) <- Inf
  mean(apply(d, 1, min))
}

nearest_neighbor_with_controls <- function(feature_df, coords_df, group_col, spot_col = "spatial_barcode_ori",
                                           n_control = 1, seed = 1) {
  set.seed(seed)
  feature_coords <- feature_df %>%
    inner_join(coords_df, by = setNames("spatial_barcode_ori", spot_col))
  all_spots <- unique(coords_df$spatial_barcode_ori)

  observed <- feature_coords %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      MeanNearestNeighborDistance = mean_nearest_neighbor(cur_data()),
      NumberOfSpots = n_distinct(.data[[spot_col]]),
      Group = "Cluster",
      .groups = "drop"
    ) %>%
    rename(Feature = all_of(group_col))

  controls <- observed %>%
    group_by(Feature) %>%
    group_modify(function(.x, .y) {
      n <- .x$NumberOfSpots[1]
      bind_rows(lapply(seq_len(n_control), function(i) {
        sampled <- sample(all_spots, n, replace = FALSE)
        sampled_coords <- coords_df %>% filter(spatial_barcode_ori %in% sampled)
        data.frame(
          Feature = .y$Feature,
          MeanNearestNeighborDistance = mean_nearest_neighbor(sampled_coords),
          NumberOfSpots = n,
          Group = "Control"
        )
      }))
    }) %>%
    ungroup()

  bind_rows(observed, controls)
}

jaccard_index <- function(x, y) {
  union_n <- sum(x > 0 | y > 0)
  if (union_n == 0) {
    return(0)
  }
  sum(x > 0 & y > 0) / union_n
}

presence_by_group <- function(assay_matrix, spot_to_group) {
  long <- as.data.frame(as.table(as.matrix(assay_matrix)), stringsAsFactors = FALSE)
  colnames(long) <- c("Clonotype", "Spot", "Value")
  long %>%
    filter(Value > 0) %>%
    inner_join(spot_to_group, by = "Spot") %>%
    distinct(Group, Clonotype) %>%
    mutate(Present = 1) %>%
    pivot_wider(names_from = Group, values_from = Present, values_fill = 0)
}

compute_shared_clone_matrix <- function(presence_df) {
  mat <- as.matrix(presence_df[, setdiff(names(presence_df), "Clonotype"), drop = FALSE])
  rownames(mat) <- presence_df$Clonotype
  shared <- t(mat) %*% mat
  shared
}

compute_jaccard_matrix <- function(presence_df) {
  mat <- as.matrix(presence_df[, setdiff(names(presence_df), "Clonotype"), drop = FALSE])
  groups <- colnames(mat)
  out <- matrix(0, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) {
      out[i, j] <- jaccard_index(mat[, i], mat[, j])
    }
  }
  out
}

classify_gc_transition <- function(parent_region, node_region) {
  parent_is_gc <- grepl("GC", parent_region)
  node_is_gc <- grepl("GC", node_region)
  parent_gc <- str_extract(parent_region, "GC\\d+")
  node_gc <- str_extract(node_region, "GC\\d+")
  case_when(
    parent_is_gc & node_is_gc & parent_gc == node_gc ~ "InnerGC",
    parent_is_gc & node_is_gc & parent_gc != node_gc ~ "InterGC",
    parent_is_gc & !node_is_gc ~ "Out GC",
    !parent_is_gc & node_is_gc ~ "To GC",
    TRUE ~ "Non-GC"
  )
}

collapse_zone <- function(region) {
  case_when(
    grepl("Follicular LZ", region) ~ "Follicular LZ",
    grepl("Inter LZ-DZ", region) ~ "Inter LZ-DZ",
    grepl("Follicular DZ", region) ~ "Follicular DZ",
    TRUE ~ "Other"
  )
}
