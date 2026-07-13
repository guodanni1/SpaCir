# Figure 3D-G: IGH clonefamily size, CSR, CDR3 length, and spatial distance.
#
# IGH clonefamily thresholds, per-sample clonefamily counts, CSR nodes, and
# random-control rules are recorded in:
#   config/clonefamily_parameters.csv

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggraph)
  library(igraph)
})

read_clonefamily_table <- function(sample_id) {
  cfg <- sample_row(sample_id)
  read_csv_if_exists(file.path(repo_root, cfg$clonefamily_table), stringsAsFactors = FALSE)
}

plot_clonefamily_size <- function(sample_id, clonefamily_df = read_clonefamily_table(sample_id)) {
  size_df <- clonefamily_df %>%
    filter(!is.na(Clonefamily), Clonefamily != "None") %>%
    count(Clonefamily, name = "CloneSize")

  p <- ggplot(size_df, aes(x = Clonefamily, y = CloneSize, fill = CloneSize)) +
    geom_col(width = 0.8) +
    scale_fill_viridis_c(option = "turbo") +
    labs(x = "Clonefamily", y = "Clonesize", fill = "Clonesize") +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  save_panel(p, paste0("figure3D_", sample_id, "_clonefamily_size.pdf"), 3, 3)
  save_table(size_df, paste0("figure3D_", sample_id, "_clonefamily_size.csv"))
  invisible(size_df)
}

plot_clonefamily_csr <- function(sample_id, clonefamily_df = read_clonefamily_table(sample_id)) {
  edge_df <- clonefamily_df %>%
    filter(!is.na(Clonefamily), Clonefamily != "None", !is.na(c_call)) %>%
    distinct(Clonefamily, junction_aa, c_call) %>%
    group_by(Clonefamily) %>%
    arrange(c_call, .by_group = TRUE) %>%
    summarise(CGenes = list(unique(c_call)), .groups = "drop") %>%
    mutate(Edges = map(CGenes, ~ if (length(.x) > 1) data.frame(from = head(.x, -1), to = tail(.x, -1)) else data.frame())) %>%
    select(Clonefamily, Edges) %>%
    unnest(Edges) %>%
    count(from, to, name = "weight")

  if (nrow(edge_df) == 0) {
    message("No CSR edges for ", sample_id)
    return(invisible(edge_df))
  }

  graph <- graph_from_data_frame(edge_df, directed = TRUE)
  p <- ggraph(graph, layout = "circle") +
    geom_edge_link(aes(width = weight), alpha = 0.45, arrow = arrow(length = unit(2, "mm"))) +
    geom_node_point(size = 3, color = "#F2B33D") +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
    theme_void()

  save_panel(p, paste0("figure3E_", sample_id, "_clonefamily_csr.pdf"), 3, 3)
  save_table(edge_df, paste0("figure3E_", sample_id, "_clonefamily_csr_edges.csv"))
  invisible(edge_df)
}

plot_cdr3_length_distribution <- function(sample_id, clonefamily_df = read_clonefamily_table(sample_id)) {
  dat <- clonefamily_df %>%
    filter(!is.na(junction_aa)) %>%
    mutate(
      CDR3Length = nchar(junction_aa),
      ClusterStatus = ifelse(is.na(Clonefamily) | Clonefamily == "None", "Unclustered", "Clustered")
    )

  p <- ggplot(dat, aes(x = CDR3Length, fill = ClusterStatus, color = ClusterStatus)) +
    geom_density(alpha = 0.35) +
    scale_fill_manual(values = c("Clustered" = "#A5A832", "Unclustered" = "#63BE7B")) +
    scale_color_manual(values = c("Clustered" = "#A5A832", "Unclustered" = "#63BE7B")) +
    labs(x = "CDR3 length", y = NULL, fill = NULL, color = NULL) +
    theme_classic(base_size = 7)

  save_panel(p, paste0("figure3F_", sample_id, "_cdr3_length.pdf"), 3, 2.5)
  invisible(dat)
}

plot_clonefamily_distance <- function(sample_id, clonefamily_df = read_clonefamily_table(sample_id),
                                      n_control = 1, seed = 1, output_tag = "figure3G") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  dat <- clonefamily_df %>%
    filter(!is.na(Clonefamily), Clonefamily != "None") %>%
    distinct(Clonefamily, spatial_barcode_ori)

  distance_df <- nearest_neighbor_with_controls(dat, coords, group_col = "Clonefamily",
                                                spot_col = "spatial_barcode_ori",
                                                n_control = n_control, seed = seed)

  p <- ggplot(distance_df, aes(x = Group, y = MeanNearestNeighborDistance, fill = Group)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.25) +
    scale_fill_manual(values = c("Cluster" = "#903586", "Control" = "#345D82")) +
    labs(x = NULL, y = "Mean Nearest Neighbor Distance") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")

  save_panel(p, paste0(output_tag, "_", sample_id, "_clonefamily_distance.pdf"), 2.4, 3)
  save_table(distance_df, paste0(output_tag, "_", sample_id, "_clonefamily_distance.csv"))
  invisible(distance_df)
}

# Examples:
# plot_clonefamily_size("LN2")
# plot_clonefamily_csr("LN2")
# plot_cdr3_length_distribution("LN2")
# plot_clonefamily_distance("LN2")
