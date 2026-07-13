# Figure 3P-R and supplementary TRB pathology panels.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

read_trb_cluster_table <- function(sample_id) {
  cfg <- sample_row(sample_id)
  read_csv_if_exists(file.path(repo_root, cfg$trb_cluster_table), stringsAsFactors = FALSE)
}

annotate_trb_pathology <- function(trb_df, pathology_db) {
  # Match by CDR3 first, and use V/J genes when they are available. The Word code
  # used VDJdb/McPAS-like tables; this function keeps that matching logic
  # configurable through pathology_database_csv.
  db <- pathology_db %>%
    rename_with(~ "CDR3b", any_of(c("CDR3", "CDR3.beta.aa", "aaSeqCDR3_most"))) %>%
    rename_with(~ "TRBV", any_of(c("bestVGene", "v_call"))) %>%
    rename_with(~ "TRBJ", any_of(c("bestJGene", "j_call")))

  trb_clean <- trb_df %>%
    rename_with(~ "CDR3b", any_of(c("aaSeqCDR3_most", "junction_aa"))) %>%
    rename_with(~ "TRBV", any_of(c("bestVGene", "v_call"))) %>%
    rename_with(~ "TRBJ", any_of(c("bestJGene", "j_call")))
  join_cols <- intersect(c("CDR3b", "TRBV", "TRBJ"), intersect(names(trb_clean), names(db)))
  trb_clean %>%
    left_join(db %>% select(any_of(c("CDR3b", "TRBV", "TRBJ", "Pathology"))) %>% distinct(),
              by = join_cols) %>%
    mutate(Pathology = ifelse(is.na(Pathology), "Unknown", Pathology))
}

plot_trb_cluster_distance <- function(sample_id, trb_df = NULL, n_control = 1, seed = 1,
                                      output_tag = "figure3P") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  if (is.null(trb_df)) {
    trb_df <- read_trb_cluster_table(sample_id)
  }
  cluster_col <- if ("Cluster" %in% names(trb_df)) "Cluster" else "clusterid_new"

  dat <- trb_df %>%
    filter(!is.na(.data[[cluster_col]])) %>%
    distinct(Cluster = .data[[cluster_col]], spatial_barcode_ori)
  distance_df <- nearest_neighbor_with_controls(dat, coords, group_col = "Cluster",
                                                n_control = n_control, seed = seed)

  p <- ggplot(distance_df, aes(x = Group, y = MeanNearestNeighborDistance, fill = Group)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.25) +
    scale_fill_manual(values = c("Cluster" = "#903586", "Control" = "#345D82")) +
    labs(x = NULL, y = "Mean Nearest Neighbor Distance") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")

  save_panel(p, paste0(output_tag, "_", sample_id, "_trb_cluster_distance.pdf"), 2.4, 3)
  save_table(distance_df, paste0(output_tag, "_", sample_id, "_trb_cluster_distance.csv"))
  invisible(distance_df)
}

plot_pathology_bar <- function(sample_id, annotated_trb_df, output_tag = "figure3Q") {
  dat <- annotated_trb_df %>%
    filter(Pathology != "Unknown", Pathology != "HomoSapiens") %>%
    count(Pathology, name = "Count") %>%
    mutate(Percentage = 100 * Count / sum(Count))

  p <- ggplot(dat, aes(x = "", y = Percentage, fill = reorder(Pathology, -Count))) +
    geom_col(width = 0.65) +
    labs(x = NULL, y = "TRB clusters Pathology annotation Percentage (%)", fill = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  save_panel(p, paste0(output_tag, "_", sample_id, "_trb_pathology_bar.pdf"), 3, 3)
  save_table(dat, paste0(output_tag, "_", sample_id, "_trb_pathology_bar.csv"))
  invisible(dat)
}

plot_pathology_spatial <- function(sample_id, annotated_trb_df, selected_pathologies = NULL,
                                   output_tag = "figure3R") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  coords <- get_spatial_coordinates(obj, cfg$image_name)
  cluster_col <- if ("Cluster" %in% names(annotated_trb_df)) "Cluster" else "clusterid_new"

  dat <- annotated_trb_df %>%
    filter(Pathology != "Unknown") %>%
    { if (!is.null(selected_pathologies)) filter(., Pathology %in% selected_pathologies) else . } %>%
    distinct(spatial_barcode_ori, Pathology, Cluster = .data[[cluster_col]]) %>%
    inner_join(coords, by = "spatial_barcode_ori")

  p <- ggplot(coords, aes(x = imagecol, y = imagerow)) +
    geom_point(color = "grey88", size = 0.2) +
    geom_point(data = dat, aes(fill = Cluster), shape = 21, size = 0.6, color = "white", linewidth = 0.1) +
    facet_wrap(~ Pathology) +
    scale_fill_viridis_c(option = "turbo") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void(base_size = 7) +
    theme(legend.position = "bottom")

  save_panel(p, paste0(output_tag, "_", sample_id, "_trb_pathology_spatial.pdf"), 6, 4)
  invisible(dat)
}

run_trb_pathology_panels <- function(sample_id, output_prefix = "figure3") {
  cfg <- sample_row(sample_id)
  trb_df <- read_trb_cluster_table(sample_id)
  pathology_db <- read_csv_if_exists(file.path(repo_root, cfg$pathology_database_csv), stringsAsFactors = FALSE)
  annotated <- annotate_trb_pathology(trb_df, pathology_db)
  save_table(annotated, paste0(output_prefix, "_", sample_id, "_trb_pathology_annotated.csv"))
  plot_trb_cluster_distance(sample_id, annotated, output_tag = paste0(output_prefix, "P"))
  plot_pathology_bar(sample_id, annotated, output_tag = paste0(output_prefix, "Q"))
  plot_pathology_spatial(sample_id, annotated, output_tag = paste0(output_prefix, "R"))
  invisible(annotated)
}

# Example:
# run_trb_pathology_panels("Tonsil2")
