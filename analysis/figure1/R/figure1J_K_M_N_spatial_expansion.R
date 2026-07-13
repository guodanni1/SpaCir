# Figure 1J/K/M/N: spatial clone expansion and nearest-neighbor distance.

source("analysis/figure1/R/00_setup.R")
suppressPackageStartupMessages({
  library(ggsignif)
  library(ggExtra)
})

calc_mnnd <- function(df_xy) {
  if (nrow(df_xy) <= 1) return(0)
  dmat <- as.matrix(dist(df_xy))
  diag(dmat) <- Inf
  mean(apply(dmat, 1, min))
}

classify_clone_size <- function(umi_per_clone) {
  cut(
    umi_per_clone,
    breaks = c(1, 2, 6, 20, Inf),
    labels = c("Not Expanded", "Small", "Medium", "Hyperexpanded"),
    right = FALSE
  )
}

calculate_clone_distance <- function(obj, sample_name, chain) {
  assay_name <- paste0("aaSeqCDR3_Lv1_UMI3_", chain, "_20G")
  mat <- as.matrix(obj@assays[[assay_name]]@data)
  nz <- which(mat > 0, arr.ind = TRUE)
  mapping <- data.frame(
    Clonotype = rownames(mat)[nz[, 1]],
    Spot = colnames(mat)[nz[, 2]],
    UMI = mat[nz]
  )

  img_key <- names(obj@images)[1]
  coords <- obj@images[[img_key]]@coordinates
  coords$barcode <- rownames(coords)

  mapping <- mapping |> inner_join(coords, by = c("Spot" = "barcode"))
  n_spots <- mapping |> distinct(Clonotype, Spot) |> count(Clonotype, name = "NumberOfCells")

  umi_per_clone <- rowSums(mat)
  dist_df <- mapping |>
    group_by(Clonotype) |>
    summarise(MeanNearestNeighborDistance = calc_mnnd(select(cur_data_all(), imagecol, imagerow)),
              .groups = "drop") |>
    left_join(n_spots, by = "Clonotype")

  dist_df$Classification <- classify_clone_size(umi_per_clone)[match(dist_df$Clonotype, names(umi_per_clone))]
  dist_df$UMI_counts <- umi_per_clone[match(dist_df$Clonotype, names(umi_per_clone))]
  dist_df$Sample <- sample_name
  dist_df$Chain <- chain
  dist_df
}

all_dist <- lapply(seq_len(nrow(sample_config)), function(i) {
  obj <- get(sample_config$object_name[i])
  bind_rows(
    calculate_clone_distance(obj, sample_config$display_sample[i], "IGH"),
    calculate_clone_distance(obj, sample_config$display_sample[i], "TRB")
  )
}) |> bind_rows()

write.csv(all_dist, file.path(output_dir, "figure1_clone_distance_table.csv"), row.names = FALSE)

plot_scatter <- function(chain, filename) {
  df <- all_dist |> filter(Chain == chain, MeanNearestNeighborDistance > 0)
  p <- ggplot(df, aes(x = NumberOfCells, y = MeanNearestNeighborDistance, color = Classification)) +
    geom_point(alpha = 0.95, size = 1.5) +
    scale_color_manual(values = c("Not Expanded" = "#99C945", "Small" = "orange",
                                  "Medium" = "blue", "Hyperexpanded" = "red")) +
    labs(title = paste(chain, "clonotypes"), x = "Number of Unique Spots",
         y = "Nearest Cell-to-Cell Distance") +
    theme_minimal() +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA))
  p_marg <- ggMarginal(p, type = "density", margins = "both", groupFill = TRUE, size = 6)
  save_panel(p_marg, filename, width = 6, height = 5)
}

plot_box <- function(chain, filename) {
  df <- all_dist |> filter(Chain == chain, Classification != "Not Expanded")
  p <- ggplot(df, aes(x = Classification, y = MeanNearestNeighborDistance, fill = Classification)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_signif(comparisons = list(c("Small", "Medium"), c("Medium", "Hyperexpanded")),
                map_signif_level = TRUE) +
    scale_fill_manual(values = c("Small" = "orange", "Medium" = "blue", "Hyperexpanded" = "red")) +
    labs(title = paste(chain, "clonotypes"), x = "", y = "Mean Nearest Neighbor Distance") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  save_panel(p, filename, width = 4, height = 4)
}

plot_scatter("IGH", "figure1J_IGH_clonotype_spatial_scatter.pdf")
plot_box("IGH", "figure1K_IGH_nearest_neighbor_distance.pdf")
plot_scatter("TRB", "figure1M_TRB_clonotype_spatial_scatter.pdf")
plot_box("TRB", "figure1N_TRB_nearest_neighbor_distance.pdf")
