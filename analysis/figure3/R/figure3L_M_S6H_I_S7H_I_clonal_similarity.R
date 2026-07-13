# Figure 3L/M and supplementary clonal similarity panels.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggdendro)
})

clonotype_presence_from_seurat <- function(sample_id, assay_name, group_column) {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  mat <- as.matrix(obj@assays[[assay_name]]@data)
  spot_to_group <- obj@meta.data %>%
    mutate(Spot = rownames(.), Group = .data[[group_column]]) %>%
    select(Spot, Group) %>%
    filter(!is.na(Group))
  presence_by_group(mat, spot_to_group)
}

plot_jaccard_tree <- function(sample_id, presence_df, output_tag = "figure3L") {
  jac <- compute_jaccard_matrix(presence_df)
  dist_matrix <- as.dist(1 - jac)
  hc <- hclust(dist_matrix)
  dendro <- ggdendro::dendro_data(hc, type = "rectangle")

  p <- ggplot() +
    geom_segment(data = dendro$segments, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.25) +
    coord_flip() +
    scale_x_continuous(breaks = seq_along(dendro$labels$label), labels = dendro$labels$label) +
    labs(x = NULL, y = "Jaccard Distance", title = "BCR Clonal Similarity Tree") +
    theme_classic(base_size = 7)

  save_panel(p, paste0(output_tag, "_", sample_id, "_jaccard_tree.pdf"), 4, 5)
  invisible(jac)
}

plot_shared_clone_bubble <- function(sample_id, presence_df, output_tag = "figure3M") {
  shared <- compute_shared_clone_matrix(presence_df)
  shared[upper.tri(shared)] <- NA
  dat <- as.data.frame(as.table(shared), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    rename(Group1 = Var1, Group2 = Var2, SharedClonotypes = Freq)

  p <- ggplot(dat, aes(x = Group1, y = Group2, size = SharedClonotypes, fill = SharedClonotypes)) +
    geom_point(shape = 21, color = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "#E5F5E0", high = "#2B8CBE", name = "Shared Clones") +
    scale_size_continuous(name = "Shared Clones") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_panel(p, paste0(output_tag, "_", sample_id, "_shared_clone_bubble.pdf"), 5, 4)
  save_table(dat, paste0(output_tag, "_", sample_id, "_shared_clone_counts.csv"))
  invisible(dat)
}

# Example:
# presence <- clonotype_presence_from_seurat("Tonsil2", "aaSeqCDR3_Lv1_UMI3_IGH_20G", "GC_all")
# plot_jaccard_tree("Tonsil2", presence)
# plot_shared_clone_bubble("Tonsil2", presence)
