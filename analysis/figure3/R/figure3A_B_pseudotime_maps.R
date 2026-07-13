# Figure 3A-B: Monocle3 pseudotime on UMAP and tissue coordinates.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(Matrix)
  library(patchwork)
})

build_monocle_cds <- function(seurat_obj, assay = "Spatial", reduction = "UMAP") {
  expression_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
  cell_metadata <- seurat_obj@meta.data
  gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

  cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
  reducedDims(cds)$UMAP <- Embeddings(seurat_obj, reduction = reduction)[colnames(cds), ]
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)
  cds
}

add_pseudotime_to_seurat <- function(seurat_obj, cds, column = "Pseudotime") {
  pt <- pseudotime(cds)
  seurat_obj[[column]] <- NA_real_
  seurat_obj[[column]][names(pt), 1] <- pt
  seurat_obj
}

plot_pseudotime_maps <- function(sample_id, region_filter = NULL, pseudotime_column = "Pseudotime",
                                 output_tag = "figure3A") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  if (!is.null(region_filter)) {
    keep <- rownames(obj@meta.data)[obj@meta.data[[cfg$gc_region_column]] %in% region_filter]
    obj <- subset(obj, cells = keep)
  }

  cds <- build_monocle_cds(obj)
  obj <- add_pseudotime_to_seurat(obj, cds, pseudotime_column)

  # UMAP displays the learned trajectory; SpatialFeaturePlot checks that the
  # same pseudotime trend is anatomically coherent in the tissue.
  p_umap <- plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "pseudotime") +
    ggtitle(sample_id)
  p_spatial <- SpatialFeaturePlot(obj, features = pseudotime_column, image.alpha = 0.35) +
    ggtitle(sample_id)

  save_panel(p_umap + p_spatial, paste0(output_tag, "_", sample_id, "_pseudotime.pdf"), 7.4, 3.7)
  saveRDS(obj, file.path(output_dir, paste0(output_tag, "_", sample_id, "_seurat_with_pseudotime.rds")))
  invisible(list(seurat = obj, cds = cds))
}

# Figure 3A:
# plot_pseudotime_maps("LN2", output_tag = "figure3A")
#
# Figure 3B:
# plot_pseudotime_maps("LN2", region_filter = c("Follicular DZ", "Follicular LZ", "Inter LZ-DZ"),
#                      pseudotime_column = "Pseudotime_GC", output_tag = "figure3B")
