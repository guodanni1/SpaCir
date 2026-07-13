# Figure 5C/D, Figure S11A, and Figure S12C/D: SHM-associated gene programs.
#
# Genes and functional categories are declared in config/shm_gene_sets.csv.
# The AUC panel uses Seurat FindAllMarkers(test.use = "roc") on SHM_VDJ groups.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

plot_auc_lollipop <- function(roc_table, genes = ordered_shm_genes) {
  auc_df <- roc_table %>%
    filter(.data$gene %in% genes) %>%
    group_by(.data$gene) %>%
    slice_max(order_by = .data$myAUC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      gene = factor(.data$gene, levels = rev(genes)),
      set_name = shm_gene_table$set_name[match(as.character(.data$gene), shm_gene_table$gene)]
    )
  ggplot(auc_df, aes(x = .data$myAUC, y = .data$gene, color = .data$set_name)) +
    geom_segment(aes(x = 0, xend = .data$myAUC, y = .data$gene, yend = .data$gene), linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("High_Mutation" = "#E3170D", "Low_Mutation" = "#8B58A4")) +
    labs(x = "AUC", y = NULL, color = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_text(face = "italic"))
}

run_shm_gene_program_panels <- function(sample_id = "Tonsil2", output_tag = "figure5") {
  obj_info <- load_sample_object(sample_id, filtered = TRUE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  shm_meta <- prepare_spot_mutation_metadata(seurat_obj, cfg)
  seurat_obj <- attach_shm_group(seurat_obj, shm_meta)
  seurat_obj <- subset(seurat_obj, subset = SHM_VDJ %in% shm_levels)
  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  Idents(seurat_obj) <- "SHM_VDJ"

  roc_table <- FindAllMarkers(
    object = seurat_obj,
    test.use = "roc",
    group.by = "SHM_VDJ",
    only.pos = FALSE,
    return.thresh = 0,
    min.pct = 0,
    logfc.threshold = 0
  ) %>%
    filter(!grepl("^MT-", .data$gene))

  p_auc <- plot_auc_lollipop(roc_table)
  save_panel(p_auc, paste0(output_tag, "C_", cfg$output_prefix, "_shm_gene_auc.pdf"), width = 3.6, height = 5.2)

  genes_present <- ordered_shm_genes[ordered_shm_genes %in% rownames(seurat_obj)]
  p_dot <- DotPlot(seurat_obj, features = rev(genes_present), group.by = "SHM_VDJ") +
    coord_flip() +
    scale_color_gradient2(high = "#A52A2A", mid = "grey95", low = "#191970") +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_text(face = "italic", size = 8))
  save_panel(p_dot, paste0(output_tag, "D_", cfg$output_prefix, "_shm_gene_dotplot.pdf"), width = 4.4, height = 5.2)

  save_table(roc_table, paste0(output_tag, "_", cfg$output_prefix, "_shm_gene_roc_table.csv"))
  invisible(list(object = seurat_obj, shm_meta = shm_meta, roc = roc_table))
}

run_s11_celltype_gene_dotplot <- function(sample_id = "Tonsil2") {
  obj_info <- load_sample_object(sample_id, filtered = FALSE)
  cfg <- obj_info$cfg
  seurat_obj <- obj_info$object
  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  group_col <- "celltype"
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("Expected a celltype metadata column for Figure S11A.", call. = FALSE)
  }
  genes_present <- ordered_shm_genes[ordered_shm_genes %in% rownames(seurat_obj)]
  p <- DotPlot(seurat_obj, features = rev(genes_present), group.by = group_col) +
    coord_flip() +
    scale_color_gradient2(high = "#A52A2A", mid = "grey95", low = "#191970") +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(face = "italic", size = 8))
  save_panel(p, paste0("figureS11A_", cfg$output_prefix, "_shm_gene_celltype_dotplot.pdf"), width = 8.8, height = 5.8)
  invisible(p)
}

# Examples:
# run_shm_gene_program_panels("Tonsil2", output_tag = "figure5")
# run_shm_gene_program_panels("LN2", output_tag = "figureS12")
# run_s11_celltype_gene_dotplot("Tonsil2")
