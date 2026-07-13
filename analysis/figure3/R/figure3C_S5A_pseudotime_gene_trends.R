# Figure 3C and Figure S5A: gene expression trends along pseudotime.
#
# Monocle3 and gene-module settings are recorded in:
#   config/sample_specific_parameters.csv
# The plotted gene lists are recorded in:
#   config/pseudotime_gene_sets.csv

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

figure3_genes <- c(
  "AICDA", "JPT1", "DUT", "HMGN2",
  "BATF", "HLA-A", "FDCSP", "CXCL13",
  "MZB1", "XBP1", "PTGDS", "HSP90B1",
  "SSR4", "IGHA1", "IGHG1", "IGKC"
)

figureS5_genes <- c(
  "MCM3", "PCNA", "TYMS", "TOP2A", "PARP1", "BIRC5", "CCNB2", "TK1",
  "ATP5MG", "H2AFZ", "IL32", "CD3D", "CD83", "SRGN", "MS4A1", "XBP1",
  "CLU", "DERL3", "UBE2J1", "CR2", "C1S", "CD63", "SEC11C", "FKBP11"
)

plot_pseudotime_gene_trends <- function(sample_id, seurat_rds_with_pseudotime,
                                        genes = figure3_genes,
                                        pseudotime_column = "Pseudotime_GC",
                                        region_column = "GC_NewClusters",
                                        output_tag = "figure3C") {
  obj <- readRDS(seurat_rds_with_pseudotime)
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) == 0) {
    stop("No requested genes were found in the Seurat object.", call. = FALSE)
  }

  dat <- FetchData(obj, vars = c(genes, pseudotime_column, region_column)) %>%
    mutate(Spot = rownames(.)) %>%
    pivot_longer(all_of(genes), names_to = "Gene", values_to = "Expression") %>%
    filter(!is.na(.data[[pseudotime_column]]))

  # A loess line matches the original trend-style panels while keeping the raw
  # spot-level points visible.
  p <- ggplot(dat, aes(x = .data[[pseudotime_column]], y = Expression, color = .data[[region_column]])) +
    geom_point(size = 0.35, alpha = 0.75) +
    geom_smooth(aes(group = Gene), method = "loess", se = FALSE, color = "black", linewidth = 0.3) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = region_colors, na.value = "#BDBDBD") +
    labs(x = "Pseudotime", y = "Expression", color = NULL) +
    theme_classic(base_size = 7) +
    theme(legend.position = "bottom")

  save_panel(p, paste0(output_tag, "_", sample_id, "_pseudotime_gene_trends.pdf"), 8, 6)
  invisible(dat)
}

# Examples:
# plot_pseudotime_gene_trends("LN2", "analysis/figure3/outputs/figure3B_LN2_seurat_with_pseudotime.rds")
# plot_pseudotime_gene_trends("LN2", "analysis/figure3/outputs/figure3B_LN2_seurat_with_pseudotime.rds",
#                             genes = figureS5_genes, output_tag = "figureS5A")
