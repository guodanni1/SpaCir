# Figure 1C: correlation between spatial gene expression and VDJ UMI counts.

source("analysis/figure1/R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat)
})

merge_summary_metadata <- function(obj, summary_csv, prefix) {
  df <- read_csv_if_exists(summary_csv, stringsAsFactors = FALSE)
  colnames(df)[1] <- "spatial_bc"
  all_spots <- rownames(obj@meta.data)
  full <- data.frame(spatial_bc = all_spots) |>
    left_join(df, by = "spatial_bc")
  full[is.na(full)] <- 0
  rownames(full) <- full$spatial_bc
  full$spatial_bc <- NULL
  colnames(full) <- paste(prefix, colnames(full), sep = "_")
  full
}

average_expression <- function(obj, genes) {
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0) stop("No requested genes found in Seurat object.", call. = FALSE)
  Matrix::colMeans(obj@assays$Spatial@data[present, , drop = FALSE])
}

plot_correlation <- function(plot_data, x_col, y_col, title, x_lab, y_lab, filename) {
  plot_data$Sample <- factor(plot_data$Sample, levels = sample_levels)
  p <- ggplot(plot_data, aes(x = .data[[x_col]], y = .data[[y_col]], color = Sample)) +
    geom_point(size = 1.2, alpha = 0.55) +
    geom_smooth(method = "lm", se = TRUE, size = 0.8) +
    scale_color_manual(values = sample_colors, drop = FALSE) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_classic()
  save_panel(p, filename, width = 5, height = 4)
}

trbc_genes <- c("TRBC1", "TRBC2")
igh_genes <- c("IGHEP2", "IGHMBP2", "IGHA2", "IGHE", "IGHG4",
               "IGHG2", "IGHGP", "IGHA1", "IGHEP1", "IGHG1",
               "IGHG3", "IGHD", "IGHM")

cor_data <- lapply(seq_len(nrow(sample_config)), function(i) {
  object_name <- sample_config$object_name[i]
  obj <- get(object_name)
  label <- sample_config$display_sample[i]

  trb_meta <- merge_summary_metadata(obj, sample_config$trb_summary_csv[i], "TRB")
  igh_meta <- merge_summary_metadata(obj, sample_config$igh_summary_csv[i], "IGH")

  trb_umi_col <- grep("TRB_.*UMI_count|TRB_.*match_count", colnames(trb_meta), value = TRUE)[1]
  igh_umi_col <- grep("IGH_.*UMI_count|IGH_.*match_count", colnames(igh_meta), value = TRUE)[1]

  data.frame(
    Sample = label,
    TRBC_expression = average_expression(obj, trbc_genes),
    IGH_expression = average_expression(obj, igh_genes),
    TRB_UMI = as.numeric(trb_meta[[trb_umi_col]]),
    IGH_UMI = as.numeric(igh_meta[[igh_umi_col]])
  )
}) |> bind_rows()

write.csv(cor_data, file.path(output_dir, "figure1C_correlation_input_table.csv"), row.names = FALSE)

plot_correlation(
  cor_data,
  "TRB_UMI",
  "TRBC_expression",
  "Correlation of TRBC Gene Expression vs Spatial-TRB UMI count",
  "Spatial-TRB UMI counts",
  "TRBC Gene UMI counts",
  "figure1C_TRBC_expression_vs_TRB_UMI.pdf"
)

plot_correlation(
  cor_data,
  "IGH_UMI",
  "IGH_expression",
  "Correlation of IGH Gene Expression vs Spatial-IGH UMI count",
  "Spatial-IGH UMI counts",
  "IGHC Gene UMI counts",
  "figure1C_IGH_expression_vs_IGH_UMI.pdf"
)
