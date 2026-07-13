# Figure 1F/G: TRB J gene usage and IGH J gene usage by sample.

source("analysis/figure1/R/00_setup.R")

prepare_gene_usage <- function(data, sample_name, gene_col, top_n = NULL) {
  gene_count <- table(data[[gene_col]])
  gene_data <- as.data.frame(gene_count)
  colnames(gene_data) <- c("Gene", "Count")
  gene_data <- gene_data |>
    filter(!is.na(Gene), Gene != "", Gene != "nan") |>
    mutate(
      Percentage = Count / sum(Count) * 100,
      Sample = sample_name
    ) |>
    arrange(Sample, desc(Percentage))

  if (!is.null(top_n)) {
    keep <- gene_data |> slice_head(n = top_n) |> pull(Gene)
    gene_data <- gene_data |>
      mutate(Gene = ifelse(Gene %in% keep, as.character(Gene), "Others")) |>
      group_by(Sample, Gene) |>
      summarise(Count = sum(Count), Percentage = sum(Percentage), .groups = "drop")
  }
  gene_data
}

plot_gene_usage <- function(all_data, title, fill_title, filename, palette) {
  all_data$Sample <- factor(all_data$Sample, levels = sample_levels)
  p <- ggplot(all_data, aes(x = Sample, y = Percentage, fill = Gene)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = palette) +
    labs(title = title, x = "Sample", y = "Percentage (%)", fill = fill_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_blank())
  save_panel(p, filename, width = 4.2, height = 3.2)
}

trb_data <- lapply(seq_len(nrow(sample_config)), function(i) {
  df <- read_csv_if_exists(sample_config$trb_expanded_csv[i]) |> distinct()
  prepare_gene_usage(df, sample_config$display_sample[i], "bestJGene")
}) |> bind_rows()

igh_data <- lapply(seq_len(nrow(sample_config)), function(i) {
  df <- read_csv_if_exists(sample_config$igh_expanded_csv[i]) |> distinct()
  prepare_gene_usage(df, sample_config$display_sample[i], "bestJGene")
}) |> bind_rows()

palette_j <- c("#ca6354", "#8B58A4", "#fbd94c", "#709ed1", "#68c46e", "#ab9fb6",
               "#bfd49b", "#a5d7e8", "#b3b3b3", "#ded6b0", "#f4a261", "#b8b0a0")

plot_gene_usage(trb_data, "TRB JGene Usage by Sample", "TRB JGene", "figure1F_TRB_JGene_usage.pdf", palette_j)
plot_gene_usage(igh_data, "IGH JGene Usage by Sample", "IGH JGene", "figure1G_IGH_JGene_usage.pdf", palette_j)
