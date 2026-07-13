# Figure S1D/E: IGH and TRB V gene usage by sample.

source("analysis/figure1/R/00_setup.R")

prepare_gene_usage <- function(data, sample_name, gene_col, top_n = 15) {
  gene_count <- table(data[[gene_col]])
  gene_data <- as.data.frame(gene_count)
  colnames(gene_data) <- c("Gene", "Count")
  gene_data <- gene_data |>
    filter(!is.na(Gene), Gene != "", Gene != "nan") |>
    mutate(Percentage = Count / sum(Count) * 100) |>
    arrange(desc(Percentage))
  keep <- gene_data |> slice_head(n = top_n) |> pull(Gene)
  gene_data |>
    mutate(Gene = ifelse(Gene %in% keep, as.character(Gene), "Others"),
           Sample = sample_name) |>
    group_by(Sample, Gene) |>
    summarise(Count = sum(Count), Percentage = sum(Percentage), .groups = "drop")
}

plot_v_usage <- function(data, title, fill_title, filename) {
  data$Sample <- factor(data$Sample, levels = sample_levels)
  p <- ggplot(data, aes(x = Sample, y = Percentage, fill = Gene)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = title, x = "Sample", y = "Percentage (%)", fill = fill_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          panel.grid = element_blank())
  save_panel(p, filename, width = 5.5, height = 4)
}

igh_v <- lapply(seq_len(nrow(sample_config)), function(i) {
  read_csv_if_exists(sample_config$igh_expanded_csv[i]) |>
    distinct() |>
    prepare_gene_usage(sample_config$display_sample[i], "bestVGene")
}) |> bind_rows()

trb_v <- lapply(seq_len(nrow(sample_config)), function(i) {
  read_csv_if_exists(sample_config$trb_expanded_csv[i]) |>
    distinct() |>
    prepare_gene_usage(sample_config$display_sample[i], "bestVGene")
}) |> bind_rows()

plot_v_usage(igh_v, "IGH VGene Usage by Sample", "IGH VGene", "figureS1D_IGH_VGene_usage_top15.pdf")
plot_v_usage(trb_v, "TRB VGene Usage by Sample", "TRB VGene", "figureS1E_TRB_VGene_usage_top15.pdf")
