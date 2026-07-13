# Figure 1H: pooled IGH SHM rate by isotype.

source("analysis/figure1/R/00_setup.R")
suppressPackageStartupMessages({
  library(ggsignif)
  library(rlang)
})

read_igh <- function(path, sample_name) {
  read_csv_if_exists(path) |>
    distinct() |>
    mutate(Sample = sample_name)
}

data <- lapply(seq_len(nrow(sample_config)), function(i) {
  read_igh(sample_config$igh_expanded_csv[i], sample_config$display_sample[i])
}) |> bind_rows()

data <- data |>
  mutate(
    bestCGene = recode(
      bestCGene,
      "IGHA1" = "IGHA", "IGHA2" = "IGHA",
      "IGHE" = "IGHE", "IGHEP1" = "IGHE",
      "IGHG1" = "IGHG", "IGHG2" = "IGHG", "IGHG3" = "IGHG",
      "IGHG4" = "IGHG", "IGHGP" = "IGHG",
      .default = as.character(bestCGene)
    )
  ) |>
  filter(bestCGene != "", !is.na(bestCGene), bestCGene != "nan")

gene_order <- c("IGHA", "IGHE", "IGHG", "IGHD", "IGHM")

plot_violin <- function(data_sub, y_col, y_label, adjust_method = "BH", p_cut = 0.1) {
  present <- intersect(gene_order, unique(data_sub$bestCGene))
  data_sub <- data_sub |>
    mutate(bestCGene = factor(bestCGene, levels = present)) |>
    droplevels()
  data_sub[[y_col]] <- as.numeric(data_sub[[y_col]])
  data_sub <- data_sub |> filter(!!sym(y_col) <= 0.15)

  comparisons <- list()
  labels <- character(0)
  if ("IGHM" %in% present) {
    comparisons <- lapply(setdiff(present, "IGHM"), function(g) c("IGHM", g))
    raw_p <- sapply(comparisons, function(pair) {
      g1 <- data_sub |> filter(bestCGene == pair[1]) |> pull(!!sym(y_col))
      g2 <- data_sub |> filter(bestCGene == pair[2]) |> pull(!!sym(y_col))
      if (length(g1) > 1 && length(g2) > 1) suppressWarnings(wilcox.test(g1, g2)$p.value) else NA_real_
    })
    show_vals <- p.adjust(raw_p, method = adjust_method)
    labels <- ifelse(is.na(show_vals) | show_vals >= p_cut, "", paste0("p.adj=", format.pval(show_vals, digits = 3, eps = 2.2e-16)))
  }

  ggplot(data_sub, aes(x = bestCGene, y = !!sym(y_col), fill = bestCGene)) +
    geom_violin(trim = FALSE, alpha = 0.8, scale = "width", adjust = 1.5) +
    geom_boxplot(width = 0.1, lwd = 0.3, outlier.shape = NA, fill = "grey", color = "black") +
    stat_summary(fun = median, geom = "errorbar", width = 0.45, color = "black", size = 0.7) +
    geom_signif(comparisons = comparisons, annotations = labels, map_signif_level = FALSE,
                step_increase = 0.10, size = 0.5, textsize = 3, tip_length = 0.01) +
    scale_fill_manual(values = c("#E7298A", "#baabab", "#9aa8b2", "#d6e7f7", "#fedbc4")) +
    labs(x = "", y = y_label) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          legend.position = "none")
}

p <- plot_violin(data, "VDJ_mutation_rate", "SHM rate")
save_panel(p, "figure1H_IGH_SHM_rate_by_isotype.pdf", width = 4, height = 4)
