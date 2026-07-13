# Figure 1I/L: TRB and IGH clone size composition.

source("analysis/figure1/R/00_setup.R")

classify_clone_size <- function(umi_per_clone) {
  cut(
    umi_per_clone,
    breaks = c(1, 2, 6, 20, Inf),
    labels = c("Not Expanded", "Small", "Medium", "Hyperexpanded"),
    right = FALSE
  )
}

clone_size_composition <- function(obj, sample_name, chain) {
  assay_name <- paste0("aaSeqCDR3_Lv1_UMI3_", chain, "_20G")
  mat <- as.matrix(obj@assays[[assay_name]]@data)
  umi_per_clone <- rowSums(mat)
  tab <- as.data.frame(table(classify_clone_size(umi_per_clone)))
  colnames(tab) <- c("Classification", "Count")
  tab |>
    mutate(
      Sample = sample_name,
      Chain = chain,
      Percentage = Count / sum(Count) * 100,
      Classification = factor(Classification, levels = c("Not Expanded", "Small", "Medium", "Hyperexpanded"))
    )
}

all_comp <- lapply(seq_len(nrow(sample_config)), function(i) {
  obj <- get(sample_config$object_name[i])
  bind_rows(
    clone_size_composition(obj, sample_config$display_sample[i], "IGH"),
    clone_size_composition(obj, sample_config$display_sample[i], "TRB")
  )
}) |> bind_rows()

plot_comp <- function(chain, filename) {
  df <- all_comp |> filter(Chain == chain)
  df$Sample <- factor(df$Sample, levels = sample_levels)
  p <- ggplot(df, aes(x = Sample, y = Percentage, fill = Classification)) +
    geom_bar(stat = "identity", width = 0.9) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_stack(vjust = 0.5), color = "white", size = 3) +
    scale_fill_manual(values = c("Not Expanded" = "#99C945", "Small" = "orange",
                                 "Medium" = "blue", "Hyperexpanded" = "red")) +
    labs(title = paste(chain, "Clone Size"), x = "", y = "Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  save_panel(p, filename, width = 4, height = 4)
}

plot_comp("IGH", "figure1I_IGH_clone_size_composition.pdf")
plot_comp("TRB", "figure1L_TRB_clone_size_composition.pdf")
