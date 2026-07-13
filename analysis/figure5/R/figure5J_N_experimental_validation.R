# Figure 5J-N: experimental validation panels.
#
# Input tables are public placeholders configured in
# config/experimental_validation_parameters.csv and figure5_inputs.csv.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure5", "R", "00_setup.R"))

plot_experimental_workflow <- function() {
  steps <- data.frame(
    step = factor(
      c("Lentiviral Prep & Infection", "Cell Culture & SHM Accumulation", "Sequencing & Data Analysis", "RT-PCR Validation"),
      levels = c("Lentiviral Prep & Infection", "Cell Culture & SHM Accumulation", "Sequencing & Data Analysis", "RT-PCR Validation")
    ),
    x = 1:4,
    label = c(
      "Target and control vectors\nKD and control viruses",
      "Ramos cells\n21-day continuous culture",
      "VDJ sequencing\nSHM and clonotype analysis",
      "RT-PCR\ncandidate validation"
    )
  )
  ggplot(steps, aes(x = .data$x, y = 1)) +
    geom_segment(aes(x = .data$x - 0.42, xend = .data$x + 0.42, yend = 1),
                 linewidth = 8, color = "#DFE8F7", lineend = "round") +
    geom_text(aes(label = .data$step), y = 1.18, size = 3.2, fontface = "bold") +
    geom_text(aes(label = .data$label), y = 0.86, size = 2.8) +
    geom_segment(data = steps[-nrow(steps), ], aes(x = .data$x + 0.48, xend = .data$x + 0.76, y = 1, yend = 1),
                 arrow = arrow(length = unit(2.5, "mm")), color = "#6E8AB6") +
    xlim(0.5, 4.5) +
    theme_void()
}

plot_validation_shm <- function(validation_shm_csv) {
  df <- read_csv_if_exists(validation_shm_csv, stringsAsFactors = FALSE)
  ggplot(df, aes(x = .data$target_gene, y = .data$shm_frequency, fill = .data$condition)) +
    geom_violin(position = position_dodge(width = 0.75), trim = FALSE, alpha = 0.75) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.14, outlier.size = 0.5) +
    labs(x = NULL, y = "SHM frequency", fill = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_qpcr_validation <- function(validation_qpcr_csv, assay_gene) {
  df <- read_csv_if_exists(validation_qpcr_csv, stringsAsFactors = FALSE) %>%
    filter(.data$assay_gene == assay_gene)
  ggplot(df, aes(x = .data$target_gene, y = .data$relative_expression, color = .data$condition)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 1.8) +
    stat_summary(aes(group = interaction(.data$target_gene, .data$condition)),
                 fun = median, geom = "crossbar", width = 0.28,
                 position = position_dodge(width = 0.55), color = "black") +
    labs(title = assay_gene, x = NULL, y = "Relative Expression", color = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

run_experimental_validation_panels <- function(sample_id = "Tonsil2") {
  cfg <- sample_row(sample_id)

  p_workflow <- plot_experimental_workflow()
  save_panel(p_workflow, "figure5J_experimental_workflow.pdf", width = 8.6, height = 2.2)

  p_shm <- plot_validation_shm(cfg$validation_shm_csv)
  save_panel(p_shm, "figure5K_shm_frequency_validation.pdf", width = 4.2, height = 3.2)

  for (gene in c("AICDA", "BCL6", "MSH6")) {
    p <- plot_qpcr_validation(cfg$validation_qpcr_csv, assay_gene = gene)
    panel <- switch(gene, AICDA = "L", BCL6 = "M", MSH6 = "N")
    save_panel(p, paste0("figure5", panel, "_", gene, "_qpcr_validation.pdf"), width = 2.7, height = 2.5)
  }

  invisible(TRUE)
}

# Example:
# run_experimental_validation_panels("Tonsil2")
