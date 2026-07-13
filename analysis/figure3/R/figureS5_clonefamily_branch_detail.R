# Figure S5B-D: selected clonefamily branch details.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))
source(file.path(figure3_dir, "R", "figure3H_I_spatial_clonefamily_branches.R"))

suppressPackageStartupMessages({
  library(ggseqlogo)
})

plot_branch_sequence_logo <- function(sequence_df, branch_col = "branch_group",
                                      sequence_col = "junction_aa", output_tag = "figureS5B") {
  # One sequence-logo panel per branch. Input should already be filtered to the
  # clonefamily branch shown in Figure S5.
  p <- ggplot(sequence_df, aes(x = .data[[sequence_col]], y = .data[[branch_col]])) +
    theme_void()
  logos <- split(sequence_df[[sequence_col]], sequence_df[[branch_col]])
  out_pdf <- file.path(output_dir, paste0(output_tag, "_sequence_logo.pdf"))
  pdf(out_pdf, width = 6, height = max(2, length(logos) * 1.2))
  for (nm in names(logos)) {
    print(ggseqlogo(logos[[nm]]) + ggtitle(nm))
  }
  dev.off()
  message("Saved: ", out_pdf)
  invisible(p)
}

plot_selected_branch_arrows <- function(sample_id, selected_branches = NULL, output_tag = "figureS5D") {
  edges <- read_tree_edges(sample_id)
  if (!is.null(selected_branches)) {
    edges <- edges %>% filter(branch_group %in% selected_branches)
  }
  plot_branch_arrows(sample_id, edges, output_tag = output_tag, parent_pattern = NULL, node_pattern = NULL)
}

# Examples:
# plot_selected_branch_arrows("LN2", selected_branches = c("Branch_9", "Branch_11", "Branch_12"))
