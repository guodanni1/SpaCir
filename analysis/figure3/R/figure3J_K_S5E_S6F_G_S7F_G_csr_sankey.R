# Figure 3J/K and supplementary Sankey panels.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))
source(file.path(figure3_dir, "R", "figure3H_I_spatial_clonefamily_branches.R"))

suppressPackageStartupMessages({
  library(networkD3)
  library(htmlwidgets)
})

prepare_sankey_edges <- function(edge_df, mode = c("region", "gc")) {
  mode <- match.arg(mode)
  edges <- prepare_branch_edges(edge_df)

  if (mode == "region") {
    edges %>%
      mutate(
        source_group = parent_GC_NewClusters,
        transition_group = GC_transition_type,
        target_group = node_GC_NewClusters
      )
  } else {
    edges %>%
      mutate(
        source_group = coalesce(str_extract(parent_GC_NewClusters, "GC\\d+"), parent_GC_NewClusters),
        transition_group = GC_transition_type,
        target_group = coalesce(str_extract(node_GC_NewClusters, "GC\\d+"), node_GC_NewClusters)
      )
  }
}

plot_transition_sankey <- function(sample_id, edge_df = NULL, mode = c("region", "gc"),
                                   output_tag = "figure3J") {
  mode <- match.arg(mode)
  if (is.null(edge_df)) {
    edge_df <- read_tree_edges(sample_id)
  }
  sankey_df <- prepare_sankey_edges(edge_df, mode = mode) %>%
    count(source_group, transition_group, target_group, name = "value") %>%
    filter(!is.na(source_group), !is.na(transition_group), !is.na(target_group), value > 0)

  nodes <- data.frame(name = unique(c(sankey_df$source_group, sankey_df$transition_group, sankey_df$target_group)))
  links <- bind_rows(
    sankey_df %>% transmute(source = match(source_group, nodes$name) - 1, target = match(transition_group, nodes$name) - 1, value),
    sankey_df %>% transmute(source = match(transition_group, nodes$name) - 1, target = match(target_group, nodes$name) - 1, value)
  )

  # networkD3 writes an HTML Sankey. This keeps the full interactive object for
  # GitHub users; export to PDF can be done manually if a static panel is needed.
  sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                          Value = "value", NodeID = "name", fontSize = 10, nodeWidth = 18)
  out_html <- file.path(output_dir, paste0(output_tag, "_", sample_id, "_", mode, "_transition_sankey.html"))
  saveWidget(sankey, out_html, selfcontained = TRUE)
  save_table(sankey_df, paste0(output_tag, "_", sample_id, "_", mode, "_transition_sankey.csv"))
  message("Saved: ", out_html)
  invisible(sankey_df)
}

summarize_clonefamily_composition <- function(clonefamily_df, clonefamily_col = "Clonefamily") {
  clonefamily_df %>%
    filter(!is.na(.data[[clonefamily_col]]), .data[[clonefamily_col]] != "None") %>%
    summarise(
      parent_c_gene_IGHG1 = mean(c_call == "IGHG1", na.rm = TRUE),
      node_c_gene_IGHG1 = mean(c_call == "IGHG1", na.rm = TRUE),
      late_gc_fraction = mean(GC_maturity == "LateGC", na.rm = TRUE),
      high_shm_fraction = mean(SHM_group == "High SHM", na.rm = TRUE)
    )
}

# Examples:
# plot_transition_sankey("Tonsil2", mode = "region", output_tag = "figure3J")
# plot_transition_sankey("Tonsil2", mode = "gc", output_tag = "figure3K")
