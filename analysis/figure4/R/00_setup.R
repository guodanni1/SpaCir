# Shared settings and helper functions for SpaCir Figure 4 and Figure S8-S10.
# Run scripts from the repository root, or set SPACIR_ROOT before sourcing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

repo_root <- Sys.getenv("SPACIR_ROOT", unset = getwd())
figure4_dir <- file.path(repo_root, "analysis", "figure4")
input_config <- file.path(figure4_dir, "config", "figure4_inputs.csv")
output_dir <- file.path(figure4_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_config <- read.csv(input_config, stringsAsFactors = FALSE)

maturity_levels <- c("EarlyGC", "LateGC")
maturity_colors <- c("EarlyGC" = "#66A61E", "LateGC" = "#903586")

# Figure 4B clone-class definition:
#   Private: the IGH clone is detected in exactly one GC site.
#   Shared: the IGH clone is detected in more than one GC site.
#   Expanded: total IGH UMI is above the sample-specific 95th percentile.
#   Unexpanded: total IGH UMI is at or below the same 95th percentile.
clone_class_levels <- c(
  "A (Private Expanded)",
  "B (Shared Expanded)",
  "C (Private Unexpanded)",
  "D (Shared Unexpanded)"
)
clone_class_colors <- c(
  "A (Private Expanded)" = "#E31A1C",
  "B (Shared Expanded)" = "#B49AC8",
  "C (Private Unexpanded)" = "#FEC44F",
  "D (Shared Unexpanded)" = "#2B8CBE"
)

gc_colors <- c(
  "GC1" = "#FF8000", "GC2" = "#E7298A", "GC3" = "#4662D7",
  "GC4" = "#BD10E0", "GC5" = "#FFFF00", "GC6" = "#FF6347",
  "GC7" = "#00FFFF", "GC8" = "#0000FF", "GC9" = "#00008B",
  "GC10" = "#FFE4B5", "GC11" = "#FF2200", "GC12" = "#ADFF2F",
  "GC13" = "#4CF5C8", "GC14" = "#FA552B"
)

# GC maturity marker sets used for Figure 4D and all EarlyGC/LateGC grouping.
# Each GC part is scored with AddModuleScore for the two marker lists below.
# The GC part is labeled EarlyGC when mean(Early score) >= mean(Late score),
# otherwise it is labeled LateGC.
early_gc_markers <- c("CD83", "CD86", "CXCR5", "GPR183", "SLAMF1", "BCL6")
late_gc_markers <- c(
  "IRF4", "PRDM1", "XBP1", "ZBTB20", "FOXP3", "DUSP2",
  "IRF8", "GADD45B", "JCHAIN", "TOX2", "SDC1"
)
gc_maturity_score_nbin <- 10
clone_expanded_umi_quantile <- 0.95

b_cell_types <- c(
  "B_naive", "B_mem", "B_activated", "B_preGC", "B_Cycling",
  "B_GC_LZ", "B_GC_DZ", "B_GC_prePB", "B_IFN", "B_plasma"
)
t_cell_types <- c(
  "T_CD4+", "T_CD4+_naive", "T_CD4+_TfH", "T_CD4+_TfH_GC",
  "T_CD8+_CD161+", "T_CD8+_cytotoxic", "T_CD8+_naive",
  "T_TfR", "T_Treg", "T_TIM3+", "NKT"
)

cell_type_palette <- c(
  "#CA6354", "#338FC4", "#FBD94C", "#709ED1", "#68C46E",
  "#AB9FB6", "#91BA90", "#F8A34B", "#D5DCA0", "#C38876",
  "#91D8F6", "#DEA4CA", "#F58888", "#64CCD3", "#FEE8E6",
  "#F2EDBA", "#D6CDE5", "#4C9E72", "#B9DDA4", "#F29987",
  "#90C3E8", "#C3BC3F", "#7ACA97", "#FEDBC4", "#8FBAC8",
  "#BAABAB", "#E7298A", "#EFB3AF", "#9AA8B2", "#C0B39E",
  "#D6E7F7", "#2B3D26", "#8B58A4", "#E6550D"
)

sample_row <- function(sample_id) {
  row <- sample_config[sample_config$sample_id == sample_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one row in figure4_inputs.csv for sample: ", sample_id, call. = FALSE)
  }
  row
}

read_csv_if_exists <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  read.csv(path, ...)
}

save_panel <- function(plot, filename, width, height) {
  out <- file.path(output_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved: ", out)
}

save_table <- function(df, filename) {
  out <- file.path(output_dir, filename)
  write.csv(df, out, row.names = FALSE)
  message("Saved: ", out)
}

load_sample_objects <- function(sample_id) {
  cfg <- sample_row(sample_id)
  list(
    cfg = cfg,
    object = if (file.exists(cfg$seurat_object_rds)) readRDS(cfg$seurat_object_rds) else NULL,
    filtered_object = if (file.exists(cfg$filtered_object_rds)) readRDS(cfg$filtered_object_rds) else NULL
  )
}

ensure_gc_maturity_factor <- function(df) {
  df %>%
    filter(.data$GC_Maturity_Label %in% maturity_levels) %>%
    mutate(GC_Maturity_Label = factor(.data$GC_Maturity_Label, levels = maturity_levels))
}

assign_gc_maturity <- function(seurat_obj,
                               gc_col = "GC_part",
                               assay = "SCT",
                               early_markers = early_gc_markers,
                               late_markers = late_gc_markers,
                               score_nbin = gc_maturity_score_nbin) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for GC maturity scoring.", call. = FALSE)
  }
  seurat_obj <- subset(seurat_obj, subset = .data[[gc_col]] != "Others")
  Seurat::DefaultAssay(seurat_obj) <- assay

  # Score each spot with the explicit EarlyGC and LateGC gene sets.
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = list(early_markers), name = "Score_EarlyGC", nbin = score_nbin)
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = list(late_markers), name = "Score_LateGC", nbin = score_nbin)

  score_df <- Seurat::FetchData(seurat_obj, vars = c("Score_EarlyGC1", "Score_LateGC1", gc_col))
  colnames(score_df)[colnames(score_df) == gc_col] <- "GC_part"
  summary_df <- score_df %>%
    group_by(.data$GC_part) %>%
    summarise(
      mean_Early = mean(.data$Score_EarlyGC1, na.rm = TRUE),
      mean_Late = mean(.data$Score_LateGC1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Match the original analysis: the larger mean module score determines
      # the maturity state assigned to the whole annotated GC part.
      GC_Maturity_Label = ifelse(.data$mean_Early >= .data$mean_Late, "EarlyGC", "LateGC"),
      GC_Maturity_Label = factor(.data$GC_Maturity_Label, levels = maturity_levels)
    )

  meta <- seurat_obj@meta.data
  meta$GC_part <- as.character(meta[[gc_col]])
  meta <- meta %>%
    left_join(summary_df %>% select(.data$GC_part, .data$GC_Maturity_Label), by = "GC_part")
  seurat_obj$GC_Maturity_Label <- meta$GC_Maturity_Label
  list(
    object = seurat_obj,
    summary = summary_df,
    parameters = data.frame(
      gc_column = gc_col,
      assay = assay,
      early_gc_markers = paste(early_markers, collapse = ";"),
      late_gc_markers = paste(late_markers, collapse = ";"),
      module_score_nbin = score_nbin,
      assignment_rule = "EarlyGC if mean_Early >= mean_Late within GC_part; otherwise LateGC",
      stringsAsFactors = FALSE
    )
  )
}

plot_gc_maturity_dotplot <- function(seurat_obj, group_col = "GC_part") {
  Seurat::DotPlot(
    seurat_obj,
    features = c("Score_EarlyGC1", "Score_LateGC1"),
    scale = FALSE,
    group.by = group_col
  ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

plot_spatial_gc_parts <- function(seurat_obj, group_col = "GC_part", label = TRUE) {
  Seurat::SpatialDimPlot(
    seurat_obj,
    group.by = group_col,
    label = label,
    label.size = 4,
    pt.size.factor = 1.6,
    cols = gc_colors
  ) +
    ggtitle(NULL)
}

celltype_columns <- function(meta_data, prefix = "q05cell_abundance_w_sf_") {
  grep(paste0("^", prefix), colnames(meta_data), value = TRUE)
}

clean_celltype_name <- function(x, prefix = "q05cell_abundance_w_sf_") {
  gsub(paste0("^", prefix), "", x)
}

build_celltype_colors <- function(cell_types) {
  setNames(rep(cell_type_palette, length.out = length(cell_types)), cell_types)
}

star_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

plot_violin_with_p <- function(df, group_col, value_col, fill_colors, y_label = NULL, title = NULL) {
  df <- df %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]])) %>%
    mutate("{group_col}" := factor(.data[[group_col]], levels = unique(.data[[group_col]])))

  groups <- levels(df[[group_col]])
  comparisons <- if (length(groups) >= 2) combn(groups, 2, simplify = FALSE) else list()
  labels <- map_chr(comparisons, function(pair) {
    p <- tryCatch(
      wilcox.test(df[[value_col]][df[[group_col]] == pair[1]], df[[value_col]][df[[group_col]] == pair[2]])$p.value,
      error = function(e) NA_real_
    )
    sprintf("p=%.3g", p)
  })

  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "grey90", color = "black") +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = y_label %||% value_col, title = title) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

  if (length(comparisons) && requireNamespace("ggsignif", quietly = TRUE)) {
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons,
      annotations = labels,
      map_signif_level = FALSE,
      step_increase = 0.1,
      tip_length = 0.01,
      textsize = 3
    )
  }
  p
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

extract_mutation_rate <- function(alignment_value) {
  if (!is.character(alignment_value) || is.na(alignment_value) || alignment_value == "") {
    return(NA_real_)
  }
  first_segment <- strsplit(gsub("^\\{|\\}$", "", alignment_value), "\\}\\{")[[1]][1]
  parts <- strsplit(first_segment, "\\|")[[1]]
  if (length(parts) < 6) {
    return(NA_real_)
  }
  start <- suppressWarnings(as.numeric(parts[1]))
  end <- suppressWarnings(as.numeric(parts[2]))
  if (!is.finite(start) || !is.finite(end) || end < start) {
    return(NA_real_)
  }
  mutation_count <- length(unlist(regmatches(parts[6], gregexpr("[0-9]+", parts[6]))))
  mutation_count / (end - start + 1)
}

prepare_shm_table <- function(alignment_csv, metadata_df, group_col) {
  read_csv_if_exists(alignment_csv, stringsAsFactors = FALSE) %>%
    distinct() %>%
    left_join(metadata_df, by = "spatial_barcode_ori") %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Others") %>%
    mutate(
      SHM = if ("VDJ_mutation_rate" %in% colnames(.)) {
        as.numeric(.data$VDJ_mutation_rate)
      } else if ("allVAlignments" %in% colnames(.)) {
        vapply(.data$allVAlignments, extract_mutation_rate, numeric(1))
      } else {
        NA_real_
      }
    )
}

matrix_to_long <- function(mat, feature_col = "Clonotype", value_col = "UMI") {
  out <- as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE)
  colnames(out) <- c(feature_col, "spot", value_col)
  out %>% filter(.data[[value_col]] > 0)
}

get_assay_matrix <- function(seurat_obj, assay_name) {
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop("Assay not found: ", assay_name, call. = FALSE)
  }
  as.matrix(seurat_obj@assays[[assay_name]]@data)
}

classify_igh_clones <- function(seurat_obj,
                                assay_name,
                                gc_col = "GC_part",
                                umi_quantile = clone_expanded_umi_quantile) {
  mat <- get_assay_matrix(seurat_obj, assay_name)
  gc_map <- data.frame(
    spot = rownames(seurat_obj@meta.data),
    GC_part = seurat_obj@meta.data[[gc_col]],
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(.data$GC_part), .data$GC_part != "Others")

  long <- matrix_to_long(mat, feature_col = "IGH_clonotype", value_col = "UMI_count") %>%
    left_join(gc_map, by = "spot") %>%
    filter(!is.na(.data$GC_part))

  summary <- long %>%
    group_by(.data$IGH_clonotype) %>%
    summarise(
      num_GC_sites = n_distinct(.data$GC_part),
      total_UMI = sum(.data$UMI_count),
      .groups = "drop"
    )
  # The expansion cutoff is sample-specific: calculate the requested UMI
  # percentile across all IGH clonotypes detected in the current sample.
  cutoff <- quantile(summary$total_UMI, umi_quantile, na.rm = TRUE)
  summary <- summary %>%
    mutate(
      class = case_when(
        .data$num_GC_sites == 1 & .data$total_UMI > cutoff ~ "A (Private Expanded)",
        .data$num_GC_sites > 1 & .data$total_UMI > cutoff ~ "B (Shared Expanded)",
        .data$num_GC_sites == 1 & .data$total_UMI <= cutoff ~ "C (Private Unexpanded)",
        .data$num_GC_sites > 1 & .data$total_UMI <= cutoff ~ "D (Shared Unexpanded)",
        TRUE ~ NA_character_
      ),
      class = factor(.data$class, levels = clone_class_levels),
      total_UMI_log = log10(.data$total_UMI + 1),
      expanded_cutoff_quantile = umi_quantile,
      expanded_cutoff_total_UMI = as.numeric(cutoff)
    )
  parameters <- data.frame(
    assay = assay_name,
    gc_column = gc_col,
    private_rule = "num_GC_sites == 1",
    shared_rule = "num_GC_sites > 1",
    expanded_rule = paste0("total_UMI > sample ", umi_quantile * 100, "th percentile"),
    unexpanded_rule = paste0("total_UMI <= sample ", umi_quantile * 100, "th percentile"),
    expanded_cutoff_quantile = umi_quantile,
    expanded_cutoff_total_UMI = as.numeric(cutoff),
    class_A = "Private Expanded",
    class_B = "Shared Expanded",
    class_C = "Private Unexpanded",
    class_D = "Shared Unexpanded",
    stringsAsFactors = FALSE
  )
  list(summary = summary, long = long, cutoff = cutoff, parameters = parameters)
}

plot_clone_detection <- function(clone_summary, title = "IGH Clone Detection in GC Sites") {
  max_site <- max(clone_summary$num_GC_sites, na.rm = TRUE)
  fill_values <- colorRampPalette(c("#7B3294", "#C2A5CF", "#F7F7F7", "#A6DBA0", "#008837"))(max_site)
  ggplot(clone_summary, aes(x = factor(.data$num_GC_sites), y = .data$total_UMI_log, fill = .data$num_GC_sites)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.85, color = "black") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_gradientn(colors = fill_values, name = "GC sites") +
    labs(title = title, x = "Number of GC sites IGH clone detected", y = "log10(Total UMI + 1)") +
    theme_classic()
}

plot_clone_class_bar_by_gc <- function(igh_long, clone_summary) {
  igh_long %>%
    left_join(clone_summary %>% select(.data$IGH_clonotype, .data$class), by = "IGH_clonotype") %>%
    filter(!is.na(.data$class)) %>%
    group_by(.data$GC_part, .data$class) %>%
    summarise(total_UMI = sum(.data$UMI_count), .groups = "drop") %>%
    group_by(.data$GC_part) %>%
    mutate(percentage = 100 * .data$total_UMI / sum(.data$total_UMI)) %>%
    ungroup() %>%
    ggplot(aes(x = .data$GC_part, y = .data$percentage, fill = .data$class)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = clone_class_colors, drop = FALSE) +
    labs(x = NULL, y = "Percentage (%)", fill = "Class") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_cgene_bar_by_gc <- function(alignment_csv, metadata_df) {
  read_csv_if_exists(alignment_csv, stringsAsFactors = FALSE) %>%
    distinct() %>%
    left_join(metadata_df, by = "spatial_barcode_ori") %>%
    filter(!is.na(.data$GC_part), .data$GC_part != "Others", !is.na(.data$bestCGene), .data$bestCGene != "nan") %>%
    group_by(.data$GC_part, .data$bestCGene) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data$GC_part) %>%
    mutate(percentage = 100 * .data$n / sum(.data$n)) %>%
    ungroup() %>%
    ggplot(aes(x = .data$GC_part, y = .data$percentage, fill = .data$bestCGene)) +
    geom_col(width = 0.85) +
    labs(x = NULL, y = "Percentage (%)", fill = "C gene") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

calc_abundance_totals <- function(meta_df, prefix = "q05cell_abundance_w_sf_") {
  t_cols <- intersect(paste0(prefix, t_cell_types), colnames(meta_df))
  b_cols <- intersect(paste0(prefix, b_cell_types), colnames(meta_df))
  meta_df %>%
    mutate(spot = rownames(meta_df)) %>%
    rowwise() %>%
    mutate(
      T_total = if (length(t_cols)) sum(c_across(all_of(t_cols)), na.rm = TRUE) else NA_real_,
      B_total = if (length(b_cols)) sum(c_across(all_of(b_cols)), na.rm = TRUE) else NA_real_
    ) %>%
    ungroup()
}

calc_diversity_norm <- function(seurat_obj, assay_name, meta_df, norm_col) {
  mat <- get_assay_matrix(seurat_obj, assay_name)
  common_spots <- intersect(colnames(mat), meta_df$spot)
  mat <- mat[, common_spots, drop = FALSE]
  meta_use <- meta_df[match(common_spots, meta_df$spot), , drop = FALSE]
  long <- matrix_to_long(mat, feature_col = "Clonotype", value_col = "UMI") %>%
    left_join(meta_use %>% select(.data$spot, .data$GC_Maturity_Label, all_of(norm_col)), by = "spot") %>%
    ensure_gc_maturity_factor()

  long %>%
    group_by(.data$GC_Maturity_Label, .data$spot, .data$Clonotype) %>%
    summarise(UMI = sum(.data$UMI), .groups = "drop") %>%
    group_by(.data$GC_Maturity_Label, .data$spot) %>%
    summarise(
      Shannon = {
        p <- .data$UMI / sum(.data$UMI)
        -sum(p * log(p), na.rm = TRUE)
      },
      UniqueClonotypes = n_distinct(.data$Clonotype),
      .groups = "drop"
    ) %>%
    left_join(meta_use %>% select(.data$spot, all_of(norm_col)), by = "spot") %>%
    rename(NormFactor = all_of(norm_col)) %>%
    mutate(
      Shannon_norm = ifelse(.data$NormFactor > 0, .data$Shannon / .data$NormFactor, NA_real_),
      UniqueClonotypes_norm = ifelse(.data$NormFactor > 0, .data$UniqueClonotypes / .data$NormFactor, NA_real_)
    )
}

plot_diversity_norm_pair <- function(igh_metrics, trb_metrics) {
  bind_rows(
    igh_metrics %>% transmute(Chain = "IGH", GC_Maturity_Label, Shannon_norm),
    trb_metrics %>% transmute(Chain = "TRB", GC_Maturity_Label, Shannon_norm)
  ) %>%
    ensure_gc_maturity_factor() %>%
    ggplot(aes(x = .data$GC_Maturity_Label, y = .data$Shannon_norm, fill = .data$GC_Maturity_Label)) +
    geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "grey90", color = "black") +
    facet_wrap(~ Chain, scales = "free_y") +
    scale_fill_manual(values = maturity_colors) +
    labs(x = NULL, y = "Normalized Shannon Diversity") +
    theme_classic() +
    theme(legend.position = "none")
}

plot_cumulative_clone_frequency <- function(seurat_obj, assay_name) {
  mat <- get_assay_matrix(seurat_obj, assay_name)
  meta <- seurat_obj@meta.data %>%
    mutate(spot = rownames(seurat_obj@meta.data)) %>%
    select(.data$spot, .data$GC_Maturity_Label) %>%
    ensure_gc_maturity_factor()
  matrix_to_long(mat, feature_col = "Clonotype", value_col = "UMI") %>%
    left_join(meta, by = "spot") %>%
    filter(!is.na(.data$GC_Maturity_Label)) %>%
    group_by(.data$GC_Maturity_Label, .data$Clonotype) %>%
    summarise(UMI = sum(.data$UMI), .groups = "drop") %>%
    group_by(.data$GC_Maturity_Label) %>%
    mutate(Frequency = .data$UMI / sum(.data$UMI)) %>%
    arrange(.data$GC_Maturity_Label, desc(.data$Frequency)) %>%
    mutate(clone_rank = row_number(), cumulative_frequency = cumsum(.data$Frequency)) %>%
    ungroup() %>%
    ggplot(aes(x = .data$clone_rank, y = .data$cumulative_frequency, color = .data$GC_Maturity_Label)) +
    geom_line(linewidth = 1.1) +
    scale_x_log10() +
    scale_color_manual(values = maturity_colors) +
    labs(x = "Number of IGH clones", y = "Cumulative Frequency", color = NULL) +
    theme_classic()
}

plot_late_clone_composition <- function(seurat_obj, assay_name, top_n = 100) {
  mat <- get_assay_matrix(seurat_obj, assay_name)
  meta <- seurat_obj@meta.data %>%
    mutate(spot = rownames(seurat_obj@meta.data)) %>%
    select(.data$spot, .data$GC_Maturity_Label) %>%
    ensure_gc_maturity_factor()
  clone_freq <- matrix_to_long(mat, feature_col = "Clonotype", value_col = "UMI") %>%
    left_join(meta, by = "spot") %>%
    group_by(.data$GC_Maturity_Label, .data$Clonotype) %>%
    summarise(UMI = sum(.data$UMI), .groups = "drop") %>%
    group_by(.data$GC_Maturity_Label) %>%
    mutate(Percent = .data$UMI / sum(.data$UMI)) %>%
    ungroup()
  top_clones <- clone_freq %>%
    filter(.data$GC_Maturity_Label == "LateGC") %>%
    arrange(desc(.data$Percent)) %>%
    slice_head(n = top_n) %>%
    pull(.data$Clonotype)
  clone_freq %>%
    mutate(clone_group = ifelse(.data$Clonotype %in% top_clones, .data$Clonotype, "Other")) %>%
    group_by(.data$GC_Maturity_Label, .data$clone_group) %>%
    summarise(Percent = sum(.data$Percent), .groups = "drop") %>%
    ggplot(aes(x = .data$GC_Maturity_Label, y = .data$Percent, fill = .data$clone_group)) +
    geom_col(width = 0.85) +
    labs(x = NULL, y = "Relative UMI Percent") +
    theme_classic() +
    theme(legend.position = "none")
}

plot_two_state_line_with_stars <- function(point_df, raw_df, group_col, value_col, line_group_col,
                                           y_label, color_values = NULL, title = NULL) {
  group_col <- rlang::ensym(group_col)
  value_col <- rlang::ensym(value_col)
  line_group_col <- rlang::ensym(line_group_col)

  point_df <- point_df %>%
    filter(!!group_col %in% maturity_levels) %>%
    mutate(!!group_col := factor(!!group_col, levels = maturity_levels))

  pvals <- raw_df %>%
    filter(!!group_col %in% maturity_levels) %>%
    group_by(!!line_group_col) %>%
    summarise(
      p.value = tryCatch(
        wilcox.test(!!value_col ~ !!group_col)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(star = star_label(.data$p.value))

  label_df <- point_df %>%
    filter(!!group_col == "LateGC") %>%
    left_join(pvals, by = rlang::as_name(line_group_col)) %>%
    mutate(label_y = !!value_col + max(point_df[[rlang::as_name(value_col)]], na.rm = TRUE) * 0.06)

  p <- ggplot(point_df, aes(x = !!group_col, y = !!value_col, group = !!line_group_col, color = !!line_group_col)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    geom_point(size = 2) +
    geom_text(
      data = label_df %>% filter(.data$star != ""),
      aes(y = .data$label_y, label = .data$star),
      show.legend = FALSE,
      size = 3
    ) +
    labs(x = "GC Maturity Stage", y = y_label, color = NULL, title = title) +
    theme_classic()
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values)
  }
  p
}

prepare_cgene_trend <- function(alignment_csv, metadata_df) {
  raw <- read_csv_if_exists(alignment_csv, stringsAsFactors = FALSE) %>%
    distinct() %>%
    left_join(metadata_df, by = "spatial_barcode_ori") %>%
    filter(.data$GC_Maturity_Label %in% maturity_levels, !is.na(.data$bestCGene), .data$bestCGene != "nan") %>%
    group_by(.data$spatial_barcode_ori, .data$bestCGene, .data$GC_Maturity_Label) %>%
    summarise(UMI_count = n_distinct(.data$UMI), .groups = "drop")
  trend <- raw %>%
    group_by(.data$bestCGene, .data$GC_Maturity_Label) %>%
    summarise(total_UMI = sum(.data$UMI_count), .groups = "drop") %>%
    group_by(.data$GC_Maturity_Label) %>%
    mutate(UMI_percent = 100 * .data$total_UMI / sum(.data$total_UMI)) %>%
    ungroup()
  list(raw = raw, trend = trend)
}

prepare_celltype_trend <- function(seurat_obj) {
  meta <- seurat_obj@meta.data
  meta$spot <- rownames(meta)
  cols <- celltype_columns(meta)
  long <- meta %>%
    select(.data$spot, .data$GC_Maturity_Label, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "CellType", values_to = "Abundance") %>%
    mutate(
      CellType = clean_celltype_name(.data$CellType),
      GC_Maturity_Label = factor(.data$GC_Maturity_Label, levels = maturity_levels),
      CellCategory = case_when(
        .data$CellType %in% b_cell_types ~ "B cell",
        .data$CellType %in% t_cell_types ~ "T cell",
        TRUE ~ "Other cell"
      )
    ) %>%
    filter(.data$GC_Maturity_Label %in% maturity_levels)
  totals <- long %>% group_by(.data$spot) %>% summarise(total = sum(.data$Abundance, na.rm = TRUE), .groups = "drop")
  raw <- long %>%
    left_join(totals, by = "spot") %>%
    mutate(Fraction = ifelse(.data$total > 0, 100 * .data$Abundance / .data$total, NA_real_))
  trend <- raw %>%
    group_by(.data$GC_Maturity_Label, .data$CellType, .data$CellCategory) %>%
    summarise(MeanFraction = mean(.data$Fraction, na.rm = TRUE), .groups = "drop")
  list(raw = raw, trend = trend, colors = build_celltype_colors(sort(unique(raw$CellType))))
}

plot_celltype_category_trend <- function(celltype_data, category) {
  raw <- celltype_data$raw %>% filter(.data$CellCategory == category)
  trend <- celltype_data$trend %>% filter(.data$CellCategory == category)
  plot_two_state_line_with_stars(
    point_df = trend,
    raw_df = raw,
    group_col = GC_Maturity_Label,
    value_col = MeanFraction,
    line_group_col = CellType,
    y_label = "Relative Abundance per spot (%)",
    color_values = celltype_data$colors,
    title = paste("Cell Type Trend -", sub(" cell$", " cells", category))
  )
}

plot_go_dotplot <- function(go_table, description_order = NULL, group_col = "group") {
  df <- go_table %>%
    mutate(
      p.adjust = as.numeric(.data$p.adjust),
      GeneRatioValue = ifelse(grepl("/", .data$GeneRatio),
                              map_dbl(strsplit(.data$GeneRatio, "/"), ~ as.numeric(.x[1]) / as.numeric(.x[2])),
                              as.numeric(.data$GeneRatio)),
      "{group_col}" := factor(.data[[group_col]], levels = maturity_levels)
    )
  if (!is.null(description_order)) {
    df$Description <- factor(df$Description, levels = rev(description_order))
  }
  ggplot(df, aes(x = .data[[group_col]], y = .data$Description)) +
    geom_point(aes(size = .data$GeneRatioValue, color = .data$p.adjust)) +
    scale_color_gradient(low = "#CC3333", high = "#6699CC", limits = c(0, 0.05), oob = scales::squish) +
    labs(x = "GC Stage", y = "GO Term", size = "Gene Ratio", color = "p.adjust") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_wgcna_hub_go_dotplot <- function(go_table) {
  module_color_map <- c(
    "black" = "grey30", "blue" = "#3399FF", "green" = "#88CC00",
    "pink" = "pink", "red" = "#D7263D", "turquoise" = "#03C1C1",
    "yellow" = "yellow", "brown" = "brown", "magenta" = "#FF00FF",
    "salmon" = "#FA8072", "royalblue" = "#4169E1", "lightcyan" = "#E0FFFF"
  )
  df <- go_table %>%
    mutate(
      combined_score = as.numeric(.data$combined_score),
      module = factor(.data$module, levels = unique(.data$module)),
      Description = factor(.data$Description, levels = rev(unique(.data$Description)))
    )
  ggplot(df, aes(x = .data$module, y = .data$Description)) +
    geom_point(aes(size = .data$combined_score, color = .data$module)) +
    scale_color_manual(values = module_color_map, drop = FALSE) +
    scale_size_continuous(name = "Combined Score") +
    labs(x = "Module", y = "GO Term") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
