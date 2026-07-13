# Shared settings and helper functions for SpaCir Figure 5 and Figure S11-S12.
# Run scripts from the repository root, or set SPACIR_ROOT before sourcing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_root <- Sys.getenv("SPACIR_ROOT", unset = getwd())
figure5_dir <- file.path(repo_root, "analysis", "figure5")
input_config <- file.path(figure5_dir, "config", "figure5_inputs.csv")
gene_config <- file.path(figure5_dir, "config", "shm_gene_sets.csv")
output_dir <- file.path(figure5_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_config <- read.csv(input_config, stringsAsFactors = FALSE)
shm_gene_table <- read.csv(gene_config, stringsAsFactors = FALSE)

shm_levels <- c("No_Mutation", "Low_Mutation", "High_Mutation")
shm_colors <- c(
  "No_Mutation" = "#4E79A7",
  "Low_Mutation" = "#8B58A4",
  "High_Mutation" = "#E3170D"
)

high_shm_genes <- shm_gene_table %>%
  filter(.data$set_name == "High_Mutation") %>%
  pull(.data$gene)
low_shm_genes <- shm_gene_table %>%
  filter(.data$set_name == "Low_Mutation") %>%
  pull(.data$gene)
ordered_shm_genes <- shm_gene_table$gene

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
    stop("Expected exactly one row in figure5_inputs.csv for sample: ", sample_id, call. = FALSE)
  }
  row
}

read_csv_if_exists <- function(path, ...) {
  full_path <- if (file.exists(path)) path else file.path(repo_root, path)
  if (!file.exists(full_path)) {
    stop("Input file not found: ", full_path, call. = FALSE)
  }
  read.csv(full_path, ...)
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

load_sample_object <- function(sample_id, filtered = TRUE) {
  cfg <- sample_row(sample_id)
  path <- if (filtered && file.exists(cfg$filtered_object_rds)) cfg$filtered_object_rds else cfg$seurat_object_rds
  full_path <- if (file.exists(path)) path else file.path(repo_root, path)
  if (!file.exists(full_path)) {
    stop("Seurat object not found: ", full_path, call. = FALSE)
  }
  list(cfg = cfg, object = readRDS(full_path))
}

mutation_rate_columns <- c(
  "allVAlignments_mutation_rate", "allDAlignments_mutation_rate",
  "allJAlignments_mutation_rate", "allCAlignments_mutation_rate",
  "VDJ_mutation_rate", "VDJC_mutation_rate", "VJ_mutation_rate"
)

prepare_spot_mutation_metadata <- function(seurat_obj, cfg, focus_region = cfg$shm_focus_region) {
  alignments <- read_csv_if_exists(cfg$igh_alignment_csv, stringsAsFactors = FALSE) %>%
    distinct()
  if (!cfg$spot_column %in% colnames(alignments)) {
    stop("Spot column not found in IGH alignment table: ", cfg$spot_column, call. = FALSE)
  }

  present_rate_cols <- intersect(mutation_rate_columns, colnames(alignments))
  alignments <- alignments %>%
    mutate(across(all_of(present_rate_cols), as.numeric))

  mutation_avgs <- alignments %>%
    group_by(.data[[cfg$spot_column]]) %>%
    summarise(
      avg_15G_VDJ_mutation_rate = mean(.data$VDJ_mutation_rate, na.rm = TRUE),
      avg_15G_VDJC_mutation_rate = mean(.data$VDJC_mutation_rate, na.rm = TRUE),
      avg_15G_VJ_mutation_rate = mean(.data$VJ_mutation_rate, na.rm = TRUE),
      .groups = "drop"
    )
  colnames(mutation_avgs)[1] <- "spatial_barcode_ori"

  meta <- seurat_obj@meta.data %>%
    mutate(spatial_barcode_ori = rownames(seurat_obj@meta.data)) %>%
    left_join(mutation_avgs, by = "spatial_barcode_ori")

  if (!is.na(focus_region) && nzchar(focus_region) && cfg$region_column %in% colnames(meta)) {
    meta <- meta %>% filter(.data[[cfg$region_column]] == focus_region)
  }

  meta %>%
    mutate(
      SHM_VDJ = case_when(
        is.na(.data$avg_15G_VDJ_mutation_rate) ~ "No_Mutation",
        .data$avg_15G_VDJ_mutation_rate <= 0 ~ "No_Mutation",
        .data$avg_15G_VDJ_mutation_rate <= 0.05 ~ "Low_Mutation",
        TRUE ~ "High_Mutation"
      ),
      SHM_VDJ = factor(.data$SHM_VDJ, levels = shm_levels)
    )
}

attach_shm_group <- function(seurat_obj, shm_meta) {
  seurat_obj$SHM_VDJ <- NA_character_
  shared <- intersect(rownames(seurat_obj@meta.data), shm_meta$spatial_barcode_ori)
  seurat_obj$SHM_VDJ[shared] <- as.character(shm_meta$SHM_VDJ[match(shared, shm_meta$spatial_barcode_ori)])
  seurat_obj$SHM_VDJ[is.na(seurat_obj$SHM_VDJ)] <- "Other"
  seurat_obj$SHM_VDJ <- factor(seurat_obj$SHM_VDJ, levels = c(shm_levels, "Other"))
  seurat_obj
}

plot_ranked_shm <- function(shm_meta, sample_label) {
  df <- shm_meta %>%
    filter(!is.na(.data$avg_15G_VDJ_mutation_rate)) %>%
    arrange(desc(.data$avg_15G_VDJ_mutation_rate)) %>%
    mutate(rank = row_number())
  vline_rank <- df %>%
    filter(.data$avg_15G_VDJ_mutation_rate < 0.05) %>%
    slice(1) %>%
    pull(.data$rank)
  if (!length(vline_rank)) vline_rank <- NA_integer_

  ggplot(df, aes(x = .data$rank, y = .data$avg_15G_VDJ_mutation_rate, color = .data$SHM_VDJ)) +
    geom_point(size = 1.4, alpha = 0.85) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "#B07AA1") +
    {if (!is.na(vline_rank)) geom_vline(xintercept = vline_rank, linetype = "dashed", color = "#B07AA1")} +
    {if (!is.na(vline_rank)) annotate("text", x = vline_rank + 5, y = 0.05, label = vline_rank, hjust = 0, vjust = -0.5, color = "#4B0082", size = 3)} +
    scale_color_manual(values = shm_colors, drop = FALSE) +
    labs(title = "Average VDJ Mutation Rate", x = "Ranked Spots", y = NULL, color = NULL) +
    theme_classic()
}

add_shm_module_scores <- function(seurat_obj, assay = "SCT", nbin = 2) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for AddModuleScore.", call. = FALSE)
  }
  Seurat::DefaultAssay(seurat_obj) <- assay
  high_genes <- intersect(high_shm_genes, rownames(seurat_obj))
  low_genes <- intersect(low_shm_genes, rownames(seurat_obj))
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = list(high_genes), name = "Score_High", nbin = nbin)
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = list(low_genes), name = "Score_Low", nbin = nbin)
  seurat_obj
}

plot_module_correlation <- function(seurat_obj, shm_meta, score_col, color_value, x_label) {
  score_df <- Seurat::FetchData(seurat_obj, vars = score_col) %>%
    mutate(spatial_barcode_ori = rownames(.)) %>%
    left_join(shm_meta %>% select(.data$spatial_barcode_ori, .data$avg_15G_VDJ_mutation_rate), by = "spatial_barcode_ori") %>%
    filter(!is.na(.data[[score_col]]), !is.na(.data$avg_15G_VDJ_mutation_rate))

  test <- suppressWarnings(cor.test(score_df[[score_col]], score_df$avg_15G_VDJ_mutation_rate, method = "spearman"))
  label <- paste0("R = ", signif(unname(test$estimate), 3), "\nP = ", format.pval(test$p.value, digits = 3))

  ggplot(score_df, aes(x = .data[[score_col]], y = .data$avg_15G_VDJ_mutation_rate)) +
    geom_point(size = 0.8, alpha = 0.75, color = color_value) +
    geom_smooth(method = "lm", se = TRUE, color = color_value, linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.3, size = 3) +
    labs(x = x_label, y = "SHM rate") +
    theme_classic()
}

star_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

celltype_columns <- function(meta_data, prefix = "q05cell_abundance_w_sf_") {
  grep(paste0("^", prefix), colnames(meta_data), value = TRUE)
}

clean_celltype_name <- function(x, prefix = "q05cell_abundance_w_sf_") {
  gsub(paste0("^", prefix), "", x)
}

prepare_celltype_fraction <- function(meta_data, shm_meta, prefix = "q05cell_abundance_w_sf_") {
  cols <- celltype_columns(meta_data, prefix)
  meta_data %>%
    mutate(spatial_barcode_ori = rownames(meta_data)) %>%
    left_join(shm_meta %>% select(.data$spatial_barcode_ori, .data$SHM_VDJ), by = "spatial_barcode_ori") %>%
    filter(.data$SHM_VDJ %in% shm_levels) %>%
    select(.data$spatial_barcode_ori, .data$SHM_VDJ, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "CellType", values_to = "Abundance") %>%
    mutate(CellType = clean_celltype_name(.data$CellType)) %>%
    group_by(.data$spatial_barcode_ori) %>%
    mutate(TotalAbundance = sum(.data$Abundance, na.rm = TRUE),
           Fraction = ifelse(.data$TotalAbundance > 0, .data$Abundance / .data$TotalAbundance * 100, NA_real_)) %>%
    ungroup() %>%
    mutate(CellCategory = case_when(
      .data$CellType %in% b_cell_types ~ "B cell",
      .data$CellType %in% t_cell_types ~ "T cell",
      TRUE ~ "Other cell"
    ))
}

plot_celltype_trend <- function(cell_frac, category, group_levels = shm_levels) {
  dat <- cell_frac %>%
    filter(.data$CellCategory == category, .data$SHM_VDJ %in% group_levels) %>%
    mutate(SHM_VDJ = factor(.data$SHM_VDJ, levels = group_levels))
  means_df <- dat %>%
    group_by(.data$SHM_VDJ, .data$CellType) %>%
    summarise(MeanFrac = mean(.data$Fraction, na.rm = TRUE), .groups = "drop")
  pair_df <- map_dfr(unique(dat$CellType), function(cell_type) {
    sub <- dat %>% filter(.data$CellType == cell_type)
    present_groups <- group_levels[group_levels %in% as.character(unique(sub$SHM_VDJ))]
    if (length(present_groups) < 2) return(data.frame())
    map_dfr(seq_along(combn(present_groups, 2, simplify = FALSE)), function(i) {
      pair <- combn(present_groups, 2, simplify = FALSE)[[i]]
      x <- (sub %>% filter(.data$SHM_VDJ == pair[1]))$Abundance
      y <- (sub %>% filter(.data$SHM_VDJ == pair[2]))$Abundance
      p <- tryCatch(suppressWarnings(wilcox.test(x, y)$p.value), error = function(e) NA_real_)
      target_y <- means_df %>%
        filter(.data$CellType == cell_type, .data$SHM_VDJ == pair[2]) %>%
        pull(.data$MeanFrac)
      data.frame(
        CellType = cell_type,
        SHM_VDJ = pair[2],
        comparison = paste(pair, collapse = "_vs_"),
        p = p,
        star = star_label(p),
        y = ifelse(length(target_y), target_y, 0) + max(means_df$MeanFrac, na.rm = TRUE) * (0.05 + 0.04 * (i - 1)),
        stringsAsFactors = FALSE
      )
    })
  }) %>%
    filter(.data$star != "")
  colors <- setNames(rep(cell_type_palette, length.out = length(unique(dat$CellType))), sort(unique(dat$CellType)))

  ggplot(means_df, aes(x = .data$SHM_VDJ, y = .data$MeanFrac, group = .data$CellType, color = .data$CellType)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    geom_point(size = 1.6) +
    geom_text(data = pair_df, aes(y = .data$y, label = .data$star), show.legend = FALSE, size = 2.6) +
    scale_color_manual(values = colors) +
    labs(title = paste("Cell Type Trend -", category), x = NULL, y = "Relative Abundance (%)", color = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
}

get_assay_matrix <- function(seurat_obj, assay_name) {
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop("Assay not found: ", assay_name, call. = FALSE)
  }
  as.matrix(seurat_obj@assays[[assay_name]]@data)
}

calc_abundance_totals <- function(meta_data, prefix = "q05cell_abundance_w_sf_") {
  t_cols <- intersect(paste0(prefix, t_cell_types), colnames(meta_data))
  b_cols <- intersect(paste0(prefix, b_cell_types), colnames(meta_data))
  meta_data %>%
    mutate(Spot = rownames(meta_data)) %>%
    rowwise() %>%
    mutate(
      T_total = if (length(t_cols)) sum(c_across(all_of(t_cols)), na.rm = TRUE) else NA_real_,
      B_total = if (length(b_cols)) sum(c_across(all_of(b_cols)), na.rm = TRUE) else NA_real_
    ) %>%
    ungroup()
}

calc_diversity_norm <- function(seurat_obj, assay_name, meta_df, norm_col) {
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("The stringdist package is required for Levenshtein distance.", call. = FALSE)
  }
  mat <- get_assay_matrix(seurat_obj, assay_name)
  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(long) <- c("Clonotype", "Spot", "UMI")
  long <- long %>% filter(.data$UMI > 0)

  by_spot <- long %>%
    group_by(.data$Spot, .data$Clonotype) %>%
    summarise(UMIs = sum(.data$UMI), .groups = "drop")

  by_spot %>%
    group_by(.data$Spot) %>%
    summarise(
      Shannon = {
        total <- sum(.data$UMIs)
        p <- .data$UMIs / total
        -sum(p * log(p), na.rm = TRUE)
      },
      MeanLevenshtein = {
        seqs <- unique(.data$Clonotype)
        if (length(seqs) > 1) {
          dm <- stringdist::stringdistmatrix(seqs, seqs, method = "lv")
          mean(dm[upper.tri(dm)], na.rm = TRUE)
        } else {
          NA_real_
        }
      },
      UniqueClonotypes = n_distinct(.data$Clonotype),
      .groups = "drop"
    ) %>%
    left_join(meta_df %>% select(.data$Spot, all_of(norm_col)), by = "Spot") %>%
    rename(NormFactor = all_of(norm_col)) %>%
    mutate(
      Shannon_norm = ifelse(.data$NormFactor > 0, .data$Shannon / .data$NormFactor, NA_real_),
      UniqueClonotypes_norm = ifelse(.data$NormFactor > 0, .data$UniqueClonotypes / .data$NormFactor, NA_real_)
    )
}

pairwise_wilcox_bh <- function(df, group_col, value_col) {
  df <- df %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]]))
  groups <- levels(factor(as.character(df[[group_col]])))
  if (length(groups) < 2) return(data.frame())
  map_dfr(combn(groups, 2, simplify = FALSE), function(pair) {
    x <- (df %>% filter(.data[[group_col]] == pair[1]))[[value_col]]
    y <- (df %>% filter(.data[[group_col]] == pair[2]))[[value_col]]
    data.frame(group1 = pair[1], group2 = pair[2], p = suppressWarnings(wilcox.test(x, y)$p.value))
  }) %>%
    mutate(p_adj = p.adjust(.data$p, method = "BH"))
}

plot_violin_metric <- function(df, group_col, value_col, y_label, colors = shm_colors) {
  df <- df %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]]))
  comps <- pairwise_wilcox_bh(df, group_col, value_col)
  subtitle <- if (nrow(comps)) paste(paste0(comps$group1, " vs ", comps$group2, ": q=", signif(comps$p_adj, 3)), collapse = "; ") else NULL
  ggplot(df, aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, alpha = 0.9, scale = "width") +
    geom_boxplot(width = 0.11, outlier.size = 0.5, fill = "grey85", color = "black") +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(x = NULL, y = y_label, subtitle = subtitle) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

prepare_cgene_trend <- function(igh_alignment_csv, shm_meta, group_levels = shm_levels) {
  alignments <- read_csv_if_exists(igh_alignment_csv, stringsAsFactors = FALSE) %>%
    distinct() %>%
    filter(!is.na(.data$bestCGene), .data$bestCGene != "nan") %>%
    left_join(shm_meta %>% select(.data$spatial_barcode_ori, .data$SHM_VDJ), by = "spatial_barcode_ori") %>%
    filter(.data$SHM_VDJ %in% group_levels)
  umi_summary <- alignments %>%
    group_by(.data$spatial_barcode_ori, .data$bestCGene, .data$SHM_VDJ) %>%
    summarise(UMI_count = n_distinct(.data$UMI), .groups = "drop")
  trend <- umi_summary %>%
    group_by(.data$bestCGene, .data$SHM_VDJ) %>%
    summarise(total_UMI = sum(.data$UMI_count), .groups = "drop") %>%
    group_by(.data$SHM_VDJ) %>%
    mutate(UMI_percent = 100 * .data$total_UMI / sum(.data$total_UMI)) %>%
    ungroup() %>%
    mutate(SHM_VDJ = factor(.data$SHM_VDJ, levels = group_levels))
  list(raw = umi_summary, trend = trend)
}

plot_cgene_trend <- function(cgene_data) {
  genes <- sort(unique(cgene_data$trend$bestCGene))
  colors <- setNames(rep(cell_type_palette, length.out = length(genes)), genes)
  ggplot(cgene_data$trend, aes(x = .data$SHM_VDJ, y = .data$UMI_percent, group = .data$bestCGene, color = .data$bestCGene)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = colors) +
    labs(x = NULL, y = "Relative Abundance (%)", color = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
