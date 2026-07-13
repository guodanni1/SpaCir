# Shared settings for SpaCir Figure 2 and Figure S2-S4 downstream plotting.
# Run scripts from the repository root, or set SPACIR_ROOT before sourcing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

repo_root <- Sys.getenv("SPACIR_ROOT", unset = getwd())
figure2_dir <- file.path(repo_root, "analysis", "figure2")
input_config <- file.path(figure2_dir, "config", "figure2_inputs.csv")
output_dir <- file.path(figure2_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_config <- read.csv(input_config, stringsAsFactors = FALSE)

# Display labels used in figures and output filenames.
sample_label_map <- setNames(sample_config$display_sample, sample_config$sample_id)
object_label_map <- setNames(sample_config$display_sample, sample_config$object_name)

display_sample <- function(x) {
  y <- as.character(x)
  mapped_sample <- sample_label_map[y]
  mapped_object <- object_label_map[y]
  out <- ifelse(!is.na(mapped_sample), mapped_sample, y)
  out <- ifelse(!is.na(mapped_object), mapped_object, out)
  out
}

sample_levels <- c("LN1", "Tonsil2", "Tonsil1", "LN2")
sample_colors <- c(
  "LN1" = "#CB5D17",
  "Tonsil2" = "#823379",
  "Tonsil1" = "#E9CC54",
  "LN2" = "#0C6CA7"
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Region colors follow the cleaned Figure 2 palette. Names are reused by
# Figure S2-S4; scripts subset the vector to the regions present in each sample.
region_colors <- c(
  "Follicular DZ" = "#E4572E",
  "Follicular LZ" = "#E9C94F",
  "Inter LZ-DZ" = "#9BCB66",
  "Mantle zone" = "#F29B7A",
  "Vascular" = "#F3E999",
  "Vascular/Mesenchymal" = "#C46693",
  "Muscle/Mesenchymal" = "#F0B5D1",
  "Inter-Follicular CD5+" = "#7C66A2",
  "Inter-Follicular CD5-" = "#6FC7B6",
  "Inter-Follicular SDC1+" = "#68BFAE",
  "Inter-Follicular low SDC1+" = "#B8D989",
  "Inter-Follicular IL7R+" = "#BCD58A",
  "Inter-Follicular CXCL10+" = "#8A6FA6",
  "Epithelial/Muscle" = "#9A7466",
  "Subepithelial/Mesenchymal/Vascular" = "#6C6E70",
  "Subepithelial/Mesenchymal" = "#817070",
  "Epithelium-like structures" = "#F0B3C7",
  "Inter-Follicular CD6+" = "#7B669A",
  "Inter-Follicular XBP1+" = "#8DD3C7"
)

cell_type_colors <- c(
  "B Cycling" = "#C9655B",
  "B preGC" = "#E5A17B",
  "B mem" = "#80B48C",
  "B naive" = "#4C89B8",
  "B GC DZ" = "#E4572E",
  "B GC LZ" = "#E9C94F",
  "T CD4+ naive" = "#C9C9A5",
  "Endo" = "#9BBFA3",
  "Mast" = "#D84C8A",
  "B plasma" = "#DCA338",
  "FDC" = "#B9B7A6",
  "Macrophages M1" = "#50A87A",
  "NK" = "#7C7C7C",
  "B IFN" = "#9DD5C7",
  "T CD4+ TfH GC" = "#F08C3A",
  "T TIM3+" = "#4E6E91",
  "DC CCR7+" = "#8ED2E8",
  "DC cDC1" = "#B783C8",
  "ILC" = "#9B7BB0",
  "B activated" = "#97A0C3",
  "Monocytes" = "#E58D84",
  "DC cDC2" = "#D890A4",
  "T CD8+ naive" = "#A8C4B0",
  "DC pDC" = "#82B4B2",
  "Macrophages M2" = "#4F7A63",
  "T Tfr" = "#4F4F4F",
  "B GC prePB" = "#C9B064",
  "T CD4+" = "#9A7AB3",
  "T CD8+ cytotoxic" = "#6B7EA7",
  "T Treg" = "#C489B5",
  "T CD8+ CD161+" = "#F0B0AC",
  "NKT" = "#9BD481"
)

read_csv_if_exists <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  read.csv(path, ...)
}

sample_row <- function(sample_id) {
  row <- sample_config[sample_config$sample_id == sample_id | sample_config$object_name == sample_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one row in figure2_inputs.csv for sample: ", sample_id, call. = FALSE)
  }
  row
}

save_panel <- function(plot, filename, width, height) {
  out <- file.path(output_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved: ", out)
}

clean_celltype_name <- function(x) {
  x %>%
    str_remove("^q05cell_abundance_w_sf_") %>%
    str_replace_all("[._]", " ") %>%
    str_squish()
}

cell_abundance_long <- function(meta_df) {
  abundance_cols <- grep("^q05cell_abundance_w_sf_", names(meta_df), value = TRUE)
  if (length(abundance_cols) == 0) {
    stop("No cell2location abundance columns found. Expected q05cell_abundance_w_sf_*.", call. = FALSE)
  }
  meta_df %>%
    mutate(Spot = rownames(meta_df)) %>%
    pivot_longer(all_of(abundance_cols), names_to = "CellType", values_to = "Abundance") %>%
    mutate(CellType = clean_celltype_name(CellType))
}

split_braced <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  strsplit(gsub("^\\{|\\}$", "", x), "\\}\\{", fixed = FALSE)[[1]]
}

# UMI-weighted Shannon diversity used throughout the Word code:
# collapse identical clonotypes within a spot, sum their UMI counts, then compute
# -sum(p * log(p)). This keeps high-UMI clonotypes weighted as in the original.
weighted_shannon <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  total <- sum(counts)
  if (length(counts) == 0 || total == 0) {
    return(NA_real_)
  }
  p <- counts / total
  -sum(p * log(p))
}

vdj_spot_metrics <- function(vdj_df, spot_col = "spatial_barcode_ori",
                             clonotype_col = "aaSeqCDR3_most",
                             umi_col = "UMI_count",
                             levenshtein_col = "MeanLevenshtein") {
  required <- c(spot_col, clonotype_col)
  missing <- setdiff(required, names(vdj_df))
  if (length(missing) > 0) {
    stop("Missing VDJ columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!umi_col %in% names(vdj_df)) {
    vdj_df[[umi_col]] <- 1
  }
  if (!levenshtein_col %in% names(vdj_df)) {
    vdj_df[[levenshtein_col]] <- NA_real_
  }
  vdj_df %>%
    transmute(
      Spot = .data[[spot_col]],
      Clonotype = .data[[clonotype_col]],
      UMI = suppressWarnings(as.numeric(.data[[umi_col]])),
      MeanLevenshtein = suppressWarnings(as.numeric(.data[[levenshtein_col]]))
    ) %>%
    filter(!is.na(Spot), !is.na(Clonotype)) %>%
    group_by(Spot, Clonotype) %>%
    summarise(UMI = sum(replace_na(UMI, 1)), MeanLevenshtein = mean(MeanLevenshtein, na.rm = TRUE), .groups = "drop") %>%
    group_by(Spot) %>%
    summarise(
      Shannon = weighted_shannon(UMI),
      UniqueClonotypes = n_distinct(Clonotype),
      MeanLevenshtein = mean(MeanLevenshtein, na.rm = TRUE),
      .groups = "drop"
    )
}
