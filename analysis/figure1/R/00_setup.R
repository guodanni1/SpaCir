# Shared settings for SpaCir Figure 1 and Figure S1 downstream plotting.
# Run scripts from the repository root, or set SPACIR_ROOT before sourcing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(Matrix)
})

repo_root <- Sys.getenv("SPACIR_ROOT", unset = getwd())
figure1_dir <- file.path(repo_root, "analysis", "figure1")
input_config <- file.path(figure1_dir, "config", "figure1_inputs.csv")
output_dir <- file.path(figure1_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_config <- read.csv(input_config, stringsAsFactors = FALSE)

# Display labels used in figures and output filenames.
sample_label_map <- setNames(sample_config$display_sample, sample_config$sample_id)
object_label_map <- setNames(sample_config$display_sample, sample_config$object_name)

display_sample <- function(x) {
  y <- as.character(x)
  mapped <- sample_label_map[y]
  mapped_obj <- object_label_map[y]
  out <- ifelse(!is.na(mapped), mapped, y)
  out <- ifelse(!is.na(mapped_obj), mapped_obj, out)
  out
}

sample_levels <- c("LN1", "Tonsil2", "Tonsil1", "LN2")
sample_colors <- c(
  "LN1" = "#CB5D17",
  "Tonsil2" = "#823379",
  "Tonsil1" = "#E9CC54",
  "LN2" = "#0C6CA7"
)

read_csv_if_exists <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  read.csv(path, ...)
}

save_panel <- function(plot, filename, width, height) {
  out <- file.path(output_dir, filename)
  ggsave(out, plot, width = width, height = height)
  message("Saved: ", out)
}
