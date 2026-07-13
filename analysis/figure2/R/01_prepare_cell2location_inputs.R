# Prepare Seurat-derived input files for the cell2location mapping step.
# This corresponds to the first "cell type annotation" block in the Word code.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

prepare_cell2location_inputs <- function(sample_id, assay = "Spatial", output_subdir = "cell2location_inputs") {
  cfg <- sample_row(sample_id)
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))
  DefaultAssay(obj) <- assay

  # cell2location expects a raw count matrix plus metadata. Matrix Market keeps
  # the count export sparse, which is important for Visium-scale objects.
  out_dir <- file.path(output_dir, output_subdir, cfg$sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  counts <- GetAssayData(obj, assay = assay, slot = "counts")
  Matrix::writeMM(counts, file.path(out_dir, "spatial_counts.mtx"))
  writeLines(rownames(counts), file.path(out_dir, "genes.tsv"))
  writeLines(colnames(counts), file.path(out_dir, "barcodes.tsv"))
  write.csv(obj@meta.data, file.path(out_dir, "metadata.csv"), quote = FALSE)

  message("Prepared cell2location inputs for ", display_sample(cfg$sample_id), ": ", out_dir)
  invisible(out_dir)
}

# Example:
# prepare_cell2location_inputs("Tonsil2")
