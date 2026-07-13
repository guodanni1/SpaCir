# Add VDJ clonotype matrices as Seurat assays.
# This is a preparation step used by Figure 1E and Figure 1I-N/S1G.

source("analysis/figure1/R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat)
})

add_clonotype_assay <- function(obj, matrix_csv, assay_name) {
  # Input matrix rows are spatial barcodes and columns are clonotypes.
  # Seurat assays require features x spots, so the matrix is completed to all
  # spots in the object and then transposed before CreateAssayObject().
  read_count_matrix <- read_csv_if_exists(matrix_csv, row.names = 1)
  all_spatial_bc <- rownames(obj@meta.data)

  complete_matrix <- matrix(0, nrow = length(all_spatial_bc), ncol = ncol(read_count_matrix))
  rownames(complete_matrix) <- all_spatial_bc
  colnames(complete_matrix) <- colnames(read_count_matrix)

  existing_barcodes <- intersect(rownames(read_count_matrix), all_spatial_bc)
  complete_matrix[existing_barcodes, ] <- as.matrix(read_count_matrix[existing_barcodes, ])

  sparse_matrix <- as(t(complete_matrix), "dgCMatrix")
  obj[[assay_name]] <- CreateAssayObject(counts = sparse_matrix)
  obj
}

for (i in seq_len(nrow(sample_config))) {
  object_name <- sample_config$object_name[i]
  if (!exists(object_name)) {
    message("Skip missing Seurat object: ", object_name)
    next
  }

  obj <- get(object_name)
  obj <- add_clonotype_assay(obj, sample_config$trb_matrix_csv[i], "aaSeqCDR3_Lv1_UMI3_TRB_20G")
  obj <- add_clonotype_assay(obj, sample_config$igh_matrix_csv[i], "aaSeqCDR3_Lv1_UMI3_IGH_20G")
  assign(object_name, obj, envir = .GlobalEnv)

  message("Added TRB/IGH clonotype assays to ", object_name, " (display label: ",
          sample_config$display_sample[i], ")")
}
