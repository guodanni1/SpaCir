# Figure 2L: CellChat interaction strength across tissue regions.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure2", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(CellChat)
})

run_cellchat_by_region <- function(sample_id, group_col = NULL, output_tag = "figure2L") {
  cfg <- sample_row(sample_id)
  group_col <- group_col %||% cfg$region_column
  obj <- readRDS(file.path(repo_root, cfg$seurat_rds))

  # CellChat uses the region labels as sender/receiver groups, matching the
  # Word-code block where createCellChat(..., group.by = "NewClusters") was used.
  cellchat <- createCellChat(object = obj, meta = obj@meta.data, group.by = group_col)
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  out_rds <- file.path(output_dir, paste0(output_tag, "_", display_sample(cfg$sample_id), "_cellchat.rds"))
  saveRDS(cellchat, out_rds)

  pdf(file.path(output_dir, paste0(output_tag, "_", display_sample(cfg$sample_id), "_interaction_strength.pdf")), width = 7, height = 4)
  netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "GnBu")
  dev.off()

  message("Saved CellChat object: ", out_rds)
  invisible(cellchat)
}

# Example:
# run_cellchat_by_region("Tonsil2")
