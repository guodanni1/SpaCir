# Figure S6: LN2 supplementary clonefamily and TRB panels.

source(file.path(Sys.getenv("SPACIR_ROOT", unset = getwd()), "analysis", "figure3", "R", "00_setup.R"))
source(file.path(figure3_dir, "R", "figure3D_E_F_G_clonefamily_summary.R"))
source(file.path(figure3_dir, "R", "figure3H_I_spatial_clonefamily_branches.R"))
source(file.path(figure3_dir, "R", "figure3J_K_S5E_S6F_G_S7F_G_csr_sankey.R"))
source(file.path(figure3_dir, "R", "figure3L_M_S6H_I_S7H_I_clonal_similarity.R"))
source(file.path(figure3_dir, "R", "figure3N_O_S6J_K_S7J_intergc_distance.R"))
source(file.path(figure3_dir, "R", "figure3P_Q_R_S6O_P_Q_S7N_O_P_trb_pathology.R"))

run_figureS6_clonefamily_summary <- function(sample_id = "LN2") {
  plot_clonefamily_size(sample_id)
  plot_clonefamily_csr(sample_id)
  plot_cdr3_length_distribution(sample_id)
  plot_clonefamily_distance(sample_id, output_tag = "figureS6D")
}

run_figureS6_spatial_similarity <- function(sample_id = "LN2", assay_name = "aaSeqCDR3_Lv1_UMI3_IGH_20G") {
  plot_clonefamily_spatial_map(sample_id, output_tag = "figureS6E")
  plot_transition_sankey(sample_id, mode = "region", output_tag = "figureS6F")
  plot_transition_sankey(sample_id, mode = "gc", output_tag = "figureS6G")
  presence <- clonotype_presence_from_seurat(sample_id, assay_name, "GC_all")
  plot_jaccard_tree(sample_id, presence, output_tag = "figureS6H")
  plot_shared_clone_bubble(sample_id, presence, output_tag = "figureS6I")
  plot_intergc_origin_counts(sample_id, output_tag = "figureS6J")
  plot_intergc_distance(sample_id, output_tag = "figureS6K")
}

run_figureS6_trb_pathology <- function(sample_id = "LN2") {
  run_trb_pathology_panels(sample_id, output_prefix = "figureS6")
}

# Examples:
# run_figureS6_clonefamily_summary()
# run_figureS6_spatial_similarity()
# run_figureS6_trb_pathology()
