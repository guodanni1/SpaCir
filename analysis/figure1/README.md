# Figure 1 and Figure S1 R Code

This folder contains cleaned, panel-oriented R scripts for Figure 1 and Figure S1.

## Sample Labels

Sample names used in figures and output files are defined in `config/figure1_inputs.csv` and applied by `R/00_setup.R`.

## Scripts By Panel

- `R/01_add_clonotype_assays.R`: helper step used before plotting; adds TRB/IGH clonotype matrices into Seurat objects.
- `R/figure1C_correlation_gene_expression_vs_umi.R`: Figure 1C, correlation between spatial gene expression and VDJ UMI counts.
- `R/figure1D_E_clonotype_density.R`: Figure 1D table values and Figure 1E clonotype density per mm2.
- `R/figure1F_G_gene_usage.R`: Figure 1F TRB J gene usage and Figure 1G IGH J gene usage.
- `R/figure1H_shm_rate.R`: Figure 1H pooled IGH SHM rate violin plot.
- `R/figure1I_L_clone_size_composition.R`: Figure 1I and 1L clone size composition.
- `R/figure1J_K_M_N_spatial_expansion.R`: Figure 1J, 1K, 1M, and 1N spatial expansion and nearest-neighbor distance plots.
- `R/figureS1C_coverage_curves.R`: Figure S1C IGH/TRB sample coverage curves.
- `R/figureS1D_E_vgene_usage.R`: Figure S1D/S1E V gene usage panels.
- `R/figureS1F_G_per_sample_shm_distance.R`: Figure S1F/S1G per-sample SHM and clone-distance panels.

Outputs are written to `analysis/figure1/outputs/`.
