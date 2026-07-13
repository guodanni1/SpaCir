# Figure 5 and Figure S11-S12 Code

This folder contains cleaned, panel-oriented code for Figure 5 and Figure S11-S12.

## Scripts By Panel

- `config/panel_map.csv`: panel-to-script mapping for Figure 5, Figure S11, and Figure S12.
- `config/figure5_inputs.csv`: sample-level object paths, assay names, and input tables.
- `config/shm_group_parameters.csv`: SHM group thresholds and diversity policy.
- `config/shm_gene_sets.csv`: High- and Low-Mutation SHM-associated genes and functional categories.
- `config/experimental_validation_parameters.csv`: input table expectations for validation panels.
- `R/00_setup.R`: shared helpers for SHM grouping, gene programs, cell-type trends, validation plots, and normalized clonotype diversity.
- `R/figure5A_B_S12A_B_shm_groups.R`: Figure 5A/B and Figure S12A/B ranked SHM and spatial SHM-group panels.
- `R/figure5C_D_S11_S12C_D_shm_gene_programs.R`: Figure 5C/D, Figure S11A, and Figure S12C/D SHM gene AUC and dotplot panels.
- `R/figure5E_F_S12E_F_module_correlations.R`: Figure 5E/F and Figure S12E/F module-score correlation panels.
- `R/figure5G_I_S12G_I_celltype_trends.R`: Figure 5G-I and Figure S12G-I cell-type trend panels.
- `R/figure5J_N_experimental_validation.R`: Figure 5J-N experimental workflow and validation panels.
- `R/figure5O_R_S12J_L_clonotype_norm_diversity.R`: Figure 5O-R and Figure S12J-L normalized clonotype diversity and C-gene panels.

## Key Parameters

SHM groups are assigned from `avg_15G_VDJ_mutation_rate`: `High_Mutation` is greater than 0.05, `Low_Mutation` is greater than 0 and less than or equal to 0.05, and `No_Mutation` is zero or missing. These thresholds are recorded in `config/shm_group_parameters.csv`.

The SHM-associated genes used for AUC ranking, dotplots, and module-score correlations are listed in `config/shm_gene_sets.csv`. The module-score panels use `Seurat::AddModuleScore(..., nbin = 2)`.

Shannon diversity panels export normalized Shannon only. IGH Shannon is normalized by B-cell abundance and TRB Shannon is normalized by T-cell abundance. Raw Shannon figure exports are intentionally omitted.

## Running

Fill in the public placeholder paths in `config/figure5_inputs.csv`, then source the required panel script from the repository root. For example:

```r
source("analysis/figure5/R/figure5A_B_S12A_B_shm_groups.R")
run_shm_group_panels("Tonsil2", output_tag = "figure5")
```

Outputs are written to `analysis/figure5/outputs/`.
