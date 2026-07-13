# Figure 3 and Figure S5-S7 Code

This folder contains cleaned, panel-oriented code for Figure 3 and Figure S5-S7.

## Sample Labels

Sample names used in figures and output files are defined in `config/figure3_inputs.csv` and applied by `R/00_setup.R`.

## Scripts By Panel

- `config/panel_map.csv`: detailed panel-to-script mapping for Figure 3 and Figure S5-S7.
- `config/supplementary_panel_parameters.csv`: supplement-specific sample IDs, marker lists, and reused code.
- `config/sample_specific_parameters.csv`: per-sample Monocle3 and pseudotime parameters.
- `config/pseudotime_gene_sets.csv`: genes plotted in Figure 3C and Figure S5A.
- `config/clonefamily_parameters.csv`: IGH clonefamily and TRB cluster thresholds and controls.
- `config/branch_arrow_parameters.csv`: sample-specific branch/arrow panels and retained branch types.
- `R/figure3A_B_pseudotime_maps.R`: Figure 3A/B Monocle3 pseudotime UMAP and spatial maps.
- `R/figure3C_S5A_pseudotime_gene_trends.R`: Figure 3C and Figure S5A pseudotime gene-trend panels.
- `R/figure3D_E_F_G_clonefamily_summary.R`: Figure 3D-G and related IGH clonefamily summary panels.
- `R/figure3H_I_spatial_clonefamily_branches.R`: Figure 3H/I spatial clonefamily maps and branch arrows.
- `R/figure3J_K_S5E_S6F_G_S7F_G_csr_sankey.R`: Figure 3J/K and supplementary Sankey panels.
- `R/figure3L_M_S6H_I_S7H_I_clonal_similarity.R`: Figure 3L/M and supplementary Jaccard/shared-clone panels.
- `R/figure3N_O_S6J_K_S7J_intergc_distance.R`: Figure 3N/O and supplementary InterGC count/distance panels.
- `R/figure3P_Q_R_S6O_P_Q_S7N_O_P_trb_pathology.R`: Figure 3P-R and supplementary TRB pathology panels.
- `R/figureS5_clonefamily_branch_detail.R`: Figure S5B-D selected clonefamily sequence-logo/tree-arrow details.
- `R/figureS6_LN2_supplement.R`: Figure S6 panel-level entry points for LN2.
- `R/figureS7_Tonsil1_supplement.R`: Figure S7 panel-level entry points for Tonsil1.

## Notes

The Figure 3 workflow has two main analysis branches:

- IGH clonefamily branch: Monocle3 pseudotime, IGH clonefamily clustering, CSR network, CDR3 length, nearest-neighbor distance, clonal-tree branch arrows, Sankey summaries, and Jaccard/shared-clone analyses.
- TRB cluster branch: GLIPH/TRB cluster distance, pathology annotation using VDJdb/McPAS-style reference tables, pathology composition, and spatial maps of annotated TRB clusters.

## Explicit Parameters

Figure 3 uses different sample-specific clonefamily and region annotations. The public code therefore keeps the panel logic in R scripts and the reproducibility-critical values in config tables:

- `config/sample_specific_parameters.csv` records the Monocle3 q-value cutoff, gene-module resolution, and sample-specific clonefamily count shown in each panel set.
- `config/pseudotime_gene_sets.csv` records the genes used for Figure 3C and Figure S5A pseudotime trend panels.
- `config/clonefamily_parameters.csv` records the IGH `lv_distance2 <= 4` clonefamily filter, TRB cluster control logic, and CSR node set.
- `config/branch_arrow_parameters.csv` records which branch-arrow summaries are retained for the main and supplementary panels.

## Running

Fill in the paths in `config/figure3_inputs.csv`, then source the panel script from the repository root. For example:

```r
source("analysis/figure3/R/figure3D_E_F_G_clonefamily_summary.R")
plot_clonefamily_size("Tonsil2")
```

Outputs are written to `analysis/figure3/outputs/`.
