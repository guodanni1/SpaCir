# Figure 4 and Figure S8-S10 Code

This folder contains cleaned, panel-oriented R code for germinal-center maturity analyses in Figure 4 and Figure S8-S10.

## Sample Labels

Sample names, object paths, assay names, and table paths are configured in `config/figure4_inputs.csv`.
The main Figure 4 scripts use the `Tonsil2` row by default. Supplementary panels reuse the same functions with the sample-specific rows listed in `config/supplementary_panel_parameters.csv`.

## Scripts By Panel

- `config/panel_map.csv`: panel-to-script mapping for Figure 4 and Figure S8-S10.
- `config/supplementary_panel_parameters.csv`: supplement-specific sample IDs, GC lists, and reused plotting functions.
- `config/gc_maturity_parameters.csv`: EarlyGC/LateGC marker genes and AddModuleScore rule.
- `config/clone_class_parameters.csv`: Figure 4B A-D clone-class definitions and the 95th percentile UMI rule.
- `config/hdwgcna_go_parameters.csv`: hdWGCNA, module, hub-gene, and GO enrichment parameters.
- `R/00_setup.R`: shared settings, marker lists, color palettes, plotting utilities, GC maturity scoring, clonotype class assignment, and normalized diversity helpers.
- `R/figure4A_D_gc_maturity_annotation.R`: Figure 4A/D and supplementary GC spatial maps, maturity dotplots, C-gene bars, and class bars.
- `R/figure4B_C_G_H_I_J_K_L_clonotype_maturity.R`: Figure 4B/C/G/H/I/J/K/L and related Figure S9/S10 clonotype, SHM, normalized diversity, cumulative-frequency, and class/isotype trend panels.
- `R/figure4E_F_S8D_E_wgcna_go.R`: Figure 4E/F and Figure S8D/E WGCNA module, hub-gene enrichment, and GO-term panels.
- `R/figure4M_N_O_celltype_trends.R`: Figure 4M/N/O and supplementary B/T/other cell-type trend panels with EarlyGC versus LateGC significance marks.
- `R/figureS8_supplement.R`: Figure S8 panel-level entry points.
- `R/figureS9_supplement.R`: Figure S9 panel-level entry points.
- `R/figureS10_supplement.R`: Figure S10 panel-level entry points.

## Notes

The public scripts keep only the EarlyGC and LateGC maturity states used in the final analysis. Diversity panels use abundance-normalized metrics. Trend panels use the significance-marked line-plot functions from `R/00_setup.R`.

## Explicit Parameters

The key Figure 4 definitions are not hidden inside helpers:

- GC maturity uses `CD83;CD86;CXCR5;GPR183;SLAMF1;BCL6` for EarlyGC and `IRF4;PRDM1;XBP1;ZBTB20;FOXP3;DUSP2;IRF8;GADD45B;JCHAIN;TOX2;SDC1` for LateGC. The module-score `nbin` is 10, and the higher mean score within each GC part defines the EarlyGC/LateGC label.
- Clone classes use private/shared GC-site detection and the sample-specific 95th percentile of total IGH UMI. A = private expanded, B = shared expanded, C = private unexpanded, and D = shared unexpanded.
- hdWGCNA and GO selections are documented in `config/hdwgcna_go_parameters.csv`; final line plots use the significance-marked helper functions and normalized diversity only.

## Running

Fill in the paths in `config/figure4_inputs.csv`, then source the required panel script from the repository root. For example:

```r
source("analysis/figure4/R/figure4B_C_G_H_I_J_K_L_clonotype_maturity.R")
run_figure4_clonotype_panels("Tonsil2")
```

Outputs are written to `analysis/figure4/outputs/`.
