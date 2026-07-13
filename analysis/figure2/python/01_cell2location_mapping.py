#!/usr/bin/env python3
"""Run cell2location mapping for Figure 2 cell-type abundance estimation.

The script exposes all paths as command-line arguments so the workflow can be
run without embedding local project directories in the source code.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run cell2location spatial mapping.")
    parser.add_argument("--spatial-h5ad", required=True, help="Spatial AnnData file.")
    parser.add_argument("--reference-model", required=True, help="Trained RegressionModel directory.")
    parser.add_argument("--output-dir", required=True, help="Directory for mapped AnnData and plots.")
    parser.add_argument("--max-epochs", type=int, default=30000, help="Maximum training epochs.")
    parser.add_argument("--n-cells-per-location", type=int, default=30, help="Expected cells per Visium spot.")
    parser.add_argument("--detection-alpha", type=float, default=20.0, help="Cell2location detection prior.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import scanpy as sc
    import cell2location

    adata_vis = sc.read_h5ad(args.spatial_h5ad)

    # Load the reference regression model and extract inferred cell-type
    # signatures from the trained cell2location regression model.
    adata_ref = adata_vis.copy()
    regression_model = cell2location.models.RegressionModel.load(args.reference_model, adata_ref)
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"]
    else:
        inf_aver = adata_ref.var[[c for c in adata_ref.var.columns if "means_per_cluster_mu_fg" in c]]

    intersect = adata_vis.var_names.intersection(inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :]

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    # These priors are exposed as CLI parameters because they are
    # sample/platform dependent.
    model = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
    )
    model.train(max_epochs=args.max_epochs, batch_size=None, train_size=1)
    adata_vis = model.export_posterior(adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": None})

    adata_vis.write(out_dir / "cell2location_mapped.h5ad")
    adata_vis.obsm["q05_cell_abundance_w_sf"].to_csv(out_dir / "cell2location_abundance_q05.csv")
    adata_vis.obs.to_csv(out_dir / "cell2location_obs_metadata.csv")


if __name__ == "__main__":
    main()
