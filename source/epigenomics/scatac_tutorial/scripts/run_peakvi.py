#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp/learn_gencore_scatac_xdg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/learn_gencore_scatac_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/learn_gencore_scatac_numba")

import matplotlib.pyplot as plt

from peakvi_workflow import (
    add_qc_flags,
    build_peakvi_graph,
    choose_qc_thresholds,
    filter_for_peakvi,
    latent_qc_correlations,
    load_atac_dataset,
    plot_barcode_rank,
    plot_cluster_qc,
    plot_differential_accessibility,
    plot_fragment_length,
    plot_latent_qc_correlations,
    plot_peakvi_umap,
    plot_qc_distributions,
    plot_qc_relationships,
    plot_training_history,
    qc_summary,
    run_peakvi_differential_accessibility,
    save_figure,
    set_plot_style,
    train_peakvi,
    training_history_frame,
    write_cache_manifest,
    write_results,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the public PBMC scATAC tutorial through PEAKVI.")
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--max-peaks", type=int, default=50_000)
    parser.add_argument("--min-cells-per-peak", type=int, default=20)
    parser.add_argument("--resolution", type=float, default=0.6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    tutorial_dir = script_dir.parent
    repo_root = tutorial_dir.parents[2]
    raw_dir = tutorial_dir / "data" / "raw"
    processed_dir = tutorial_dir / "data" / "processed"
    results_dir = tutorial_dir / "results"
    figure_dir = repo_root / "source" / "img" / "epigenomics" / "scatac"
    model_dir = results_dir / "models" / "peakvi_pbmc"

    results_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    set_plot_style()

    adata, singlecell_all = load_atac_dataset(raw_dir, processed_dir)
    thresholds = choose_qc_thresholds(adata)
    add_qc_flags(adata, thresholds)
    qc_summary(adata).to_csv(results_dir / "qc_summary.tsv", sep="\t", index=False)

    figures = [
        ("01_barcode_rank.png", plot_barcode_rank(singlecell_all)),
        ("02_fragment_length.png", plot_fragment_length(processed_dir)),
        ("03_qc_distributions.png", plot_qc_distributions(adata, thresholds)),
        ("04_qc_relationships.png", plot_qc_relationships(adata)),
    ]
    for filename, fig in figures:
        save_figure(fig, figure_dir / filename)
        plt.close(fig)

    peakvi_adata = filter_for_peakvi(
        adata,
        min_cells_per_peak=args.min_cells_per_peak,
        max_peaks=None if args.max_peaks <= 0 else args.max_peaks,
    )
    write_cache_manifest(
        peakvi_adata,
        results_dir / "peakvi_cache_manifest.json",
        parameters={
            "max_epochs": args.max_epochs,
            "max_peaks": args.max_peaks,
            "min_cells_per_peak": args.min_cells_per_peak,
            "resolution": args.resolution,
            "n_latent": 15,
            "seed": 17,
        },
    )
    model = train_peakvi(
        peakvi_adata,
        model_dir=model_dir,
        max_epochs=args.max_epochs,
    )
    history = training_history_frame(model)
    history.to_csv(results_dir / "peakvi_training_history.tsv", sep="\t", index=False)
    build_peakvi_graph(peakvi_adata, resolution=args.resolution)
    correlations = latent_qc_correlations(peakvi_adata)
    group1 = "0"
    group2 = None
    da_results = run_peakvi_differential_accessibility(
        model,
        peakvi_adata,
        group1=group1,
        group2=group2,
        delta=0.05,
        fdr_target=0.05,
        n_samples_overall=5_000,
        seed=17,
    )
    da_results.to_csv(
        results_dir / f"peakvi_da_cluster_{group1}_vs_rest.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    figures = [
        ("05_peakvi_training.png", plot_training_history(history)),
        ("06_peakvi_latent_qc.png", plot_latent_qc_correlations(correlations)),
        ("07_peakvi_umap.png", plot_peakvi_umap(peakvi_adata)),
        ("08_peakvi_cluster_qc.png", plot_cluster_qc(peakvi_adata)),
        (
            "09_peakvi_differential_accessibility.png",
            plot_differential_accessibility(da_results, group1, group2),
        ),
    ]
    for filename, fig in figures:
        save_figure(fig, figure_dir / filename)
        plt.close(fig)

    write_results(peakvi_adata, thresholds, correlations, results_dir)
    peakvi_adata.write_h5ad(results_dir / "pbmc_peakvi.h5ad", compression="gzip")

    print(f"Input called nuclei: {adata.n_obs:,}")
    print(f"Nuclei retained for PEAKVI: {peakvi_adata.n_obs:,}")
    print(f"Peaks retained for PEAKVI: {peakvi_adata.n_vars:,}")
    print(f"PEAKVI clusters: {peakvi_adata.obs['leiden_peakvi'].nunique()}")
    print(
        f"Cluster {group1} candidate marker peaks: "
        f"{int(da_results['tutorial_marker'].sum()):,}"
    )
    print(f"Figures: {figure_dir}")
    print(f"Results: {results_dir}")


if __name__ == "__main__":
    main()
