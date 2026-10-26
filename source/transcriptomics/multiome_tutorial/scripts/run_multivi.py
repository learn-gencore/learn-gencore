#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp/learn_gencore_multiome_xdg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/learn_gencore_multiome_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/learn_gencore_multiome_numba")
# On macOS, the pip wheels of PyTorch and scikit-learn each bundle an OpenMP
# runtime. In this environment, MultiVI training stalled after Scanpy's
# neighbor and UMAP steps unless OpenMP was limited to one thread.
os.environ.setdefault("OMP_NUM_THREADS", "1")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from multiome_workflow import (
    DC_MARKERS,
    POPULATION_MARKERS,
    T_STATE_MARKERS,
    add_qc_flags,
    assign_populations,
    atac_baseline,
    build_embedding_graph,
    build_multivi_input,
    check_observed_libraries,
    choose_qc_thresholds,
    cluster_qc_table,
    clustering_agreement,
    compare_rna_groups,
    latent_qc_correlations,
    latent_qc_correlations_within,
    load_multiome_dataset,
    load_multivi,
    choose_marker_loci,
    locus_accessibility,
    marker_locus_specificity,
    mitochondrial_link,
    marker_detection,
    marker_evidence,
    normalize_rna,
    plot_atac_qc,
    plot_baselines,
    plot_cluster_agreement,
    plot_cluster_qc,
    plot_embedding_comparison,
    plot_fragment_length,
    plot_joint_cell_calling,
    plot_joint_qc,
    plot_latent_qc_correlations,
    plot_marker_evidence,
    plot_marker_loci,
    plot_multivi_umap,
    plot_rna_qc,
    plot_training_history,
    population_palette,
    resolution_check,
    qc_reason_counts,
    retention_table,
    rna_baseline,
    save_figure,
    score_doublets,
    score_populations,
    select_atac_features,
    select_rna_features,
    small_group_masks,
    split_population_by_depth,
    set_plot_style,
    subset_paired,
    train_multivi,
    training_history_frame,
    validate_cache_manifest,
    write_cache_manifest,
    write_results,
)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the public PBMC multiome tutorial through paired MultiVI.")
    parser.add_argument("--max-epochs", type=int, default=500)
    parser.add_argument("--n-top-genes", type=int, default=3_000)
    parser.add_argument("--max-peaks", type=int, default=30_000)
    parser.add_argument("--min-cells-per-peak", type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--reuse-model",
        action="store_true",
        help="Reload a saved model whose manifest matches the current input instead of training.",
    )
    return parser.parse_args()


def run_manifest(args: argparse.Namespace, model, mdata, n_peaks_detected: int) -> dict:
    packages = ["anndata", "mudata", "muon", "scanpy", "scvi-tools", "torch", "numpy", "scipy", "scikit-learn"]
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": {package: version(package) for package in packages},
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "seed": args.seed,
        "n_genes": int(mdata["rna"].n_vars),
        "n_regions": int(mdata["atac"].n_vars),
        "n_peaks_detected_min_cells": n_peaks_detected,
        "n_nuclei": int(mdata.n_obs),
        "n_hidden": int(model.module.n_hidden),
        "n_latent": int(model.module.n_latent),
        "epochs_trained": int(len(model.history["reconstruction_loss_train"])),
        "max_epochs": args.max_epochs,
    }


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    tutorial_dir = script_dir.parent
    repo_root = tutorial_dir.parents[2]
    raw_dir = tutorial_dir / "data" / "raw"
    processed_dir = tutorial_dir / "data" / "processed"
    results_dir = tutorial_dir / "results"
    figure_dir = repo_root / "source" / "img" / "transcriptomics" / "multiome"
    model_dir = results_dir / "models" / "multivi_pbmc_paired"

    results_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    set_plot_style()

    # Joint QC on every called nucleus.
    rna, atac, barcode_metrics = load_multiome_dataset(raw_dir, processed_dir)
    score_doublets(rna, seed=args.seed)
    thresholds = choose_qc_thresholds(rna.obs)
    add_qc_flags(rna.obs, thresholds)
    atac.obs = rna.obs.copy()

    figures = [
        ("01_joint_cell_calling.png", plot_joint_cell_calling(barcode_metrics)),
        ("02_fragment_length.png", plot_fragment_length(processed_dir)),
        ("03_rna_qc.png", plot_rna_qc(rna.obs, thresholds)),
        ("04_atac_qc.png", plot_atac_qc(rna.obs, thresholds)),
        ("05_joint_qc.png", plot_joint_qc(rna.obs, thresholds)),
    ]
    for filename, fig in figures:
        save_figure(fig, figure_dir / filename)
        plt.close(fig)
    ad.AnnData(obs=rna.obs).write_h5ad(results_dir / "joint_qc_obs.h5ad")

    # Paired nuclei, features, and single-modality baselines.
    rna_paired, atac_paired = subset_paired(rna, atac)
    genes = select_rna_features(rna_paired, n_top_genes=args.n_top_genes)
    peaks = select_atac_features(
        atac_paired,
        min_cells_per_peak=args.min_cells_per_peak,
        max_peaks=None if args.max_peaks <= 0 else args.max_peaks,
    )
    normalize_rna(rna_paired)
    container = ad.AnnData(obs=rna_paired.obs.copy())
    container.obsm["X_pca_rna"] = rna_baseline(rna_paired, genes, seed=args.seed)
    container.obsm["X_lsi_atac"], lsi_depth = atac_baseline(atac_paired, peaks)

    # MultiVI on raw counts of the selected features.
    mdata = build_multivi_input(rna_paired, atac_paired, genes, peaks)
    library_check = check_observed_libraries(mdata)
    parameters = {
        "max_epochs": args.max_epochs,
        "n_top_genes": args.n_top_genes,
        "max_peaks": args.max_peaks,
        "min_cells_per_peak": args.min_cells_per_peak,
        "resolution": args.resolution,
        "seed": args.seed,
    }
    manifest_path = results_dir / "multivi_cache_manifest.json"
    history_path = results_dir / "multivi_training_history.tsv"
    if args.reuse_model:
        cached = validate_cache_manifest(mdata, manifest_path)
        if cached["parameters"] != parameters:
            raise ValueError("Saved model parameters differ from the requested run. Retrain without --reuse-model.")
        model = load_multivi(mdata, model_dir)
        history = pd.read_csv(history_path, sep="\t")
    else:
        write_cache_manifest(mdata, manifest_path, parameters=parameters)
        model = train_multivi(mdata, model_dir=model_dir, max_epochs=args.max_epochs, seed=args.seed)
        history = training_history_frame(model)
        history.to_csv(history_path, sep="\t", index=False)
    container.obsm["X_multivi"] = mdata.obsm["X_multivi"]

    # Graphs, clusters, and the latent-space audit.
    for rep_key, name in [("X_pca_rna", "rna"), ("X_lsi_atac", "atac"), ("X_multivi", "multivi")]:
        build_embedding_graph(container, rep_key, name, resolution=args.resolution, seed=args.seed)
    correlations = latent_qc_correlations(container.obsm["X_multivi"], container.obs)

    # Name joint clusters from observed RNA markers.
    rna_paired.obs["leiden_multivi"] = container.obs["leiden_multivi"].to_numpy()
    scores = score_populations(rna_paired, "leiden_multivi", seed=args.seed)
    assignment = assign_populations(scores)
    population = container.obs["leiden_multivi"].map(assignment["population"]).astype(str)
    order = [name for name in [*POPULATION_MARKERS, "Unresolved"] if name in set(population)]
    container.obs["population"] = pd.Categorical(population, categories=order)
    rna_paired.obs["population"] = container.obs["population"].values
    atac_paired.obs["population"] = container.obs["population"].values
    correlations_within = latent_qc_correlations_within(container.obsm["X_multivi"], container.obs, "population")
    palette = population_palette(order)

    gene_order = list(dict.fromkeys(gene for genes_ in POPULATION_MARKERS.values() for gene in genes_))
    evidence = marker_evidence(rna_paired, atac_paired, "population", gene_order)
    locus_specificity = marker_locus_specificity(atac_paired, rna_paired.var, atac_paired.obs["population"])
    loci = [
        locus_accessibility(atac_paired, rna_paired.var, atac_paired.obs["population"], gene)
        for gene in choose_marker_loci(locus_specificity)
    ]
    agreement = clustering_agreement(container.obs, "leiden_multivi", ["leiden_rna", "leiden_atac"])
    agreement = pd.concat(
        [agreement, clustering_agreement(container.obs, "leiden_rna", ["leiden_atac"])],
        ignore_index=True,
    )

    # Inspect what the joint clusters add beyond the broad labels.
    t_like = container.obs["population"].astype(str).isin(["T cells", "Unresolved"]).to_numpy()
    t_states = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_multivi"], T_STATE_MARKERS)
    t_states_rna = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_rna"], T_STATE_MARKERS)
    t_states_atac = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_atac"], T_STATE_MARKERS)
    dc_tables = {
        key: marker_detection(rna_paired, container.obs[key], DC_MARKERS)
        for key in ["leiden_rna", "leiden_atac", "leiden_multivi"]
    }
    dc_signal = dc_tables["leiden_rna"][["FCER1A", "CD1C", "CLEC10A"]].mean(axis=1)
    dc_like = container.obs["leiden_rna"].astype(str).eq(str(dc_signal.idxmax()))
    depth_group = split_population_by_depth(container.obs, "CD14 monocytes")
    monocyte_qc = cluster_qc_table(container.obs, sorted(container.obs.loc[depth_group.index, "leiden_multivi"].astype(str).unique()))
    cluster_group = container.obs.loc[depth_group.index, "leiden_multivi"].astype(str).groupby(depth_group).unique()
    monocyte_qc["depth_group"] = [
        next(group for group, members in cluster_group.items() if cluster in members) for cluster in monocyte_qc.index
    ]
    monocyte_groups = depth_group[~dc_like.reindex(depth_group.index).to_numpy()]
    monocyte_de = compare_rna_groups(rna_paired, monocyte_groups)
    significant = monocyte_de.loc[monocyte_de["pvals_adj"] < 0.05, "names"].unique().tolist()
    mito_link = mitochondrial_link(rna_paired, significant)
    small_groups = small_group_masks(container.obs, dc_tables["leiden_rna"], t_states_rna)
    small_group_check = resolution_check(container, small_groups, [0.8, 1.2, 1.6, 2.0, 3.0], seed=args.seed)

    figures = [
        ("06_single_modality_baselines.png", plot_baselines(container)),
        ("07_multivi_training.png", plot_training_history(history)),
        ("08_multivi_latent_qc.png", plot_latent_qc_correlations(correlations)),
        ("09_multivi_umap.png", plot_multivi_umap(container)),
        ("10_multivi_cluster_qc.png", plot_cluster_qc(container.obs)),
        ("11_marker_evidence.png", plot_marker_evidence(evidence, "population", gene_order)),
        ("12_embedding_comparison.png", plot_embedding_comparison(container, "population", palette)),
        ("13_cluster_agreement.png", plot_cluster_agreement(container.obs, "population", agreement)),
        ("14_marker_loci.png", plot_marker_loci(loci)),
    ]
    for filename, fig in figures:
        save_figure(fig, figure_dir / filename)
        plt.close(fig)

    tables = {
        "retention": retention_table(rna.obs),
        "qc_flag_counts": qc_reason_counts(rna.obs),
        "observed_library_check": library_check,
        "lsi_depth_correlations": lsi_depth,
        "multivi_latent_qc_correlations": correlations,
        "multivi_latent_qc_correlations_within_population": correlations_within,
        "cluster_marker_scores": scores.join(assignment),
        "marker_evidence": evidence,
        "marker_locus_specificity": locus_specificity,
        "t_cell_state_markers": t_states,
        "t_cell_state_markers_by_rna_cluster": t_states_rna,
        "t_cell_state_markers_by_atac_cluster": t_states_atac,
        "dc_markers_by_rna_cluster": dc_tables["leiden_rna"],
        "dc_markers_by_atac_cluster": dc_tables["leiden_atac"],
        "dc_markers_by_joint_cluster": dc_tables["leiden_multivi"],
        "cd14_monocyte_cluster_qc": monocyte_qc,
        "cd14_monocyte_depth_group_rna_comparison": monocyte_de,
        "small_group_resolution_check": small_group_check,
        "depth_group_gene_mito_link": mito_link,
        "clustering_agreement": agreement,
    }
    write_results(container, thresholds, tables, results_dir)
    pd.Series(genes, name="gene").to_csv(results_dir / "model_genes.tsv", sep="\t", index=False)
    pd.Series(peaks, name="peak").to_csv(results_dir / "model_peaks.tsv", sep="\t", index=False)
    n_peaks_detected = int((atac_paired.var["n_cells"] >= args.min_cells_per_peak).sum())
    manifest = run_manifest(args, model, mdata, n_peaks_detected)
    with (results_dir / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)
    container.write_h5ad(results_dir / "pbmc_multiome_paired_embeddings.h5ad", compression="gzip")

    print(f"Called nuclei: {rna.n_obs:,}")
    print(rna.obs["qc_category"].value_counts().to_string())
    print(f"Predicted doublets: {int(rna.obs['qc_doublet'].sum()):,}")
    print(f"Nuclei in the paired model: {mdata.n_obs:,}")
    print(f"Model features: {mdata['rna'].n_vars:,} genes, {mdata['atac'].n_vars:,} peaks")
    print(f"Epochs trained: {manifest['epochs_trained']}")
    print(assignment.to_string())
    print(agreement.to_string(index=False))
    print(t_states.round(2).to_string())
    print(monocyte_qc.round(3).to_string())
    print(f"DC-like nuclei excluded from the monocyte comparison: {int(dc_like.reindex(depth_group.index).sum())}")
    print(monocyte_groups.value_counts().to_string())
    print(monocyte_de.loc[monocyte_de["pvals_adj"] < 0.05, ["higher_in", "names", "logfoldchanges", "pvals_adj"]].to_string())
    print(small_group_check.round(3).to_string(index=False))
    print(mito_link.to_string(index=False))
    print(f"Figures: {figure_dir}")
    print(f"Results: {results_dir}")


if __name__ == "__main__":
    main()
