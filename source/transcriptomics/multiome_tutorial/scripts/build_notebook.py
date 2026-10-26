#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import nbformat as nbf


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    tutorial_dir = script_dir.parent
    output_path = tutorial_dir / "notebooks" / "multiome_multivi_tutorial.ipynb"

    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {
        "kernelspec": {
            "display_name": "Python (gencore-multiome-multivi)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
    }
    notebook["cells"] = [
        markdown(
            """
# Single-cell multiome analysis with MultiVI

This tutorial analyzes the public **10x Genomics PBMC from a healthy donor, no
cell sorting (3k)** multiome dataset, in which RNA and chromatin accessibility
were measured in the same nuclei. It asks what the two measurements show
together that each shows only partly on its own.

The workflow judges each measurement separately during QC, builds an RNA-only
and an ATAC-only baseline, trains a paired MultiVI model on raw counts, audits
the joint representation for technical structure, and names broad PBMC
populations from observed RNA markers before checking them against observed
accessibility.

## Learning objectives

By the end, you should be able to:

- identify the Cell Ranger ARC files used for secondary multiome analysis;
- calculate RNA and ATAC QC metrics separately for the same nuclei;
- summarize which nuclei have usable RNA, usable ATAC, both, or neither;
- keep raw counts for modeling and normalized values for plots;
- build RNA-only (PCA) and ATAC-only (TF-IDF/LSI) baselines;
- register a MuData object and train, save, and reload paired MultiVI;
- audit the joint representation against technical covariates;
- compare RNA-only, ATAC-only, and joint clusters for the same nuclei;
- separate observed expression and observed accessibility from model output.
"""
        ),
        markdown(
            """
## Reproducibility and compute

Create the environment with `environment.yml`, then run
`scripts/download_data.sh` and `scripts/prepare_fragment_metrics.sh`.

The fragment file is about 540 MB. The rendered notebook reloads a model
trained by `scripts/run_multivi.py` on a CPU. `OMP_NUM_THREADS` is set to 1
below and Scanpy is imported before scvi-tools. On macOS the PyTorch and
scikit-learn wheels each bundle an OpenMP runtime. In this environment, MultiVI
training stalled after Scanpy's neighbor and UMAP steps without the thread
limit, and importing scvi-tools (and with it PyTorch) before Scanpy crashed
numba's parallel QC code.
"""
        ),
        code(
            """
import os
import sys
import json
import logging
import warnings
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp/learn_gencore_multiome_xdg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/learn_gencore_multiome_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/learn_gencore_multiome_numba")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_WARNINGS", "0")
warnings.filterwarnings("ignore", message="IProgress not found.*")

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
# Import Scanpy before scvi-tools. In this environment, importing scvi-tools
# (and with it PyTorch) first made numba's parallel QC code crash the kernel.
import scanpy as sc
import scvi
from IPython.display import display

cwd = Path.cwd().resolve()
candidates = [
    cwd.parent if cwd.name == "notebooks" else cwd,
    cwd / "source" / "transcriptomics" / "multiome_tutorial",
]
tutorial_dir = next(path for path in candidates if (path / "scripts").exists())
repo_root = tutorial_dir.parents[2]
raw_dir = tutorial_dir / "data" / "raw"
processed_dir = tutorial_dir / "data" / "processed"
results_dir = tutorial_dir / "results"
figure_dir = repo_root / "source" / "img" / "transcriptomics" / "multiome"
sys.path.insert(0, str(tutorial_dir / "scripts"))

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

set_plot_style()
scvi.settings.verbosity = logging.WARNING
for package in ["scanpy", "mudata", "muon", "scvi-tools", "torch"]:
    print(f"{package}: {version(package)}")
"""
        ),
        markdown(
            """
## 1. Check the public inputs

The filtered matrix contains RNA and peak counts for called nuclei only.
`per_barcode_metrics.csv` also covers background barcodes, so it shows how
Cell Ranger ARC separated nuclei from empty droplets. The fragment file is
needed for fragment periodicity and nucleosome signal.

The dataset page is:
https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0
"""
        ),
        code(
            """
required_inputs = [
    (raw_dir, "filtered_feature_bc_matrix.h5", "raw"),
    (raw_dir, "per_barcode_metrics.csv", "raw"),
    (raw_dir, "atac_peak_annotation.tsv", "raw"),
    (raw_dir, "atac_fragments.tsv.gz", "raw"),
    (raw_dir, "atac_fragments.tsv.gz.tbi", "raw"),
    (processed_dir, "fragment_metrics.tsv", "processed"),
    (processed_dir, "fragment_length_histogram.tsv", "processed"),
]
input_status = pd.DataFrame(
    {
        "location": [location for _, _, location in required_inputs],
        "file": [name for _, name, _ in required_inputs],
        "present": [(directory / name).exists() for directory, name, _ in required_inputs],
        "size_MB": [
            round((directory / name).stat().st_size / 1e6, 1)
            if (directory / name).exists()
            else None
            for directory, name, _ in required_inputs
        ],
    }
)
display(input_status)
assert input_status["present"].all(), (
    "Run scripts/download_data.sh and scripts/prepare_fragment_metrics.sh "
    "before continuing."
)
"""
        ),
        markdown(
            """
## 2. Split the matrix into RNA and ATAC objects

Cell Ranger ARC stores genes and peaks in one matrix. The loader splits them by
`feature_types`, keeps raw counts in `layers["counts"]`, and gives both objects
the same per-nucleus metadata. The identity columns (`cell_id`, `sample_id`,
`donor_id`, `modality_status`, `has_rna`, `has_atac`) are set now so that later
analyses that mix paired and single-modality nuclei can rely on them.

In this ARC 2.0.0 output, `gex_umis_count` and `gex_genes_count` are zero for
every barcode. Exonic plus intronic UMIs equal the matrix totals, so the loader
uses those columns instead.
"""
        ),
        code(
            """
rna, atac, barcode_metrics = load_multiome_dataset(raw_dir, processed_dir)
print(rna)
print(atac)
assert rna.obs_names.equals(atac.obs_names)
"""
        ),
        markdown(
            """
## 3. Review joint cell calling and fragment periodicity

Cell Ranger ARC calls a barcode a nucleus using both measurements. Called
nuclei should separate from background barcodes in RNA UMIs and ATAC fragments
together. The fragment-length curve should show a nucleosome-free component
and periodic mono- and di-nucleosomal structure.
"""
        ),
        code(
            """
fig = plot_joint_cell_calling(barcode_metrics)
save_figure(fig, figure_dir / "01_joint_cell_calling.png")
plt.show()

fig = plot_fragment_length(processed_dir)
save_figure(fig, figure_dir / "02_fragment_length.png")
plt.show()
"""
        ),
        markdown(
            """
## 4. Flag RNA and ATAC quality separately

Each measurement gets its own review flags. RNA depth and gene-count tails use
the 2.5th percentile starting rule from the scRNA tutorial, and mitochondrial
UMIs use its 97.5th percentile. ATAC FRiP and TSS fraction use the 1st
percentile and nucleosome signal the 99th percentile, as in the scATAC
tutorial. Upper depth tails use the 99.5th percentile in both. These are
review thresholds for this dataset, not universal cutoffs.

Scrublet scores droplets from RNA. A doublet droplet contains two nuclei, so a
predicted doublet is excluded from the paired model whatever its ATAC metrics.
"""
        ),
        code(
            """
score_doublets(rna, seed=17)
thresholds = choose_qc_thresholds(rna.obs)
add_qc_flags(rna.obs, thresholds)
atac.obs = rna.obs.copy()

display(pd.Series(thresholds, name="threshold").to_frame())
display(qc_reason_counts(rna.obs))

fig = plot_rna_qc(rna.obs, thresholds)
save_figure(fig, figure_dir / "03_rna_qc.png")
plt.show()

fig = plot_atac_qc(rna.obs, thresholds)
save_figure(fig, figure_dir / "04_atac_qc.png")
plt.show()
"""
        ),
        markdown(
            """
## 5. Summarize which measurements each nucleus can contribute

Combining the flags gives four groups: both measurements usable, RNA only,
ATAC only, or neither. This tutorial's paired model keeps nuclei that pass both
and are not predicted doublets. Nuclei with one usable measurement are recorded
rather than discarded, because a mixed paired and single-modality model could
use them.
"""
        ),
        code(
            """
display(retention_table(rna.obs))

low_rna = rna.obs["qc_category"].eq("ATAC only")
display(
    rna.obs.loc[low_rna, ["total_counts", "pct_counts_mt", "atac_fragments", "frip", "tss_fraction"]]
    .median()
    .rename("median, ATAC-only nuclei")
    .to_frame()
)

fig = plot_joint_qc(rna.obs, thresholds)
save_figure(fig, figure_dir / "05_joint_qc.png")
plt.show()
"""
        ),
        markdown(
            """
## 6. Select features and build single-modality baselines

The model uses the 3,000 most variable nuclear genes (Seurat v3 rule on raw
counts) and the 30,000 peaks with the largest binary accessibility variance
among peaks detected in at least 30 nuclei. Mitochondrial genes and their
nuclear MTRNR2L copies are excluded from the modeled genes: in nuclei their
counts follow mitochondrial RNA carried over with attached cytoplasm, and vary
with how cleanly each nucleus was isolated rather than with nuclear gene
expression.

The baselines show what each measurement resolves on its own: PCA on scaled,
log-normalized RNA, and TF-IDF followed by LSI on the same peaks. LSI
components that track sequencing depth (absolute correlation above 0.5 with
log10 fragments) are dropped. Normalized values are used for these baselines
and for plots only; MultiVI receives raw counts.
"""
        ),
        code(
            """
rna_paired, atac_paired = subset_paired(rna, atac)
genes = select_rna_features(rna_paired, n_top_genes=3_000)
peaks = select_atac_features(atac_paired, min_cells_per_peak=30, max_peaks=30_000)
print(f"Paired nuclei: {rna_paired.n_obs:,}")
print(f"Peaks detected in at least 30 nuclei: {int((atac_paired.var['n_cells'] >= 30).sum()):,}")
print(f"Model genes: {len(genes):,}; model peaks: {len(peaks):,}")

normalize_rna(rna_paired)
container = ad.AnnData(obs=rna_paired.obs.copy())
container.obsm["X_pca_rna"] = rna_baseline(rna_paired, genes, seed=17)
container.obsm["X_lsi_atac"], lsi_depth = atac_baseline(atac_paired, peaks)
display(lsi_depth.head())
"""
        ),
        markdown(
            """
## 7. Build the MultiVI input and check observed libraries

MultiVI decides that a nucleus is missing a modality when its counts in that
modality sum to zero. In a paired analysis every nucleus was measured in both,
so after feature selection every nucleus must still have counts in both
modalities. Otherwise MultiVI would treat it as unmeasured in one of them.
"""
        ),
        code(
            """
mdata = build_multivi_input(rna_paired, atac_paired, genes, peaks)
print(mdata)
display(check_observed_libraries(mdata))
"""
        ),
        markdown(
            """
## 8. Train or reload MultiVI

Set `RUN_TRAINING = True` to fit the model in this notebook. The rendered copy
reloads the model generated by `scripts/run_multivi.py`.

All nuclei come from one paired sample, so no `batch_key` is passed. When
paired and single-modality nuclei are combined, the scvi-tools MultiVI tutorial
asks for the main `batch_key` to be the modality, with other batches passed as
`categorical_covariate_keys`.
"""
        ),
        code(
            """
RUN_TRAINING = False
model_dir = results_dir / "models" / "multivi_pbmc_paired"
manifest_path = results_dir / "multivi_cache_manifest.json"
history_path = results_dir / "multivi_training_history.tsv"
parameters = {
    "max_epochs": 500,
    "n_top_genes": 3_000,
    "max_peaks": 30_000,
    "min_cells_per_peak": 30,
    "resolution": 0.8,
    "seed": 17,
}

if RUN_TRAINING:
    manifest = write_cache_manifest(mdata, manifest_path, parameters=parameters)
    model = train_multivi(mdata, model_dir=model_dir, max_epochs=500, seed=17)
    history = training_history_frame(model)
    history.to_csv(history_path, sep="\\t", index=False)
else:
    required_cache = [model_dir / "model.pt", history_path, manifest_path]
    missing_cache = [str(path) for path in required_cache if not path.exists()]
    if missing_cache:
        raise FileNotFoundError(
            "Cached MultiVI outputs are missing. Run `python scripts/run_multivi.py` "
            "from the tutorial directory, or set RUN_TRAINING = True. Missing: "
            + ", ".join(missing_cache)
        )
    manifest = validate_cache_manifest(mdata, manifest_path)
    model = load_multivi(mdata, model_dir)
    history = pd.read_csv(history_path, sep="\\t")

print(model)
display(pd.Series(manifest["parameters"], name="value").to_frame())
print(f"Epochs trained: {history['reconstruction_loss_train'].notna().sum()}")
fig = plot_training_history(history)
save_figure(fig, figure_dir / "07_multivi_training.png")
plt.show()
"""
        ),
        markdown(
            """
## 9. Build graphs and clusters for all three representations

The same graph, UMAP, and Leiden settings are applied to the RNA-only PCA, the
ATAC-only LSI, and the MultiVI latent space. The three UMAPs are separate
coordinate systems; compare which nuclei group together, not their positions.
"""
        ),
        code(
            """
container.obsm["X_multivi"] = mdata.obsm["X_multivi"]
for rep_key, name in [("X_pca_rna", "rna"), ("X_lsi_atac", "atac"), ("X_multivi", "multivi")]:
    build_embedding_graph(container, rep_key, name, resolution=0.8, seed=17)
    print(f"{name}: {container.obs[f'leiden_{name}'].nunique()} clusters")

fig = plot_baselines(container)
save_figure(fig, figure_dir / "06_single_modality_baselines.png")
plt.show()
"""
        ),
        markdown(
            """
## 10. Audit the joint representation

A converged model can still organize nuclei by depth or quality. Correlate each
latent dimension with RNA and ATAC QC metrics, then inspect the UMAP and the
cluster-level QC distributions.
"""
        ),
        code(
            """
correlations = latent_qc_correlations(container.obsm["X_multivi"], container.obs)
display(correlations.reindex(correlations["correlation"].abs().sort_values(ascending=False).index).head(10))

fig = plot_latent_qc_correlations(correlations)
save_figure(fig, figure_dir / "08_multivi_latent_qc.png")
plt.show()

fig = plot_multivi_umap(container)
save_figure(fig, figure_dir / "09_multivi_umap.png")
plt.show()

fig = plot_cluster_qc(container.obs)
save_figure(fig, figure_dir / "10_multivi_cluster_qc.png")
plt.show()
"""
        ),
        markdown(
            """
## 11. Name joint clusters from observed RNA markers

Each joint cluster is scored with several markers per broad PBMC population on
observed, log-normalized RNA. A cluster takes the best-scoring population when
that score is at least 0.3 and leads the runner-up by at least 0.2; otherwise
it stays `Unresolved`. The names are provisional labels supported by marker
evidence, not validated cell identities.
"""
        ),
        code(
            """
rna_paired.obs["leiden_multivi"] = container.obs["leiden_multivi"].to_numpy()
scores = score_populations(rna_paired, "leiden_multivi", seed=17)
assignment = assign_populations(scores)
display(scores.round(2).join(assignment[["population", "margin"]].round(2)))

population = container.obs["leiden_multivi"].map(assignment["population"]).astype(str)
order = [name for name in [*POPULATION_MARKERS, "Unresolved"] if name in set(population)]
container.obs["population"] = pd.Categorical(population, categories=order)
rna_paired.obs["population"] = container.obs["population"].values
atac_paired.obs["population"] = container.obs["population"].values
palette = population_palette(order)
display(container.obs["population"].value_counts().to_frame("nuclei"))

# Overall latent-QC correlations mix cell-type differences with technical ones.
# Centering within populations asks what remains among nuclei of one type.
correlations_within = latent_qc_correlations_within(container.obsm["X_multivi"], container.obs, "population")
display(
    correlations_within.reindex(correlations_within["correlation_within"].abs().sort_values(ascending=False).index).head(8)
)
"""
        ),
        markdown(
            """
## 12. Check the labels against observed accessibility

The RNA panel shows the markers used to name clusters. The ATAC panel shows,
for the same genes, the fraction of nuclei with at least one fragment in a peak
that 10x annotates as that gene's promoter (an x marks genes with no annotated
promoter peak). Both panels are observed data. Accessibility was not used to
assign labels, although it did shape the joint clusters that were labelled.

The locus plot shows, for each population, the marker locus whose most
population-specific peak (within 20 kb of the gene's annotated TSS) differs
most from every other population. `marker_locus_specificity` reports that
value for every marker, including the ones not plotted.
"""
        ),
        code(
            """
gene_order = list(dict.fromkeys(gene for genes_ in POPULATION_MARKERS.values() for gene in genes_))
evidence = marker_evidence(rna_paired, atac_paired, "population", gene_order)
fig = plot_marker_evidence(evidence, "population", gene_order)
save_figure(fig, figure_dir / "11_marker_evidence.png")
plt.show()

locus_specificity = marker_locus_specificity(atac_paired, rna_paired.var, atac_paired.obs["population"])
display(locus_specificity.round(3))
loci = [
    locus_accessibility(atac_paired, rna_paired.var, atac_paired.obs["population"], gene)
    for gene in choose_marker_loci(locus_specificity)
]
fig = plot_marker_loci(loci)
save_figure(fig, figure_dir / "14_marker_loci.png")
plt.show()
"""
        ),
        markdown(
            """
## 13. Compare what RNA, ATAC, and the joint model resolve

Coloring all three embeddings with the same labels shows which populations
each measurement separates on its own. The contingency tables show how each
labelled population is distributed across RNA-only and ATAC-only clusters.
The adjusted Rand index (ARI) summarizes agreement between two clusterings;
it is also lowered when one clustering simply has more clusters than the
other, so compare it alongside the cluster counts.
"""
        ),
        code(
            """
agreement = clustering_agreement(container.obs, "leiden_multivi", ["leiden_rna", "leiden_atac"])
agreement = pd.concat(
    [agreement, clustering_agreement(container.obs, "leiden_rna", ["leiden_atac"])],
    ignore_index=True,
)
display(agreement.round(3))

fig = plot_embedding_comparison(container, "population", palette)
save_figure(fig, figure_dir / "12_embedding_comparison.png")
plt.show()

fig = plot_cluster_agreement(container.obs, "population", agreement)
save_figure(fig, figure_dir / "13_cluster_agreement.png")
plt.show()
"""
        ),
        markdown(
            """
## 14. Inspect what the joint clusters add

The broad labels hide structure inside each population. Three checks use
observed RNA only:

- T-cell state markers in each joint T-cell cluster, including markers commonly
  used for regulatory T cells (FOXP3, IL2RA, CTLA4, IKZF2, RTKN2, TIGIT);
- conventional DC markers (FCER1A, CD1C, CLEC10A) in the clusters of each
  representation, to see whether a small population is kept or merged;
- the CD14 monocyte joint clusters, split by whether their median ATAC fragment
  count is below the population median, with a comparison of observed RNA
  between the two groups. Nuclei in the RNA-only cluster with the strongest
  DC-marker detection are left out of that comparison.
"""
        ),
        code(
            """
t_like = container.obs["population"].astype(str).isin(["T cells", "Unresolved"]).to_numpy()
t_states = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_multivi"], T_STATE_MARKERS)
t_states_rna = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_rna"], T_STATE_MARKERS)
t_states_atac = marker_detection(rna_paired[t_like], container.obs.loc[t_like, "leiden_atac"], T_STATE_MARKERS)
display(t_states.round(2))
print("Same markers, grouped by RNA-only and by ATAC-only clusters:")
display(t_states_rna.round(2))
display(t_states_atac.round(2))

dc_tables = {
    key: marker_detection(rna_paired, container.obs[key], DC_MARKERS)
    for key in ["leiden_rna", "leiden_atac", "leiden_multivi"]
}
for key, table in dc_tables.items():
    print(key)
    display(table.sort_values("CLEC10A", ascending=False).head(3).round(2))
dc_signal = dc_tables["leiden_rna"][["FCER1A", "CD1C", "CLEC10A"]].mean(axis=1)
dc_like = container.obs["leiden_rna"].astype(str).eq(str(dc_signal.idxmax()))
display(container.obs.loc[dc_like, "leiden_multivi"].value_counts().to_frame("DC-like nuclei by joint cluster"))
"""
        ),
        code(
            """
depth_group = split_population_by_depth(container.obs, "CD14 monocytes")
monocyte_clusters = sorted(container.obs.loc[depth_group.index, "leiden_multivi"].astype(str).unique())
monocyte_qc = cluster_qc_table(container.obs, monocyte_clusters)
cluster_group = container.obs.loc[depth_group.index, "leiden_multivi"].astype(str).groupby(depth_group).unique()
monocyte_qc["depth_group"] = [
    next(group for group, members in cluster_group.items() if cluster in members) for cluster in monocyte_qc.index
]
display(monocyte_qc.round(3))

monocyte_groups = depth_group[~dc_like.reindex(depth_group.index).to_numpy()]
print(f"DC-like nuclei left out: {int(dc_like.reindex(depth_group.index).sum())}")
display(monocyte_groups.value_counts().to_frame("nuclei"))
monocyte_de = compare_rna_groups(rna_paired, monocyte_groups)
display(monocyte_de.loc[monocyte_de["pvals_adj"] < 0.05, ["higher_in", "names", "logfoldchanges", "pvals_adj", "pct_nz_group", "pct_nz_reference"]])

# Do the genes that differ track mitochondrial RNA?
significant = monocyte_de.loc[monocyte_de["pvals_adj"] < 0.05, "names"].unique().tolist()
mito_link = mitochondrial_link(rna_paired, significant)
display(mito_link)
"""
        ),
        markdown(
            """
Two small groups are clearest in the RNA-only clustering: the cluster with the
strongest conventional-DC markers and the cluster with the strongest
regulatory-T markers. The next cell asks whether the joint clustering gives
either group its own cluster as the Leiden resolution increases, and how the
ATAC-only clustering handles them. `share_of_cluster` is the fraction of the
best-matching cluster made up of the group.
"""
        ),
        code(
            """
small_groups = small_group_masks(container.obs, dc_tables["leiden_rna"], t_states_rna)
small_group_check = resolution_check(container, small_groups, [0.8, 1.2, 1.6, 2.0, 3.0], seed=17)
display(small_group_check.round(2))
"""
        ),
        markdown(
            """
## 15. Save compact, reusable results

The trained model and large objects stay outside Git. The compact results
contain thresholds, the retention table, the latent audit, cluster marker
scores, marker evidence, clustering agreement, and per-nucleus metadata with
all three UMAPs.
"""
        ),
        code(
            """
tables = {
    "retention": retention_table(rna.obs),
    "qc_flag_counts": qc_reason_counts(rna.obs),
    "observed_library_check": check_observed_libraries(mdata),
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
print("Results written to", results_dir.relative_to(repo_root))
"""
        ),
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(notebook, output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
