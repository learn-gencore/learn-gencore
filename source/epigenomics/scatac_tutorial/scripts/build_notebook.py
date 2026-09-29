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
    output_path = tutorial_dir / "notebooks" / "scatac_peakvi_tutorial.ipynb"

    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {
        "kernelspec": {
            "display_name": "Python (gencore-scatac-peakvi)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
    }
    notebook["cells"] = [
        markdown(
            """
# Single-cell ATAC-seq analysis through PEAKVI

This tutorial analyzes the public **10x Genomics 10k Human PBMCs, ATAC v2,
Chromium X** dataset. It starts with Cell Ranger ATAC outputs, performs
fragment-level and per-nucleus QC, creates a filtered peak-by-nucleus matrix,
trains PEAKVI, and audits the resulting latent space with UMAP, Leiden
clustering, technical covariates, and a cluster-marker accessibility example.

The workflow preserves raw peak counts, keeps QC and modeling separate, and
audits the learned representation for technical structure. This tutorial ends
after calling candidate cluster-associated accessible peaks with PEAKVI.

## Learning objectives

By the end, you should be able to:

- identify the 10x files needed for secondary scATAC analysis;
- interpret barcode rank, fragment periodicity, FRiP, TSS-associated fragments,
  nucleosome signal, duplicates, and mitochondrial reads;
- create transparent QC flags from the observed dataset;
- preserve raw peak counts in an AnnData layer;
- choose a peak feature set suitable for a teaching-scale PEAKVI run;
- train, save, and reload a PEAKVI model;
- construct a neighborhood graph, UMAP, and Leiden clusters from `X_peakvi`;
- check whether the learned representation still follows fragment depth or QC;
- call and interpret candidate accessible peaks for one robust cluster.
"""
        ),
        markdown(
            """
## Reproducibility and compute

Create the environment with `environment.yml`, then run
`scripts/download_data.sh` and `scripts/prepare_fragment_metrics.sh`.

The complete fragment file is about 2.3 GB. PEAKVI is much faster with a GPU.
The rendered example uses 50,000 accessible regions and a cached trained model
so that the notebook can be reviewed on a CPU. For a research analysis with a
GPU, increase the feature count only after checking memory and convergence.
"""
        ),
        code(
            """
import os
import sys
import logging
import warnings
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp/learn_gencore_scatac_xdg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/learn_gencore_scatac_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/learn_gencore_scatac_numba")
os.environ.setdefault("KMP_WARNINGS", "0")
warnings.filterwarnings("ignore", message="IProgress not found.*")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvi
from IPython.display import display

cwd = Path.cwd().resolve()
candidates = [
    cwd.parent if cwd.name == "notebooks" else cwd,
    cwd / "source" / "epigenomics" / "scatac_tutorial",
]
tutorial_dir = next(path for path in candidates if (path / "scripts").exists())
repo_root = tutorial_dir.parents[2]
raw_dir = tutorial_dir / "data" / "raw"
processed_dir = tutorial_dir / "data" / "processed"
results_dir = tutorial_dir / "results"
figure_dir = repo_root / "source" / "img" / "epigenomics" / "scatac"
sys.path.insert(0, str(tutorial_dir / "scripts"))

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
    validate_cache_manifest,
    write_cache_manifest,
    write_results,
)

set_plot_style()
scvi.settings.verbosity = logging.WARNING
print("scanpy:", version("scanpy"))
print("scvi-tools:", version("scvi-tools"))
"""
        ),
        markdown(
            """
## 1. Check the public inputs

The peak matrix contains called nuclei only. `singlecell.csv` also contains
background barcodes and Cell Ranger fragment counts, which lets us reconstruct
the barcode-rank and several QC plots. The fragment file is needed for fragment
periodicity and nucleosome signal.

The dataset page is:
https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-x-2-standard
"""
        ),
        code(
            """
required_inputs = [
    (raw_dir, "filtered_peak_bc_matrix.h5", "raw"),
    (raw_dir, "singlecell.csv", "raw"),
    (raw_dir, "peak_annotation.tsv", "raw"),
    (raw_dir, "fragments.tsv.gz", "raw"),
    (raw_dir, "fragments.tsv.gz.tbi", "raw"),
    (processed_dir, "fragment_metrics.tsv", "processed"),
    (processed_dir, "fragment_length_histogram.tsv", "processed"),
]
input_status = pd.DataFrame(
    {
        "location": [location for _, _, location in required_inputs],
        "file": [name for _, name, _ in required_inputs],
        "present": [(directory / name).exists() for directory, name, _ in required_inputs],
        "size_GB": [
            round((directory / name).stat().st_size / 1e9, 3)
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
## 2. Load the peak matrix and attach QC metadata

Raw peak counts are preserved in `adata.layers["counts"]`. The matrix has
nuclei as rows and genomic regions as columns. The peak annotation file may
contain more than one candidate gene for a region, so the loader retains the
closest entry only for descriptive metadata.
"""
        ),
        code(
            """
adata, singlecell_all = load_atac_dataset(raw_dir, processed_dir)
print(adata)
print(f"Called nuclei: {adata.n_obs:,}")
print(f"Called peaks: {adata.n_vars:,}")
"""
        ),
        markdown(
            """
## 3. Inspect cell calling and fragment periodicity

The barcode-rank curve separates called nuclei from background barcodes. The
fragment-length curve should show a strong nucleosome-free component and
periodic mono- and di-nucleosomal structure. These are sample-level checks and
should be reviewed before choosing cell-level thresholds.
"""
        ),
        code(
            """
fig = plot_barcode_rank(singlecell_all)
save_figure(fig, figure_dir / "01_barcode_rank.png")
plt.show()

fig = plot_fragment_length(processed_dir)
save_figure(fig, figure_dir / "02_fragment_length.png")
plt.show()
"""
        ),
        markdown(
            """
## 4. Define transparent QC flags

This example uses the observed lower 1% tails for fragment count, FRiP, and
TSS-associated fragment fraction, the upper 0.5% tail for fragment count, and
the upper 1% tail for nucleosome signal. These are review thresholds for this
dataset, not universal scATAC cutoffs.

`tss_fraction` is the fraction of Cell Ranger high-quality fragments assigned
near annotated TSS regions. It should not be confused with the aggregate TSS
enrichment score calculated from a TSS-centered insertion profile.
"""
        ),
        code(
            """
thresholds = choose_qc_thresholds(adata)
add_qc_flags(adata, thresholds)

display(pd.Series(thresholds, name="threshold").to_frame())
display(qc_summary(adata))

fig = plot_qc_distributions(adata, thresholds)
save_figure(fig, figure_dir / "03_qc_distributions.png")
plt.show()

fig = plot_qc_relationships(adata)
save_figure(fig, figure_dir / "04_qc_relationships.png")
plt.show()
"""
        ),
        markdown(
            """
## 5. Prepare the matrix for PEAKVI

We retain nuclei that pass all tutorial QC flags and peaks detected in at least
20 retained nuclei. The CPU-rendered notebook then keeps the 50,000 peaks with
the largest binary accessibility variance. This is a compute-aware teaching
choice. It must be reported because feature selection changes the model input.

PEAKVI models whether each region is accessible, while accounting for cell and
region detection effects. It accepts binary or count matrices. Here the
unmodified Cell Ranger counts remain in the `counts` layer.
"""
        ),
        code(
            """
peakvi_adata = filter_for_peakvi(
    adata,
    min_cells_per_peak=20,
    max_peaks=50_000,
)
print(peakvi_adata)
print(f"Retained nuclei: {peakvi_adata.n_obs:,}")
print(f"Retained peaks: {peakvi_adata.n_vars:,}")
"""
        ),
        markdown(
            """
## 6. Train or reload PEAKVI

Set `RUN_TRAINING = True` to fit the model in this notebook. The rendered copy
reloads the model generated by `scripts/run_peakvi.py`. Saving the model avoids
repeating a several-minute CPU training run and makes the result auditable.

For a multi-library dataset, pass a real library or batch column as
`batch_key` during `PEAKVI.setup_anndata`. This public example has one sample,
so there is no estimable biological batch effect.
"""
        ),
        code(
            """
RUN_TRAINING = False
model_dir = results_dir / "models" / "peakvi_pbmc"
manifest_path = results_dir / "peakvi_cache_manifest.json"
history_path = results_dir / "peakvi_training_history.tsv"

if RUN_TRAINING:
    model = train_peakvi(
        peakvi_adata,
        model_dir=model_dir,
        max_epochs=100,
        n_latent=15,
        seed=17,
    )
    history = training_history_frame(model)
    history.to_csv(history_path, sep="\t", index=False)
    manifest = write_cache_manifest(
        peakvi_adata,
        manifest_path,
        parameters={
            "max_epochs": 100,
            "max_peaks": 50_000,
            "min_cells_per_peak": 20,
            "resolution": 0.6,
            "n_latent": 15,
            "seed": 17,
        },
    )
else:
    required_cache = [model_dir / "model.pt", history_path, manifest_path]
    missing_cache = [str(path) for path in required_cache if not path.exists()]
    if missing_cache:
        raise FileNotFoundError(
            "Cached PEAKVI outputs are missing. Run `python scripts/run_peakvi.py` "
            "from the tutorial directory, or set RUN_TRAINING = True. Missing: "
            + ", ".join(missing_cache)
        )
    manifest = validate_cache_manifest(peakvi_adata, manifest_path)
    scvi.model.PEAKVI.setup_anndata(peakvi_adata, layer="counts")
    model = scvi.model.PEAKVI.load(
        model_dir,
        adata=peakvi_adata,
        accelerator="cpu",
        device="auto",
    )
    peakvi_adata.obsm["X_peakvi"] = model.get_latent_representation()
    history = pd.read_csv(history_path, sep="\t")

display(pd.Series(manifest["parameters"], name="value").to_frame())
display(history.tail())
fig = plot_training_history(history)
save_figure(fig, figure_dir / "05_peakvi_training.png")
plt.show()
"""
        ),
        markdown(
            """
## 7. Audit the PEAKVI latent representation

A model can converge while retaining unwanted technical structure. We inspect
the correlation of every latent dimension with fragment depth, FRiP,
TSS-associated fragment fraction, and nucleosome signal before interpreting
the embedding.
"""
        ),
        code(
            """
correlations = latent_qc_correlations(peakvi_adata)
display(correlations.reindex(correlations["correlation"].abs().sort_values(ascending=False).index).head(12))

fig = plot_latent_qc_correlations(correlations)
save_figure(fig, figure_dir / "06_peakvi_latent_qc.png")
plt.show()
"""
        ),
        markdown(
            """
## 8. Build the neighborhood graph, UMAP, and Leiden clusters

The graph and UMAP are computed from `X_peakvi`, not directly from the sparse
peak matrix. Coloring the same embedding by QC variables makes depth-driven or
low-quality islands visible. The cluster-level boxplots provide a second check
that a cluster is not simply a technical tail.
"""
        ),
        code(
            """
build_peakvi_graph(peakvi_adata, resolution=0.6, seed=17)
print(peakvi_adata.obs["leiden_peakvi"].value_counts().sort_index())

cluster_qc = peakvi_adata.obs.groupby("leiden_peakvi", observed=True).agg(
    n_nuclei=("leiden_peakvi", "size"),
    median_fragments=("passed_filters", "median"),
    median_frip=("frip", "median"),
    median_tss_fraction=("tss_fraction", "median"),
    median_nucleosome_signal=("nucleosome_signal", "median"),
)
display(cluster_qc)

fig = plot_peakvi_umap(peakvi_adata)
save_figure(fig, figure_dir / "07_peakvi_umap.png")
plt.show()

fig = plot_cluster_qc(peakvi_adata)
save_figure(fig, figure_dir / "08_peakvi_cluster_qc.png")
plt.show()
"""
        ),
        markdown(
            """
## 9. Call candidate cluster-associated accessible peaks

PEAKVI can compare accessibility between cell groups at individual regions.
Here cluster 0, the largest cluster after per-nucleus QC and review of the
diagnostics above, is compared with all other retained nuclei. This reference
group is a heterogeneous mixture and includes clusters 15 and 16, whose QC
profiles warrant sensitivity analysis before biological interpretation.

The model uses its default posterior change threshold (`delta=0.05`) and a
posterior expected FDR target of 0.05. In the scvi-tools implementation used
here, `delta` applies to a posterior log2 accessibility change calculated with
an estimated pseudocount. It is not an absolute five-percentage-point cutoff.

For a compact teaching table, we additionally require an estimated effect of
at least 0.10 toward the target cluster, an empirical effect of at least 0.05,
and empirical detection in at least 5% of target nuclei. These extra filters
are reporting choices, not part of PEAKVI's expected-FDR calculation.

The returned `effect_size` is a separate absolute probability difference:
`comparison - target`. A negative value therefore means greater accessibility
in the target cluster. The tested universe is the 50,000 peaks supplied to
this CPU teaching model. The rendered analysis uses 5,000 posterior samples,
which can require several gigabytes of memory at this feature count. A smaller
value such as 500 is useful for a quick preview, but thresholded calls can
change and final results should use the documented full setting.
"""
        ),
        code(
            """
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

da_summary = pd.DataFrame(
    {
        "measure": [
            "target cluster",
            "target nuclei",
            "comparison nuclei",
            "tested peaks",
            "PEAKVI posterior FDR calls",
            "candidate target-cluster markers after reporting filters",
        ],
        "value": [
            group1,
            int(peakvi_adata.obs["leiden_peakvi"].astype(str).eq(group1).sum()),
            int(peakvi_adata.obs["leiden_peakvi"].astype(str).ne(group1).sum()),
            len(da_results),
            int(da_results["is_da_fdr"].sum()),
            int(da_results["tutorial_marker"].sum()),
        ],
    }
)
display(da_summary)

marker_columns = [
    "peak",
    "gene",
    "distance",
    "effect_size",
    "emp_effect",
    "prob_da",
    "est_prob1",
    "est_prob2",
    "emp_prob1",
    "emp_prob2",
]
display(da_results.loc[da_results["tutorial_marker"], marker_columns].head(20))

fig = plot_differential_accessibility(da_results, group1=group1, group2=group2)
save_figure(fig, figure_dir / "09_peakvi_differential_accessibility.png")
plt.show()
"""
        ),
        markdown(
            """
## 10. Save compact, reusable results

The trained model is stored separately from the public input data. The compact
results contain UMAP coordinates, clusters, QC metrics, latent-to-QC
correlations, and the complete peak-level comparison.
"""
        ),
        code(
            """
write_results(peakvi_adata, thresholds, correlations, results_dir)
da_results.to_csv(
    results_dir / f"peakvi_da_cluster_{group1}_vs_rest.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
print("Results written to", results_dir.relative_to(repo_root))
"""
        ),
        markdown(
            """
## Interpretation and stopping point

At this point we have a QC-audited PEAKVI representation, unsupervised
clusters, and candidate cluster-associated accessible peaks. Cluster numbers
are not cell types, and the nearest gene is a positional annotation rather
than a demonstrated regulatory target.

In this rendered run, the two smallest clusters have lower fragment depth,
FRiP, and TSS fraction than the main groups, and one also has unusually high
nucleosome signal. They remain visible here because this is exactly the kind of
post-model QC finding that should trigger a documented sensitivity analysis
before cell-type annotation. A converged model does not make those nuclei
biological populations.

This single public donor supports a methods demonstration, descriptive
clustering, and marker accessibility within this sample. It does not support
donor-level or condition-level differential accessibility, batch-effect
estimation, or population claims. Those analyses require independent
biological replicates and a sample-aware design.

Common pitfalls at this stage include filtering with copied thresholds,
discarding rare nuclei because their QC distributions differ, treating UMAP
distance as a quantitative effect size, calling clusters cell types without
orthogonal evidence, and assuming model convergence removes all technical
structure. Cluster marker discovery is also partly circular because clustering
and peak testing use the same selected accessibility matrix.
"""
        ),
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(notebook, output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
