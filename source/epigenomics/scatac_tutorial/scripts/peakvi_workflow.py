from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
from scipy import sparse


def set_plot_style() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 220,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor="white")


def _ordered_names_sha256(values) -> str:
    payload = "\n".join(map(str, values)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_cache_manifest(
    adata: ad.AnnData,
    output_path: str | Path,
    parameters: dict | None = None,
) -> dict:
    manifest = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_names_sha256": _ordered_names_sha256(adata.obs_names),
        "var_names_sha256": _ordered_names_sha256(adata.var_names),
        "parameters": parameters or {},
    }
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def validate_cache_manifest(adata: ad.AnnData, manifest_path: str | Path) -> dict:
    manifest_path = Path(manifest_path)
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    observed = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_names_sha256": _ordered_names_sha256(adata.obs_names),
        "var_names_sha256": _ordered_names_sha256(adata.var_names),
    }
    mismatches = [key for key, value in observed.items() if manifest.get(key) != value]
    if mismatches:
        raise ValueError(
            "The cached PEAKVI model does not match the current filtered matrix: "
            + ", ".join(mismatches)
            + ". Rerun scripts/run_peakvi.py."
        )
    return manifest


def load_atac_dataset(raw_dir: str | Path, processed_dir: str | Path) -> tuple[ad.AnnData, pd.DataFrame]:
    raw_dir = Path(raw_dir)
    processed_dir = Path(processed_dir)

    adata = sc.read_10x_h5(raw_dir / "filtered_peak_bc_matrix.h5", gex_only=False)
    adata.var_names_make_unique()
    adata.layers["counts"] = sparse.csr_matrix(adata.X, dtype=np.float32)

    annotation = pd.read_csv(raw_dir / "peak_annotation.tsv", sep="\t")
    annotation["peak_id"] = (
        annotation["chrom"].astype(str)
        + ":"
        + annotation["start"].astype(str)
        + "-"
        + annotation["end"].astype(str)
    )
    annotation["absolute_distance"] = annotation["distance"].abs()
    annotation = (
        annotation.sort_values(["peak_id", "absolute_distance"])
        .drop_duplicates("peak_id", keep="first")
        .set_index("peak_id")
        .reindex(adata.var_names)
    )
    if annotation.index.has_duplicates or annotation["chrom"].isna().any():
        raise ValueError("Peak annotations could not be aligned to the peak matrix.")
    for column in ["chrom", "start", "end", "gene", "distance", "peak_type"]:
        adata.var[column] = annotation[column].to_numpy()

    singlecell_all = pd.read_csv(raw_dir / "singlecell.csv")
    singlecell = singlecell_all.set_index("barcode").reindex(adata.obs_names)
    if singlecell["is__cell_barcode"].isna().any():
        raise ValueError("Some filtered peak matrix barcodes are absent from singlecell.csv.")
    adata.obs = adata.obs.join(singlecell)

    fragment_metrics_path = processed_dir / "fragment_metrics.tsv"
    if fragment_metrics_path.exists():
        fragment_metrics = pd.read_csv(fragment_metrics_path, sep="\t", index_col="barcode")
        adata.obs = adata.obs.join(fragment_metrics)

    passed = adata.obs["passed_filters"].clip(lower=1)
    total = adata.obs["total"].clip(lower=1)
    adata.obs["frip"] = adata.obs["peak_region_fragments"] / passed
    adata.obs["tss_fraction"] = adata.obs["TSS_fragments"] / passed
    adata.obs["duplicate_fraction"] = adata.obs["duplicate"] / total
    adata.obs["mitochondrial_fraction"] = adata.obs["mitochondrial"] / total
    adata.obs["log10_fragments"] = np.log10(passed)

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    return adata, singlecell_all


def choose_qc_thresholds(adata: ad.AnnData) -> dict[str, float]:
    obs = adata.obs
    thresholds = {
        "min_fragments": float(np.floor(obs["passed_filters"].quantile(0.01) / 100) * 100),
        "max_fragments": float(np.ceil(obs["passed_filters"].quantile(0.995) / 100) * 100),
        "min_frip": float(np.floor(obs["frip"].quantile(0.01) * 100) / 100),
        "min_tss_fraction": float(np.floor(obs["tss_fraction"].quantile(0.01) * 100) / 100),
    }
    if "nucleosome_signal" in obs and obs["nucleosome_signal"].notna().any():
        thresholds["max_nucleosome_signal"] = float(np.ceil(obs["nucleosome_signal"].quantile(0.99) * 10) / 10)
    return thresholds


def add_qc_flags(adata: ad.AnnData, thresholds: dict[str, float]) -> None:
    obs = adata.obs
    obs["qc_low_fragments"] = obs["passed_filters"] < thresholds["min_fragments"]
    obs["qc_high_fragments"] = obs["passed_filters"] > thresholds["max_fragments"]
    obs["qc_low_frip"] = obs["frip"] < thresholds["min_frip"]
    obs["qc_low_tss"] = obs["tss_fraction"] < thresholds["min_tss_fraction"]
    flag_columns = ["qc_low_fragments", "qc_high_fragments", "qc_low_frip", "qc_low_tss"]
    if "max_nucleosome_signal" in thresholds:
        obs["qc_high_nucleosome"] = obs["nucleosome_signal"] > thresholds["max_nucleosome_signal"]
        flag_columns.append("qc_high_nucleosome")

    obs["qc_pass"] = ~obs[flag_columns].any(axis=1)
    readable = {
        "qc_low_fragments": "low fragments",
        "qc_high_fragments": "high fragments",
        "qc_low_frip": "low FRiP",
        "qc_low_tss": "low TSS fraction",
        "qc_high_nucleosome": "high nucleosome signal",
    }
    obs["qc_reason"] = pd.Categorical(
        [
            "pass"
            if not row.any()
            else "; ".join(readable[column] for column in flag_columns if bool(row[column]))
            for _, row in obs[flag_columns].iterrows()
        ]
    )


def filter_for_peakvi(
    adata: ad.AnnData,
    min_cells_per_peak: int = 20,
    max_peaks: int | None = 50_000,
) -> ad.AnnData:
    filtered = adata[adata.obs["qc_pass"].to_numpy()].copy()
    counts = sparse.csr_matrix(filtered.layers["counts"])
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    keep = detected >= min_cells_per_peak

    if max_peaks is not None and keep.sum() > max_peaks:
        detection_rate = detected / filtered.n_obs
        binary_variance = detection_rate * (1 - detection_rate)
        eligible = np.where(keep)[0]
        selected = eligible[np.argsort(binary_variance[eligible])[::-1][:max_peaks]]
        keep = np.zeros(filtered.n_vars, dtype=bool)
        keep[selected] = True

    filtered = filtered[:, keep].copy()
    filtered.X = sparse.csr_matrix(filtered.layers["counts"], dtype=np.float32)
    filtered.var["n_cells"] = np.asarray((filtered.X > 0).sum(axis=0)).ravel()
    filtered.var["detection_rate"] = filtered.var["n_cells"] / filtered.n_obs
    return filtered


def train_peakvi(
    adata: ad.AnnData,
    model_dir: str | Path,
    max_epochs: int = 100,
    n_latent: int = 15,
    seed: int = 17,
) -> scvi.model.PEAKVI:
    scvi.settings.seed = seed
    scvi.model.PEAKVI.setup_anndata(adata, layer="counts")
    model = scvi.model.PEAKVI(adata, n_latent=n_latent)
    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        accelerator="cpu",
        devices=1,
        batch_size=128,
    )
    model_dir = Path(model_dir)
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(model_dir, overwrite=True)
    adata.obsm["X_peakvi"] = model.get_latent_representation()
    return model


def build_peakvi_graph(adata: ad.AnnData, resolution: float = 0.6, seed: int = 17) -> None:
    sc.pp.neighbors(adata, use_rep="X_peakvi", n_neighbors=30, metric="cosine", random_state=seed)
    sc.tl.umap(adata, min_dist=0.25, random_state=seed)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="In the future, the default backend for leiden.*",
            category=FutureWarning,
        )
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added="leiden_peakvi",
            random_state=seed,
            flavor="leidenalg",
            directed=True,
            n_iterations=-1,
        )


def latent_qc_correlations(adata: ad.AnnData) -> pd.DataFrame:
    latent = np.asarray(adata.obsm["X_peakvi"])
    metrics = ["log10_fragments", "frip", "tss_fraction", "nucleosome_signal"]
    rows = []
    for component in range(latent.shape[1]):
        for metric in metrics:
            if metric not in adata.obs:
                continue
            values = adata.obs[metric].to_numpy(dtype=float)
            valid = np.isfinite(values)
            correlation = np.corrcoef(latent[valid, component], values[valid])[0, 1]
            rows.append({"component": component + 1, "metric": metric, "correlation": correlation})
    return pd.DataFrame(rows)


def plot_barcode_rank(singlecell_all: pd.DataFrame) -> plt.Figure:
    data = singlecell_all.loc[singlecell_all["barcode"] != "NO_BARCODE", ["passed_filters", "is__cell_barcode"]].copy()
    data = data.loc[data["passed_filters"] > 0].sort_values("passed_filters", ascending=False).reset_index(drop=True)
    data["rank"] = np.arange(1, len(data) + 1)
    called = data["is__cell_barcode"].eq(1)

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.scatter(data.loc[~called, "rank"], data.loc[~called, "passed_filters"], s=1, color="#b8b8b8", label="background barcode")
    ax.scatter(data.loc[called, "rank"], data.loc[called, "passed_filters"], s=2, color="#20639b", label="called nucleus")
    ax.set(xscale="log", yscale="log", xlabel="Barcode rank", ylabel="High-quality fragments", title="Cell Ranger barcode rank curve")
    ax.legend(frameon=False, markerscale=4)
    fig.tight_layout()
    return fig


def plot_fragment_length(processed_dir: str | Path) -> plt.Figure:
    histogram = pd.read_csv(Path(processed_dir) / "fragment_length_histogram.tsv", sep="\t")
    histogram = histogram.loc[histogram["fragment_length"] <= 700]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(histogram["fragment_length"], histogram["n_fragments"], color="#173f5f", linewidth=1.4)
    ax.axvline(147, color="#ed553b", linestyle="--", linewidth=1, label="147 bp")
    ax.axvline(294, color="#f6d55c", linestyle="--", linewidth=1, label="294 bp")
    ax.set(xlabel="Fragment length (bp)", ylabel="Unique fragments", title="Nucleosomal fragment periodicity")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_qc_distributions(adata: ad.AnnData, thresholds: dict[str, float]) -> plt.Figure:
    metrics = [
        ("passed_filters", "High-quality fragments", (thresholds["min_fragments"], thresholds["max_fragments"])),
        ("frip", "Fraction of fragments in peaks", (thresholds["min_frip"],)),
        ("tss_fraction", "Fraction of fragments near TSS", (thresholds["min_tss_fraction"],)),
        ("nucleosome_signal", "Nucleosome signal", (thresholds.get("max_nucleosome_signal", np.nan),)),
        ("duplicate_fraction", "Duplicate fraction", ()),
        ("mitochondrial_fraction", "Mitochondrial fraction", ()),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4))
    for ax, (column, label, lines) in zip(axes.flat, metrics):
        if column not in adata.obs:
            ax.set_visible(False)
            continue
        sns.histplot(adata.obs[column], bins=50, color="#20639b", ax=ax)
        for value in lines:
            if np.isfinite(value):
                ax.axvline(value, color="#ed553b", linestyle="--", linewidth=1.2)
        ax.set(xlabel=label, ylabel="Nuclei")
    fig.suptitle("Per-nucleus quality metrics and data-derived review thresholds", y=1.01)
    fig.tight_layout()
    return fig


def plot_qc_relationships(adata: ad.AnnData) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    color = np.where(adata.obs["qc_pass"], "#20639b", "#ed553b")
    for ax, column, label in zip(
        axes,
        ["frip", "tss_fraction", "nucleosome_signal"],
        ["FRiP", "TSS fraction", "Nucleosome signal"],
    ):
        ax.scatter(adata.obs["log10_fragments"], adata.obs[column], c=color, s=4, alpha=0.35, linewidths=0)
        ax.set(xlabel="log10 high-quality fragments", ylabel=label)
    fig.suptitle("QC metrics should be reviewed jointly")
    fig.tight_layout()
    return fig


def training_history_frame(model: scvi.model.PEAKVI) -> pd.DataFrame:
    series = {}
    for key, values in model.history.items():
        array = np.asarray(values).ravel()
        series[key] = pd.Series(array)
    return pd.DataFrame(series)


def plot_training_history(history: scvi.model.PEAKVI | pd.DataFrame) -> plt.Figure:
    if isinstance(history, pd.DataFrame):
        history_frame = history
    else:
        history_frame = training_history_frame(history)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    plotted = False
    loss_keys = [
        ("reconstruction_loss_train", "training", "#20639b"),
        ("reconstruction_loss_validation", "validation", "#ed553b"),
    ]
    if not any(key in history_frame for key, _, _ in loss_keys):
        loss_keys = [
            ("elbo_train", "training ELBO", "#20639b"),
            ("elbo_validation", "validation ELBO", "#ed553b"),
        ]
    for key, label, color in loss_keys:
        if key in history_frame:
            values = history_frame[key].dropna()
            ax.plot(values.index + 1, values.values, label=label, color=color)
            plotted = True
    if not plotted:
        ax.text(0.5, 0.5, "Training history keys were unavailable", ha="center", va="center")
    ax.set(xlabel="Epoch", ylabel="Loss", title="PEAKVI convergence")
    if plotted:
        ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_latent_qc_correlations(correlations: pd.DataFrame) -> plt.Figure:
    pivot = correlations.pivot(index="metric", columns="component", values="correlation")
    fig, ax = plt.subplots(figsize=(10, 3.6))
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={"label": "Pearson correlation"})
    ax.set(title="Technical association of PEAKVI latent dimensions", xlabel="PEAKVI latent dimension", ylabel="QC metric")
    fig.tight_layout()
    return fig


def _umap_panel(ax: plt.Axes, coordinates: np.ndarray, values, title: str, categorical: bool = False) -> None:
    if categorical:
        categories = pd.Categorical(values)
        palette = sns.color_palette("tab20", n_colors=max(len(categories.categories), 1))
        for category, color in zip(categories.categories, palette):
            mask = categories == category
            ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=4, alpha=0.75, linewidths=0, color=color, label=str(category))
        ax.legend(frameon=False, markerscale=2.5, fontsize=7)
    else:
        points = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=np.asarray(values), cmap="viridis", s=4, alpha=0.75, linewidths=0)
        plt.colorbar(points, ax=ax, fraction=0.045, pad=0.02)
    ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_peakvi_umap(adata: ad.AnnData) -> plt.Figure:
    coordinates = adata.obsm["X_umap"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    _umap_panel(axes[0, 0], coordinates, adata.obs["leiden_peakvi"], "Leiden clusters", categorical=True)
    _umap_panel(axes[0, 1], coordinates, adata.obs["log10_fragments"], "Fragment depth")
    _umap_panel(axes[1, 0], coordinates, adata.obs["frip"], "FRiP")
    _umap_panel(axes[1, 1], coordinates, adata.obs["tss_fraction"], "TSS fraction")
    fig.suptitle("PEAKVI latent-space diagnostics")
    fig.tight_layout()
    return fig


def plot_cluster_qc(adata: ad.AnnData) -> plt.Figure:
    columns = ["log10_fragments", "frip", "tss_fraction", "nucleosome_signal"]
    labels = ["log10 fragments", "FRiP", "TSS fraction", "Nucleosome signal"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    data = adata.obs.reset_index(drop=True)
    for ax, column, label in zip(axes.flat, columns, labels):
        sns.boxplot(data=data, x="leiden_peakvi", y=column, color="#5b9bd5", showfliers=False, ax=ax)
        ax.set(xlabel="Leiden cluster", ylabel=label)
    fig.suptitle("QC metrics by PEAKVI cluster")
    fig.tight_layout()
    return fig


def run_peakvi_differential_accessibility(
    model: scvi.model.PEAKVI,
    adata: ad.AnnData,
    group1: str,
    group2: str | None = None,
    delta: float = 0.05,
    fdr_target: float = 0.05,
    n_samples_overall: int = 5_000,
    min_abs_estimated_effect: float = 0.10,
    min_abs_empirical_effect: float = 0.05,
    min_empirical_detection: float = 0.05,
    seed: int = 17,
) -> pd.DataFrame:
    group1 = str(group1)
    group2 = None if group2 is None else str(group2)
    if group2 is not None and group1 == group2:
        raise ValueError("Differential accessibility requires two different groups.")
    observed_groups = adata.obs["leiden_peakvi"].astype(str)
    if group1 not in set(observed_groups) or (
        group2 is not None and group2 not in set(observed_groups)
    ):
        raise ValueError("Comparison groups must exist in adata.obs['leiden_peakvi'].")

    scvi.settings.seed = seed
    results = model.differential_accessibility(
        adata=adata,
        groupby="leiden_peakvi",
        group1=[group1],
        group2=group2,
        mode="change",
        delta=delta,
        batch_size=128,
        fdr_target=fdr_target,
        batch_correction=False,
        all_stats=True,
        use_permutation=False,
        silent=True,
        n_samples_overall=n_samples_overall,
    )

    annotation_columns = [
        column
        for column in ["chrom", "start", "end", "gene", "distance", "peak_type", "n_cells", "detection_rate"]
        if column in adata.var
    ]
    results = results.join(adata.var[annotation_columns], how="left")
    results.index.name = "peak"
    results.insert(0, "peak", results.index.astype(str))
    results["mean_estimated_accessibility"] = (results["est_prob1"] + results["est_prob2"]) / 2
    results["abs_effect_size"] = results["effect_size"].abs()
    group2_label = f"cluster {group2}" if group2 is not None else "all other clusters"
    results["more_accessible_in"] = np.where(
        results["effect_size"] < 0,
        f"cluster {group1}",
        group2_label,
    )
    results["tutorial_marker"] = (
        results["is_da_fdr"].astype(bool)
        & results["effect_size"].le(-min_abs_estimated_effect)
        & results["emp_effect"].le(-min_abs_empirical_effect)
        & results["emp_prob1"].ge(min_empirical_detection)
    )
    results["comparison"] = f"cluster {group1} vs {group2_label}"
    results["delta"] = delta
    results["fdr_target"] = fdr_target
    results["n_samples_overall"] = n_samples_overall
    results["seed"] = seed
    results["min_abs_estimated_effect"] = min_abs_estimated_effect
    results["min_abs_empirical_effect"] = min_abs_empirical_effect
    results["min_empirical_detection"] = min_empirical_detection
    return results.sort_values(
        ["is_da_fdr", "abs_effect_size", "prob_da"],
        ascending=[False, False, False],
    )


def plot_differential_accessibility(
    results: pd.DataFrame,
    group1: str,
    group2: str | None = None,
    n_top_per_group: int = 6,
) -> plt.Figure:
    called = results["is_da_fdr"].astype(bool)
    marker = results["tutorial_marker"].astype(bool)
    group2_label = f"cluster {group2}" if group2 is not None else "all other clusters"
    group1_color = "#20639b"
    group2_color = "#ed553b"

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
    axes[0].scatter(
        results.loc[~called, "mean_estimated_accessibility"],
        results.loc[~called, "effect_size"],
        s=3,
        alpha=0.12,
        color="#8f8f8f",
        linewidths=0,
        rasterized=True,
    )
    for mask, color, label in [
        (called & results["effect_size"].lt(0), group1_color, f"cluster {group1}"),
        (called & results["effect_size"].gt(0), group2_color, group2_label),
    ]:
        axes[0].scatter(
            results.loc[mask, "mean_estimated_accessibility"],
            results.loc[mask, "effect_size"],
            s=4,
            alpha=0.35,
            color=color,
            label=f"more accessible in {label}",
            linewidths=0,
            rasterized=True,
        )
    axes[0].axhline(0, color="#444444", linewidth=0.8)
    axes[0].set(
        xlabel="Mean PEAKVI estimated accessibility",
        ylabel=f"Estimated effect: {group2_label} minus cluster {group1}",
        title=(
            f"{int(called.sum()):,} posterior FDR calls; "
            f"{int(marker.sum()):,} candidates after reporting filters"
        ),
    )
    axes[0].legend(frameon=False, markerscale=3)

    group1_top = results.loc[called & results["effect_size"].lt(0)].nsmallest(
        n_top_per_group, "effect_size"
    )
    group2_top = results.loc[called & results["effect_size"].gt(0)].nlargest(
        n_top_per_group, "effect_size"
    )
    top = pd.concat([group1_top, group2_top]).sort_values("effect_size")
    labels = []
    for _, row in top.iterrows():
        gene = str(row.get("gene", "")).strip()
        labels.append(f"{gene} | {row['peak']}" if gene and gene.lower() != "nan" else str(row["peak"]))
    bar_colors = np.where(top["effect_size"] < 0, group1_color, group2_color)
    y = np.arange(len(top))
    axes[1].barh(y, top["effect_size"], color=bar_colors)
    axes[1].set_yticks(y, labels=labels, fontsize=8)
    axes[1].axvline(0, color="#444444", linewidth=0.8)
    axes[1].set(
        xlabel=f"Estimated effect: {group2_label} minus cluster {group1}",
        title="Largest called accessibility differences",
    )
    fig.suptitle(f"PEAKVI differential accessibility: cluster {group1} vs {group2_label}", y=1.01)
    fig.tight_layout()
    return fig


def qc_summary(adata: ad.AnnData) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "measure": [
                "called nuclei",
                "nuclei passing tutorial QC",
                "nuclei flagged",
                "median high-quality fragments",
                "median FRiP",
                "median TSS fraction",
            ],
            "value": [
                adata.n_obs,
                int(adata.obs["qc_pass"].sum()),
                int((~adata.obs["qc_pass"]).sum()),
                float(adata.obs["passed_filters"].median()),
                float(adata.obs["frip"].median()),
                float(adata.obs["tss_fraction"].median()),
            ],
        }
    )


def write_results(
    adata: ad.AnnData,
    thresholds: dict[str, float],
    correlations: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "qc_thresholds.json").open("w") as handle:
        json.dump(thresholds, handle, indent=2)
    correlations.to_csv(output_dir / "peakvi_latent_qc_correlations.tsv", sep="\t", index=False)
    metadata = adata.obs[
        ["leiden_peakvi", "passed_filters", "frip", "tss_fraction", "nucleosome_signal"]
    ].copy()
    metadata["UMAP1"] = adata.obsm["X_umap"][:, 0]
    metadata["UMAP2"] = adata.obsm["X_umap"][:, 1]
    metadata.to_csv(output_dir / "cell_metadata.tsv.gz", sep="\t", compression="gzip")
