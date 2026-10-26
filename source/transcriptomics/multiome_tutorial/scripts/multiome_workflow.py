from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc  # keep before scvi: see README, OpenMP on macOS
import scvi
import seaborn as sns
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

SAMPLE_ID = "pbmc_unsorted_3k"
DONOR_ID = "pbmc_unsorted_3k_donor"

QC_CATEGORIES = ["pass both", "RNA only", "ATAC only", "neither"]

RNA_FLAGS = {
    "qc_rna_low_umis": "low RNA UMIs",
    "qc_rna_high_umis": "high RNA UMIs",
    "qc_rna_low_genes": "few detected genes",
    "qc_rna_high_mt": "high mitochondrial fraction",
}
ATAC_FLAGS = {
    "qc_atac_low_fragments": "low ATAC fragments",
    "qc_atac_high_fragments": "high ATAC fragments",
    "qc_atac_low_frip": "low FRiP",
    "qc_atac_low_tss": "low TSS fraction",
    "qc_atac_high_nucleosome": "high nucleosome signal",
}

# Broad PBMC populations and the RNA markers used to name clusters. Each list
# uses several genes so that no label rests on a single transcript.
POPULATION_MARKERS = {
    "T cells": ["CD3D", "CD3E", "CD3G", "IL7R"],
    "NK cells": ["GNLY", "NKG7", "KLRD1", "KLRF1"],
    "B cells": ["MS4A1", "CD79A", "BANK1", "PAX5"],
    "CD14 monocytes": ["CD14", "LYZ", "VCAN", "FCN1"],
    "FCGR3A monocytes": ["FCGR3A", "MS4A7", "CDKN1C", "LST1"],
    "Conventional DCs": ["FCER1A", "CD1C", "CLEC10A"],
    "Plasmacytoid DCs": ["LILRA4", "CLEC4C", "IL3RA"],
}

# Markers for states within the broad populations, used to inspect joint clusters
# after naming. The regulatory T-cell genes follow commonly used PBMC references.
T_STATE_MARKERS = [
    "CD4", "CD8A", "CD8B", "CCR7", "SELL", "LEF1", "IL7R",
    "GZMK", "CCL5", "NKG7", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "RTKN2", "TIGIT",
]
DC_MARKERS = ["FCER1A", "CD1C", "CLEC10A", "CD14", "LYZ"]

LATENT_QC_METRICS = [
    "log10_umis",
    "pct_counts_mt",
    "intronic_fraction",
    "doublet_score",
    "log10_fragments",
    "frip",
    "tss_fraction",
    "nucleosome_signal",
]


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


def _mudata_fingerprint(mdata: md.MuData) -> dict:
    return {
        "n_obs": int(mdata.n_obs),
        "n_genes": int(mdata["rna"].n_vars),
        "n_regions": int(mdata["atac"].n_vars),
        "obs_names_sha256": _ordered_names_sha256(mdata.obs_names),
        "gene_names_sha256": _ordered_names_sha256(mdata["rna"].var_names),
        "region_names_sha256": _ordered_names_sha256(mdata["atac"].var_names),
    }


def write_cache_manifest(
    mdata: md.MuData,
    output_path: str | Path,
    parameters: dict | None = None,
) -> dict:
    manifest = {**_mudata_fingerprint(mdata), "parameters": parameters or {}}
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def validate_cache_manifest(mdata: md.MuData, manifest_path: str | Path) -> dict:
    with Path(manifest_path).open() as handle:
        manifest = json.load(handle)
    observed = _mudata_fingerprint(mdata)
    mismatches = [key for key, value in observed.items() if manifest.get(key) != value]
    if mismatches:
        raise ValueError(
            "The cached MultiVI model does not match the current model input: "
            + ", ".join(mismatches)
            + ". Rerun scripts/run_multivi.py."
        )
    return manifest


def _parse_region(names: pd.Index) -> pd.DataFrame:
    parts = pd.Series(names, index=names).str.rsplit(":", n=1, expand=True)
    bounds = parts[1].str.split("-", expand=True)
    return pd.DataFrame(
        {
            "chrom": parts[0].to_numpy(),
            "start": bounds[0].astype(int).to_numpy(),
            "end": bounds[1].astype(int).to_numpy(),
        },
        index=names,
    )


def _annotate_peaks(atac: ad.AnnData, annotation_path: Path) -> None:
    annotation = pd.read_csv(annotation_path, sep="\t")
    annotation["peak_id"] = (
        annotation["chrom"].astype(str)
        + ":"
        + annotation["start"].astype(str)
        + "-"
        + annotation["end"].astype(str)
    )
    missing = set(atac.var_names) - set(annotation["peak_id"])
    if missing:
        raise ValueError(f"{len(missing)} peaks are absent from the peak annotation file.")

    annotation["absolute_distance"] = annotation["distance"].abs()
    nearest = (
        annotation.sort_values(["peak_id", "absolute_distance"])
        .drop_duplicates("peak_id", keep="first")
        .set_index("peak_id")
        .reindex(atac.var_names)
    )
    promoter_genes = (
        annotation.loc[annotation["peak_type"].eq("promoter")]
        .groupby("peak_id")["gene"]
        .agg(lambda genes: ";".join(sorted(set(genes.astype(str)))))
        .reindex(atac.var_names)
        .fillna("")
    )
    coordinates = _parse_region(atac.var_names)
    for column in ["chrom", "start", "end"]:
        atac.var[column] = coordinates[column].to_numpy()
    atac.var["nearest_gene"] = nearest["gene"].astype(str).to_numpy()
    atac.var["distance"] = nearest["distance"].to_numpy()
    atac.var["peak_type"] = nearest["peak_type"].astype(str).to_numpy()
    atac.var["promoter_genes"] = promoter_genes.to_numpy()


def load_multiome_dataset(
    raw_dir: str | Path,
    processed_dir: str | Path,
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    """Split the Cell Ranger ARC matrix into RNA and ATAC objects with shared QC metadata."""
    raw_dir = Path(raw_dir)
    processed_dir = Path(processed_dir)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique.*")
        combined = sc.read_10x_h5(raw_dir / "filtered_feature_bc_matrix.h5", gex_only=False)
    feature_types = combined.var["feature_types"].astype(str)
    rna = combined[:, feature_types.eq("Gene Expression").to_numpy()].copy()
    atac = combined[:, feature_types.eq("Peaks").to_numpy()].copy()
    rna.var_names_make_unique()
    for modality in (rna, atac):
        modality.X = sparse.csr_matrix(modality.X, dtype=np.float32)
        modality.layers["counts"] = modality.X.copy()

    barcode_metrics = pd.read_csv(raw_dir / "per_barcode_metrics.csv")
    called = barcode_metrics.set_index("barcode").reindex(rna.obs_names)
    if called["is_cell"].isna().any() or called["is_cell"].ne(1).any():
        raise ValueError("Every matrix barcode should be a called cell in per_barcode_metrics.csv.")

    obs = pd.DataFrame(index=rna.obs_names)
    obs["cell_id"] = SAMPLE_ID + "_" + obs.index.astype(str)
    obs["original_nucleus_id"] = obs["cell_id"]
    obs["sample_id"] = SAMPLE_ID
    obs["donor_id"] = DONOR_ID
    obs["modality_status"] = "paired"
    obs["has_rna"] = True
    obs["has_atac"] = True
    obs["atac_barcode"] = called["atac_barcode"].astype(str)
    # In this Cell Ranger ARC 2.0.0 output, gex_umis_count and gex_genes_count are
    # zero for every barcode. Exonic plus intronic UMIs equal the matrix totals.
    for column in [
        "gex_exonic_umis",
        "gex_intronic_umis",
        "atac_raw_reads",
        "atac_dup_reads",
        "atac_mitochondrial_reads",
        "atac_fragments",
        "atac_TSS_fragments",
        "atac_peak_region_fragments",
    ]:
        obs[column] = called[column].to_numpy()

    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    for column in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        obs[column] = rna.obs[column].to_numpy()
    obs["log10_umis"] = np.log10(obs["total_counts"].clip(lower=1))
    barcode_umis = (obs["gex_exonic_umis"] + obs["gex_intronic_umis"]).clip(lower=1)
    obs["intronic_fraction"] = obs["gex_intronic_umis"] / barcode_umis

    fragments = obs["atac_fragments"].clip(lower=1)
    raw_reads = obs["atac_raw_reads"].clip(lower=1)
    obs["log10_fragments"] = np.log10(fragments)
    obs["frip"] = obs["atac_peak_region_fragments"] / fragments
    obs["tss_fraction"] = obs["atac_TSS_fragments"] / fragments
    obs["duplicate_fraction"] = obs["atac_dup_reads"] / raw_reads
    obs["atac_mitochondrial_fraction"] = obs["atac_mitochondrial_reads"] / raw_reads

    fragment_metrics_path = processed_dir / "fragment_metrics.tsv"
    if fragment_metrics_path.exists():
        fragment_metrics = pd.read_csv(fragment_metrics_path, sep="\t", index_col="barcode")
        obs = obs.join(fragment_metrics[["nucleosome_signal"]])

    _annotate_peaks(atac, raw_dir / "atac_peak_annotation.tsv")

    rna.obs = obs.copy()
    atac.obs = obs.copy()
    return rna, atac, barcode_metrics


def score_doublets(rna: ad.AnnData, seed: int = 17) -> None:
    """Score droplets for RNA doublets. A doublet affects both modalities of that droplet."""
    counts = ad.AnnData(X=rna.layers["counts"].copy(), obs=pd.DataFrame(index=rna.obs_names))
    sc.pp.scrublet(counts, random_state=seed, verbose=False)
    rna.obs["doublet_score"] = counts.obs["doublet_score"].to_numpy()
    rna.obs["predicted_doublet"] = counts.obs["predicted_doublet"].astype(bool).to_numpy()


QC_QUANTILES = {
    "min_umis": 0.025,
    "max_umis": 0.995,
    "min_genes": 0.025,
    "max_pct_mt": 0.975,
    "min_fragments": 0.025,
    "max_fragments": 0.995,
    "min_frip": 0.01,
    "min_tss_fraction": 0.01,
    "max_nucleosome_signal": 0.99,
}
QC_COLUMNS = {
    "min_umis": "total_counts",
    "max_umis": "total_counts",
    "min_genes": "n_genes_by_counts",
    "max_pct_mt": "pct_counts_mt",
    "min_fragments": "atac_fragments",
    "max_fragments": "atac_fragments",
    "min_frip": "frip",
    "min_tss_fraction": "tss_fraction",
    "max_nucleosome_signal": "nucleosome_signal",
}


def choose_qc_thresholds(obs: pd.DataFrame, quantiles: dict[str, float] | None = None) -> dict[str, float]:
    """Dataset-specific review thresholds from the observed tails of each metric."""
    quantiles = {**QC_QUANTILES, **(quantiles or {})}
    steps = {"total_counts": 10, "n_genes_by_counts": 10, "atac_fragments": 10, "pct_counts_mt": 0.1}
    thresholds = {}
    for name, quantile in quantiles.items():
        column = QC_COLUMNS[name]
        if column not in obs or obs[column].isna().all():
            continue
        value = obs[column].quantile(quantile)
        step = steps.get(column, 0.01)
        rounder = np.floor if name.startswith("min_") else np.ceil
        thresholds[name] = float(np.round(rounder(value / step) * step, 4))
    return thresholds


def _reasons(flags: pd.DataFrame, readable: dict[str, str]) -> list[str]:
    return [
        "pass" if not row.any() else "; ".join(readable[column] for column in flags.columns if bool(row[column]))
        for _, row in flags.iterrows()
    ]


def add_qc_flags(obs: pd.DataFrame, thresholds: dict[str, float]) -> None:
    """Record modality-specific QC flags without removing any nucleus."""
    obs["qc_rna_low_umis"] = obs["total_counts"] < thresholds["min_umis"]
    obs["qc_rna_high_umis"] = obs["total_counts"] > thresholds["max_umis"]
    obs["qc_rna_low_genes"] = obs["n_genes_by_counts"] < thresholds["min_genes"]
    obs["qc_rna_high_mt"] = obs["pct_counts_mt"] > thresholds["max_pct_mt"]
    obs["qc_atac_low_fragments"] = obs["atac_fragments"] < thresholds["min_fragments"]
    obs["qc_atac_high_fragments"] = obs["atac_fragments"] > thresholds["max_fragments"]
    obs["qc_atac_low_frip"] = obs["frip"] < thresholds["min_frip"]
    obs["qc_atac_low_tss"] = obs["tss_fraction"] < thresholds["min_tss_fraction"]
    atac_flags = dict(ATAC_FLAGS)
    if "max_nucleosome_signal" in thresholds:
        obs["qc_atac_high_nucleosome"] = obs["nucleosome_signal"] > thresholds["max_nucleosome_signal"]
    else:
        atac_flags.pop("qc_atac_high_nucleosome")

    obs["qc_pass_rna"] = ~obs[list(RNA_FLAGS)].any(axis=1)
    obs["qc_pass_atac"] = ~obs[list(atac_flags)].any(axis=1)
    obs["qc_reason_rna"] = pd.Categorical(_reasons(obs[list(RNA_FLAGS)], RNA_FLAGS))
    obs["qc_reason_atac"] = pd.Categorical(_reasons(obs[list(atac_flags)], atac_flags))
    category = np.select(
        [
            obs["qc_pass_rna"] & obs["qc_pass_atac"],
            obs["qc_pass_rna"] & ~obs["qc_pass_atac"],
            ~obs["qc_pass_rna"] & obs["qc_pass_atac"],
        ],
        QC_CATEGORIES[:3],
        default=QC_CATEGORIES[3],
    )
    obs["qc_category"] = pd.Categorical(category, categories=QC_CATEGORIES)
    obs["qc_doublet"] = obs["predicted_doublet"].astype(bool)
    obs["qc_pass_joint"] = obs["qc_pass_rna"] & obs["qc_pass_atac"] & ~obs["qc_doublet"]


def retention_table(obs: pd.DataFrame) -> pd.DataFrame:
    table = (
        obs.groupby(["qc_category", "qc_doublet"], observed=False)
        .size()
        .rename("n_nuclei")
        .reset_index()
    )
    table["retained_for_paired_model"] = table["qc_category"].eq("pass both") & ~table["qc_doublet"]
    return table


def qc_reason_counts(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for flags in (RNA_FLAGS, ATAC_FLAGS):
        for column, label in flags.items():
            if column in obs:
                rows.append({"flag": label, "n_nuclei": int(obs[column].sum())})
    rows.append({"flag": "predicted doublet", "n_nuclei": int(obs["qc_doublet"].sum())})
    return pd.DataFrame(rows)


def subset_paired(rna: ad.AnnData, atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    keep = rna.obs["qc_pass_joint"].to_numpy()
    rna_paired = rna[keep].copy()
    atac_paired = atac[keep].copy()
    if not rna_paired.obs_names.equals(atac_paired.obs_names):
        raise ValueError("RNA and ATAC nuclei are not in the same order.")
    return rna_paired, atac_paired


def select_rna_features(rna: ad.AnnData, n_top_genes: int = 3_000) -> pd.Index:
    """Most variable nuclear genes.

    Mitochondrial genes and their nuclear MTRNR2L copies are excluded. In nuclei
    their counts follow mitochondrial RNA carried over with attached cytoplasm,
    not expression by the nucleus being profiled.
    """
    rna.var["mito_derived"] = rna.var_names.str.startswith(("MT-", "MTRNR2L"))
    candidates = rna[:, ~rna.var["mito_derived"].to_numpy()].copy()
    sc.pp.highly_variable_genes(candidates, flavor="seurat_v3", n_top_genes=n_top_genes, layer="counts")
    genes = candidates.var_names[candidates.var["highly_variable"].to_numpy()]
    rna.var["highly_variable"] = rna.var_names.isin(genes)
    return genes


def select_atac_features(
    atac: ad.AnnData,
    min_cells_per_peak: int = 30,
    max_peaks: int | None = 30_000,
) -> pd.Index:
    counts = sparse.csr_matrix(atac.layers["counts"])
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    atac.var["n_cells"] = detected
    atac.var["detection_rate"] = detected / atac.n_obs
    keep = detected >= min_cells_per_peak
    if max_peaks is not None and keep.sum() > max_peaks:
        binary_variance = atac.var["detection_rate"].to_numpy() * (1 - atac.var["detection_rate"].to_numpy())
        eligible = np.where(keep)[0]
        selected = eligible[np.argsort(binary_variance[eligible])[::-1][:max_peaks]]
        keep = np.zeros(atac.n_vars, dtype=bool)
        keep[selected] = True
    atac.var["selected_for_model"] = keep
    return atac.var_names[keep]


def normalize_rna(rna: ad.AnnData) -> None:
    """Library-size normalization and log transform for plots and the RNA baseline only."""
    rna.X = rna.layers["counts"].copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)


def rna_baseline(rna: ad.AnnData, genes: pd.Index, n_pcs: int = 30, seed: int = 17) -> np.ndarray:
    subset = rna[:, genes].copy()
    sc.pp.scale(subset, max_value=10)
    sc.pp.pca(subset, n_comps=n_pcs, random_state=seed)
    return subset.obsm["X_pca"]


def atac_baseline(
    atac: ad.AnnData,
    peaks: pd.Index,
    n_comps: int = 30,
    max_depth_correlation: float = 0.5,
) -> tuple[np.ndarray, pd.DataFrame]:
    """TF-IDF and LSI on the modeled peaks, dropping components that mainly track depth."""
    import muon as mu

    subset = atac[:, peaks].copy()
    subset.X = sparse.csr_matrix(subset.layers["counts"], dtype=np.float32)
    mu.atac.pp.tfidf(subset, scale_factor=1e4)
    mu.atac.tl.lsi(subset, n_comps=n_comps)
    lsi = np.asarray(subset.obsm["X_lsi"])
    depth = subset.obs["log10_fragments"].to_numpy(dtype=float)
    correlations = pd.DataFrame(
        {
            "component": np.arange(1, lsi.shape[1] + 1),
            "depth_correlation": [np.corrcoef(lsi[:, index], depth)[0, 1] for index in range(lsi.shape[1])],
        }
    )
    correlations["used"] = correlations["depth_correlation"].abs() <= max_depth_correlation
    return lsi[:, correlations["used"].to_numpy()], correlations


def build_multivi_input(
    rna: ad.AnnData,
    atac: ad.AnnData,
    genes: pd.Index,
    peaks: pd.Index,
) -> md.MuData:
    """Raw counts for the selected features. Normalized matrices never enter the model."""
    modalities = {}
    for name, source, features in [("rna", rna, genes), ("atac", atac, peaks)]:
        counts = sparse.csr_matrix(source[:, features].layers["counts"], dtype=np.float32)
        modalities[name] = ad.AnnData(
            X=counts,
            obs=pd.DataFrame(index=source.obs_names),
            var=source.var.loc[features].copy(),
        )
    mdata = md.MuData(modalities)
    mdata.obs = rna.obs.copy().reindex(mdata.obs_names)
    return mdata


def check_observed_libraries(mdata: md.MuData) -> pd.DataFrame:
    """MultiVI treats a cell with zero counts in a modality as unmeasured in that modality."""
    rows = []
    for name in ["rna", "atac"]:
        library = np.asarray(mdata[name].X.sum(axis=1)).ravel()
        rows.append(
            {
                "modality": name,
                "features": mdata[name].n_vars,
                "nuclei": mdata[name].n_obs,
                "nuclei_with_zero_counts": int((library == 0).sum()),
                "min_counts": float(library.min()),
                "median_counts": float(np.median(library)),
            }
        )
    table = pd.DataFrame(rows)
    if table["nuclei_with_zero_counts"].gt(0).any():
        raise ValueError(
            "Some paired nuclei have no counts in a modality after feature selection. "
            "MultiVI would treat that modality as unmeasured for them.\n" + table.to_string(index=False)
        )
    return table


def setup_multivi(mdata: md.MuData) -> None:
    scvi.model.MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna", "atac_layer": "atac"},
    )


def train_multivi(
    mdata: md.MuData,
    model_dir: str | Path,
    max_epochs: int = 500,
    seed: int = 17,
) -> scvi.model.MULTIVI:
    scvi.settings.seed = seed
    setup_multivi(mdata)
    model = scvi.model.MULTIVI(
        mdata,
        n_genes=mdata["rna"].n_vars,
        n_regions=mdata["atac"].n_vars,
    )
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
    mdata.obsm["X_multivi"] = model.get_latent_representation()
    return model


def load_multivi(mdata: md.MuData, model_dir: str | Path) -> scvi.model.MULTIVI:
    setup_multivi(mdata)
    model = scvi.model.MULTIVI.load(model_dir, adata=mdata, accelerator="cpu", device="auto")
    mdata.obsm["X_multivi"] = model.get_latent_representation()
    return model


def training_history_frame(model: scvi.model.MULTIVI) -> pd.DataFrame:
    series = {}
    for key, values in model.history.items():
        series[key] = pd.Series(np.asarray(values).ravel())
    return pd.DataFrame(series)


def build_embedding_graph(
    container: ad.AnnData,
    rep_key: str,
    name: str,
    resolution: float = 0.8,
    n_neighbors: int = 20,
    seed: int = 17,
) -> None:
    """Neighbor graph, UMAP, and Leiden clusters from one representation, stored under `name`."""
    sc.pp.neighbors(
        container,
        use_rep=rep_key,
        n_neighbors=n_neighbors,
        metric="cosine",
        key_added=name,
        random_state=seed,
    )
    sc.tl.umap(container, neighbors_key=name, min_dist=0.25, random_state=seed)
    container.obsm[f"X_umap_{name}"] = container.obsm.pop("X_umap")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="In the future, the default backend for leiden.*",
            category=FutureWarning,
        )
        sc.tl.leiden(
            container,
            resolution=resolution,
            neighbors_key=name,
            key_added=f"leiden_{name}",
            random_state=seed,
            flavor="leidenalg",
            directed=True,
            n_iterations=-1,
        )


def latent_qc_correlations(latent: np.ndarray, obs: pd.DataFrame) -> pd.DataFrame:
    latent = np.asarray(latent)
    rows = []
    for component in range(latent.shape[1]):
        for metric in LATENT_QC_METRICS:
            if metric not in obs:
                continue
            values = obs[metric].to_numpy(dtype=float)
            valid = np.isfinite(values)
            correlation = np.corrcoef(latent[valid, component], values[valid])[0, 1]
            rows.append({"component": component + 1, "metric": metric, "correlation": correlation})
    return pd.DataFrame(rows)


def latent_qc_correlations_within(latent: np.ndarray, obs: pd.DataFrame, groupby: str) -> pd.DataFrame:
    """Latent-QC correlations after removing each group's mean from both variables.

    Cell types differ in depth and quality, so an overall correlation can reflect
    cell type. Centering within groups asks whether a dimension tracks QC among
    nuclei of the same population.
    """
    latent = pd.DataFrame(np.asarray(latent), index=obs.index)
    groups = obs[groupby].astype(str)
    centered_latent = latent - latent.groupby(groups).transform("mean")
    rows = []
    for metric in LATENT_QC_METRICS:
        if metric not in obs:
            continue
        values = obs[metric].astype(float)
        centered = values - values.groupby(groups).transform("mean")
        valid = centered.notna().to_numpy()
        for component in latent.columns:
            correlation = np.corrcoef(centered_latent.loc[valid, component], centered[valid])[0, 1]
            rows.append({"component": component + 1, "metric": metric, "correlation_within": correlation})
    return pd.DataFrame(rows)


def score_populations(rna: ad.AnnData, groupby: str, seed: int = 17) -> pd.DataFrame:
    """Mean marker score per cluster from observed, log-normalized RNA."""
    scored = rna.copy()
    columns = []
    for population, genes in POPULATION_MARKERS.items():
        present = [gene for gene in genes if gene in scored.var_names]
        column = f"score_{population}"
        sc.tl.score_genes(scored, present, score_name=column, random_state=seed)
        columns.append(column)
    scores = scored.obs.groupby(groupby, observed=True)[columns].mean()
    scores.columns = list(POPULATION_MARKERS)
    return scores


def assign_populations(scores: pd.DataFrame, min_score: float = 0.3, min_margin: float = 0.2) -> pd.DataFrame:
    """Name each cluster by its highest marker score, or leave it unresolved."""
    ordered = np.sort(scores.to_numpy(), axis=1)
    table = pd.DataFrame(index=scores.index)
    table["best_population"] = scores.idxmax(axis=1)
    table["best_score"] = ordered[:, -1]
    table["margin"] = ordered[:, -1] - ordered[:, -2]
    table["population"] = np.where(
        (table["best_score"] >= min_score) & (table["margin"] >= min_margin),
        table["best_population"],
        "Unresolved",
    )
    return table


def marker_detection(rna: ad.AnnData, groups: pd.Series, genes: list[str]) -> pd.DataFrame:
    """Fraction of nuclei in each group with at least one UMI for each gene."""
    present = [gene for gene in genes if gene in rna.var_names]
    detected = pd.DataFrame(
        sparse.csr_matrix(rna[:, present].layers["counts"]).toarray() > 0,
        index=rna.obs_names,
        columns=present,
    )
    groups = pd.Series(np.asarray(groups), index=rna.obs_names)
    table = detected.groupby(groups, observed=True).mean()
    table.insert(0, "n_nuclei", groups.value_counts().reindex(table.index).to_numpy())
    return table


def cluster_qc_table(obs: pd.DataFrame, clusters: list[str], cluster_key: str = "leiden_multivi") -> pd.DataFrame:
    columns = ["total_counts", "pct_counts_mt", "atac_fragments", "frip", "tss_fraction", "nucleosome_signal", "doublet_score"]
    subset = obs.loc[obs[cluster_key].astype(str).isin(clusters)]
    table = subset.groupby(subset[cluster_key].astype(str), observed=True)[columns].median()
    table.insert(0, "n_nuclei", subset[cluster_key].astype(str).value_counts().reindex(table.index).to_numpy())
    return table.sort_values("atac_fragments")


def split_population_by_depth(
    obs: pd.DataFrame,
    population: str,
    cluster_key: str = "leiden_multivi",
    label_key: str = "population",
) -> pd.Series:
    """Label a population's joint clusters as lower or higher ATAC depth.

    A cluster is 'lower ATAC depth' when its median fragment count is below the
    median of all nuclei in the population.
    """
    members = obs.loc[obs[label_key].astype(str).eq(population)]
    population_median = members["atac_fragments"].median()
    cluster_medians = members.groupby(members[cluster_key].astype(str), observed=True)["atac_fragments"].median()
    lower = set(cluster_medians.index[cluster_medians < population_median])
    return pd.Series(
        np.where(members[cluster_key].astype(str).isin(lower), "lower ATAC depth", "higher ATAC depth"),
        index=members.index,
        name="depth_group",
    )


def compare_rna_groups(rna: ad.AnnData, groups: pd.Series) -> pd.DataFrame:
    """Wilcoxon comparison of observed, log-normalized RNA between two groups of nuclei.

    Mitochondrial genes are left out so the comparison asks about nuclear
    expression. MTRNR2L genes are kept and reported, as a check on whether the
    groups still differ in mitochondrial-derived RNA.
    """
    subset = rna[groups.index, ~rna.var["mt"].to_numpy()].copy()
    subset.obs["comparison_group"] = pd.Categorical(groups.to_numpy())
    sc.tl.rank_genes_groups(subset, "comparison_group", method="wilcoxon", pts=True)
    results = []
    for group in subset.obs["comparison_group"].cat.categories:
        table = sc.get.rank_genes_groups_df(subset, group)
        table.insert(0, "higher_in", group)
        results.append(table)
    return pd.concat(results, ignore_index=True)


def mitochondrial_link(rna: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    """Spearman correlation between each gene's share of a nucleus's UMIs and its mitochondrial percentage."""
    from scipy.stats import spearmanr

    counts = sparse.csr_matrix(rna.layers["counts"])
    totals = np.asarray(counts.sum(axis=1)).ravel()
    rows = []
    for gene in genes:
        share = np.asarray(counts[:, rna.var_names.get_loc(gene)].todense()).ravel() / totals
        rho, p_value = spearmanr(share, rna.obs["pct_counts_mt"].to_numpy())
        rows.append({"gene": gene, "spearman_rho_with_pct_mt": rho, "p_value": p_value, "n_nuclei": rna.n_obs})
    return pd.DataFrame(rows)


def best_cluster_share(mask: pd.Series, clusters: pd.Series) -> dict:
    """How well one clustering isolates a group of nuclei.

    Returns the cluster in which the group makes up the largest fraction, that
    fraction, and how many of the group's nuclei it holds.
    """
    clusters = clusters.astype(str)
    share = mask.groupby(clusters).mean()
    best = share.idxmax()
    return {
        "best_cluster": best,
        "share_of_cluster": float(share[best]),
        "group_nuclei_in_cluster": int(mask[clusters.eq(best)].sum()),
        "group_size": int(mask.sum()),
    }


def small_group_masks(obs: pd.DataFrame, dc_by_rna: pd.DataFrame, t_states_by_rna: pd.DataFrame, min_nuclei: int = 10) -> dict[str, pd.Series]:
    """RNA-only clusters with the strongest conventional-DC and regulatory-T marker detection."""
    dc_cluster = dc_by_rna.loc[dc_by_rna["n_nuclei"] >= min_nuclei, ["FCER1A", "CD1C", "CLEC10A"]].mean(axis=1).idxmax()
    treg = t_states_by_rna.loc[t_states_by_rna["n_nuclei"] >= min_nuclei, ["FOXP3", "IL2RA", "IKZF2", "RTKN2"]]
    treg_cluster = treg.mean(axis=1).idxmax()
    rna_clusters = obs["leiden_rna"].astype(str)
    return {
        "DC-like (RNA-only cluster)": rna_clusters.eq(str(dc_cluster)),
        "Treg-like (RNA-only cluster)": rna_clusters.eq(str(treg_cluster)),
    }


def resolution_check(
    container: ad.AnnData,
    groups: dict[str, pd.Series],
    resolutions: list[float],
    seed: int = 17,
) -> pd.DataFrame:
    """Whether each small group gets its own joint cluster as Leiden resolution increases."""
    scratch = ad.AnnData(obs=container.obs[[]].copy())
    scratch.obsp["multivi_connectivities"] = container.obsp["multivi_connectivities"]
    scratch.obsp["multivi_distances"] = container.obsp["multivi_distances"]
    scratch.uns["multivi"] = container.uns["multivi"]
    rows = []
    for name, mask in groups.items():
        rows.append({"group": name, "clustering": "ATAC-only", "resolution": np.nan, "n_clusters": container.obs["leiden_atac"].nunique(), **best_cluster_share(mask, container.obs["leiden_atac"])})
    for resolution in resolutions:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            sc.tl.leiden(
                scratch,
                resolution=resolution,
                neighbors_key="multivi",
                key_added="clusters",
                random_state=seed,
                flavor="leidenalg",
                directed=True,
                n_iterations=-1,
            )
        for name, mask in groups.items():
            rows.append(
                {
                    "group": name,
                    "clustering": "joint",
                    "resolution": resolution,
                    "n_clusters": scratch.obs["clusters"].nunique(),
                    **best_cluster_share(mask, scratch.obs["clusters"]),
                }
            )
    return pd.DataFrame(rows)


def cluster_contingency(obs: pd.DataFrame, row_key: str, col_key: str) -> pd.DataFrame:
    table = pd.crosstab(obs[row_key], obs[col_key])
    return table.div(table.sum(axis=1), axis=0)


def clustering_agreement(obs: pd.DataFrame, reference_key: str, other_keys: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "comparison": [f"{reference_key} vs {key}" for key in other_keys],
            "adjusted_rand_index": [adjusted_rand_score(obs[reference_key], obs[key]) for key in other_keys],
        }
    )


def promoter_accessibility(atac: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    """Per-nucleus indicator of any fragment in a gene's annotated promoter peaks."""
    counts = sparse.csc_matrix(atac.layers["counts"])
    promoter_lists = atac.var["promoter_genes"].str.split(";")
    columns = {}
    for gene in genes:
        mask = promoter_lists.apply(lambda values, gene=gene: gene in values).to_numpy()
        if mask.any():
            columns[gene] = np.asarray(counts[:, mask].sum(axis=1)).ravel() > 0
        else:
            columns[gene] = np.full(atac.n_obs, np.nan)
    return pd.DataFrame(columns, index=atac.obs_names)


def marker_evidence(rna: ad.AnnData, atac: ad.AnnData, groupby: str, genes: list[str]) -> pd.DataFrame:
    """Observed RNA detection and observed promoter accessibility for each group and gene."""
    present = [gene for gene in genes if gene in rna.var_names]
    counts = pd.DataFrame(
        sparse.csr_matrix(rna[:, present].layers["counts"]).toarray() > 0,
        index=rna.obs_names,
        columns=present,
    )
    expression = pd.DataFrame(
        sparse.csr_matrix(rna[:, present].X).toarray(),
        index=rna.obs_names,
        columns=present,
    )
    promoters = promoter_accessibility(atac, present)
    groups = rna.obs[groupby]
    rows = []
    for group in groups.cat.categories if hasattr(groups, "cat") else sorted(groups.unique()):
        mask = (groups == group).to_numpy()
        for gene in present:
            rows.append(
                {
                    groupby: group,
                    "gene": gene,
                    "n_nuclei": int(mask.sum()),
                    "rna_fraction_detected": float(counts.loc[mask, gene].mean()),
                    "rna_mean_log_expression": float(expression.loc[mask, gene].mean()),
                    "promoter_fraction_accessible": float(promoters.loc[mask, gene].mean()),
                }
            )
    return pd.DataFrame(rows)


def locus_accessibility(
    atac: ad.AnnData,
    rna_var: pd.DataFrame,
    groups: pd.Series,
    gene: str,
    flank: int = 20_000,
) -> tuple[pd.DataFrame, dict]:
    """Fraction of nuclei per group with any fragment in each peak near one gene.

    The 10x gene `interval` spans the gene's annotated transcription start sites,
    not its full body, so the window is centered on the TSS region.
    """
    chrom, bounds = str(rna_var.loc[gene, "interval"]).rsplit(":", 1)
    gene_start, gene_end = (int(value) for value in bounds.split("-"))
    window = (gene_start - flank, gene_end + flank)
    in_window = (
        atac.var["chrom"].eq(chrom)
        & atac.var["end"].ge(window[0])
        & atac.var["start"].le(window[1])
    ).to_numpy()
    counts = sparse.csc_matrix(atac.layers["counts"])[:, in_window]
    detected = counts > 0
    peaks = atac.var.loc[in_window, ["start", "end"]].copy()
    groups = pd.Series(pd.Categorical(groups), index=atac.obs_names)
    fractions = {}
    for group in groups.cat.categories:
        mask = (groups == group).to_numpy()
        fractions[group] = np.asarray(detected[mask].mean(axis=0)).ravel()
    table = pd.DataFrame(fractions, index=peaks.index).T
    locus = {"gene": gene, "chrom": chrom, "gene_start": gene_start, "gene_end": gene_end, "window": window, "peaks": peaks}
    return table, locus


def marker_locus_specificity(
    atac: ad.AnnData,
    rna_var: pd.DataFrame,
    groups: pd.Series,
    flank: int = 20_000,
) -> pd.DataFrame:
    """For each population's marker genes, the most population-specific peak near its TSS.

    Specificity is the fraction of the population's nuclei with a fragment in the
    peak minus the highest fraction in any other population.
    """
    rows = []
    for population, genes in POPULATION_MARKERS.items():
        if population not in set(groups.astype(str)):
            continue
        for gene in genes:
            if gene not in rna_var.index:
                continue
            table, _ = locus_accessibility(atac, rna_var, groups, gene, flank=flank)
            if table.shape[1] == 0:
                rows.append({"population": population, "gene": gene, "n_peaks": 0})
                continue
            others = table.drop(index=population).max(axis=0)
            specificity = table.loc[population] - others
            best = specificity.idxmax()
            rows.append(
                {
                    "population": population,
                    "gene": gene,
                    "n_peaks": table.shape[1],
                    "best_peak": best,
                    "fraction_in_population": float(table.loc[population, best]),
                    "max_fraction_elsewhere": float(others[best]),
                    "specificity": float(specificity[best]),
                }
            )
    return pd.DataFrame(rows)


def choose_marker_loci(specificity: pd.DataFrame) -> list[str]:
    """One marker locus per population: the one with the most specific peak."""
    ranked = specificity.dropna(subset=["specificity"]).sort_values("specificity", ascending=False)
    return ranked.drop_duplicates("population")["gene"].tolist()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def population_palette(populations: list[str]) -> dict[str, tuple]:
    colors = sns.color_palette("tab10", n_colors=10)
    palette = {}
    for index, population in enumerate(p for p in populations if p != "Unresolved"):
        palette[population] = colors[index % len(colors)]
    palette["Unresolved"] = (0.6, 0.6, 0.6)
    return palette


QC_COLORS = {
    "pass both": "#20639b",
    "RNA only": "#3caea3",
    "ATAC only": "#f6d55c",
    "neither": "#ed553b",
}


def plot_joint_cell_calling(barcode_metrics: pd.DataFrame) -> plt.Figure:
    data = barcode_metrics[["atac_fragments", "is_cell"]].copy()
    data["gex_umis"] = barcode_metrics["gex_exonic_umis"] + barcode_metrics["gex_intronic_umis"]
    data = data.loc[(data["gex_umis"] > 0) & (data["atac_fragments"] > 0)]
    called = data["is_cell"].eq(1)
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    ax.scatter(
        data.loc[~called, "atac_fragments"],
        data.loc[~called, "gex_umis"],
        s=1,
        color="#b8b8b8",
        alpha=0.3,
        linewidths=0,
        rasterized=True,
        label="background barcode",
    )
    ax.scatter(
        data.loc[called, "atac_fragments"],
        data.loc[called, "gex_umis"],
        s=3,
        color="#20639b",
        alpha=0.6,
        linewidths=0,
        rasterized=True,
        label="called nucleus",
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="ATAC fragments per barcode",
        ylabel="RNA UMIs per barcode",
        title="Cell Ranger ARC calls nuclei from both measurements",
    )
    ax.legend(frameon=False, markerscale=5, loc="upper left")
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


def _qc_histograms(obs: pd.DataFrame, metrics: list[tuple], title: str, log_columns: set[str]) -> plt.Figure:
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 3.8))
    for ax, (column, label, lines) in zip(np.atleast_1d(axes), metrics):
        if column not in obs:
            ax.set_visible(False)
            continue
        sns.histplot(obs[column], bins=50, color="#20639b", ax=ax, log_scale=column in log_columns)
        for value in lines:
            if value is not None and np.isfinite(value):
                ax.axvline(value, color="#ed553b", linestyle="--", linewidth=1.2)
        ax.set(xlabel=label, ylabel="Nuclei")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


def plot_rna_qc(obs: pd.DataFrame, thresholds: dict[str, float]) -> plt.Figure:
    metrics = [
        ("total_counts", "RNA UMIs", (thresholds["min_umis"], thresholds["max_umis"])),
        ("n_genes_by_counts", "Detected genes", (thresholds["min_genes"],)),
        ("pct_counts_mt", "Mitochondrial UMIs (%)", (thresholds["max_pct_mt"],)),
        ("intronic_fraction", "Intronic UMI fraction", ()),
        ("doublet_score", "Scrublet doublet score", ()),
    ]
    return _qc_histograms(obs, metrics, "RNA quality metrics and data-derived review thresholds", {"total_counts", "n_genes_by_counts"})


def plot_atac_qc(obs: pd.DataFrame, thresholds: dict[str, float]) -> plt.Figure:
    metrics = [
        ("atac_fragments", "ATAC fragments", (thresholds["min_fragments"], thresholds["max_fragments"])),
        ("frip", "Fraction of fragments in peaks", (thresholds["min_frip"],)),
        ("tss_fraction", "Fraction of fragments near TSS", (thresholds["min_tss_fraction"],)),
        ("nucleosome_signal", "Nucleosome signal", (thresholds.get("max_nucleosome_signal"),)),
    ]
    return _qc_histograms(obs, metrics, "ATAC quality metrics and data-derived review thresholds", {"atac_fragments"})


def plot_joint_qc(obs: pd.DataFrame, thresholds: dict[str, float]) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0), gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    for category in QC_CATEGORIES:
        mask = obs["qc_category"].eq(category) & ~obs["qc_doublet"]
        axes[0].scatter(
            obs.loc[mask, "atac_fragments"],
            obs.loc[mask, "total_counts"],
            s=5,
            alpha=0.7,
            linewidths=0,
            color=QC_COLORS[category],
            label=f"{category} ({int(mask.sum()):,})",
        )
    doublets = obs["qc_doublet"]
    axes[0].scatter(
        obs.loc[doublets, "atac_fragments"],
        obs.loc[doublets, "total_counts"],
        s=14,
        facecolors="none",
        edgecolors="black",
        linewidths=0.6,
        label=f"predicted doublet ({int(doublets.sum()):,})",
    )
    for value in (thresholds["min_fragments"], thresholds["max_fragments"]):
        axes[0].axvline(value, color="#444444", linestyle="--", linewidth=0.8)
    for value in (thresholds["min_umis"], thresholds["max_umis"]):
        axes[0].axhline(value, color="#444444", linestyle="--", linewidth=0.8)
    axes[0].set(xscale="log", yscale="log", xlabel="ATAC fragments", ylabel="RNA UMIs", title="Depth in each measurement")
    axes[0].legend(frameon=False, fontsize=8, markerscale=1.5)

    colors = obs["qc_category"].map(QC_COLORS).astype(str)
    axes[1].scatter(obs["pct_counts_mt"], obs["tss_fraction"], s=5, alpha=0.7, linewidths=0, c=colors)
    axes[1].axvline(thresholds["max_pct_mt"], color="#444444", linestyle="--", linewidth=0.8)
    axes[1].axhline(thresholds["min_tss_fraction"], color="#444444", linestyle="--", linewidth=0.8)
    axes[1].set(xlabel="Mitochondrial UMIs (%)", ylabel="TSS fraction", title="RNA quality against ATAC quality")

    counts = obs["qc_category"].value_counts().reindex(QC_CATEGORIES)
    axes[2].barh(QC_CATEGORIES[::-1], counts.reindex(QC_CATEGORIES[::-1]).to_numpy(), color=[QC_COLORS[c] for c in QC_CATEGORIES[::-1]])
    for y, value in enumerate(counts.reindex(QC_CATEGORIES[::-1]).to_numpy()):
        axes[2].text(value, y, f" {int(value):,}", va="center", fontsize=9)
    axes[2].set(xlabel="Nuclei", title="Which measurements pass QC")
    fig.suptitle("Joint QC: each nucleus is judged separately in each modality", y=1.02)
    fig.tight_layout()
    return fig


def plot_training_history(history: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    plotted = False
    for key, label, color in [
        ("reconstruction_loss_train", "training", "#20639b"),
        ("reconstruction_loss_validation", "validation", "#ed553b"),
    ]:
        if key in history:
            values = history[key].dropna()
            ax.plot(values.index + 1, values.values, label=label, color=color)
            plotted = True
    if not plotted:
        ax.text(0.5, 0.5, "Training history keys were unavailable", ha="center", va="center")
    ax.set(xlabel="Epoch", ylabel="Reconstruction loss", title="MultiVI convergence")
    if plotted:
        ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_latent_qc_correlations(correlations: pd.DataFrame) -> plt.Figure:
    pivot = correlations.pivot(index="metric", columns="component", values="correlation")
    pivot = pivot.reindex([metric for metric in LATENT_QC_METRICS if metric in pivot.index])
    fig, ax = plt.subplots(figsize=(11, 4.4))
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={"label": "Pearson correlation"})
    ax.set(title="Technical association of MultiVI latent dimensions", xlabel="MultiVI latent dimension", ylabel="QC metric")
    fig.tight_layout()
    return fig


def _umap_panel(ax: plt.Axes, coordinates: np.ndarray, values, title: str, categorical: bool = False, palette=None) -> None:
    if categorical:
        categories = pd.Categorical(values)
        if palette is None:
            colors = sns.color_palette("tab20", n_colors=max(len(categories.categories), 1))
            palette = dict(zip(categories.categories, colors))
        for category in categories.categories:
            mask = np.asarray(categories == category)
            ax.scatter(
                coordinates[mask, 0],
                coordinates[mask, 1],
                s=5,
                alpha=0.8,
                linewidths=0,
                color=palette[category],
                label=str(category),
            )
        ax.legend(frameon=False, markerscale=2.5, fontsize=7, loc="best")
    else:
        points = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=np.asarray(values), cmap="viridis", s=5, alpha=0.8, linewidths=0)
        plt.colorbar(points, ax=ax, fraction=0.045, pad=0.02)
    ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_baselines(container: ad.AnnData) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    _umap_panel(axes[0, 0], container.obsm["X_umap_rna"], container.obs["leiden_rna"], "RNA only: Leiden clusters", categorical=True)
    _umap_panel(axes[0, 1], container.obsm["X_umap_rna"], container.obs["log10_umis"], "RNA only: log10 UMIs")
    _umap_panel(axes[1, 0], container.obsm["X_umap_atac"], container.obs["leiden_atac"], "ATAC only: Leiden clusters", categorical=True)
    _umap_panel(axes[1, 1], container.obsm["X_umap_atac"], container.obs["log10_fragments"], "ATAC only: log10 fragments")
    fig.suptitle("Single-modality baselines on the same nuclei")
    fig.tight_layout()
    return fig


def plot_multivi_umap(container: ad.AnnData) -> plt.Figure:
    coordinates = container.obsm["X_umap_multivi"]
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    _umap_panel(axes[0, 0], coordinates, container.obs["leiden_multivi"], "Leiden clusters", categorical=True)
    _umap_panel(axes[0, 1], coordinates, container.obs["log10_umis"], "log10 RNA UMIs")
    _umap_panel(axes[0, 2], coordinates, container.obs["pct_counts_mt"], "Mitochondrial UMIs (%)")
    _umap_panel(axes[1, 0], coordinates, container.obs["log10_fragments"], "log10 ATAC fragments")
    _umap_panel(axes[1, 1], coordinates, container.obs["frip"], "FRiP")
    _umap_panel(axes[1, 2], coordinates, container.obs["tss_fraction"], "TSS fraction")
    fig.suptitle("MultiVI joint embedding: clusters and technical gradients")
    fig.tight_layout()
    return fig


def plot_cluster_qc(obs: pd.DataFrame, groupby: str = "leiden_multivi") -> plt.Figure:
    columns = ["log10_umis", "pct_counts_mt", "doublet_score", "log10_fragments", "frip", "tss_fraction"]
    labels = ["log10 RNA UMIs", "Mitochondrial UMIs (%)", "Doublet score", "log10 ATAC fragments", "FRiP", "TSS fraction"]
    fig, axes = plt.subplots(2, 3, figsize=(17, 8.5))
    data = obs.reset_index(drop=True)
    for ax, column, label in zip(axes.flat, columns, labels):
        sns.boxplot(data=data, x=groupby, y=column, color="#5b9bd5", showfliers=False, ax=ax)
        ax.set(xlabel="MultiVI Leiden cluster", ylabel=label)
    fig.suptitle("QC metrics by MultiVI cluster")
    fig.tight_layout()
    return fig


def plot_marker_evidence(evidence: pd.DataFrame, groupby: str, gene_order: list[str]) -> plt.Figure:
    groups = list(pd.unique(evidence[groupby]))
    fig, axes = plt.subplots(1, 2, figsize=(16, 0.42 * len(groups) + 2.6), sharey=True)
    for ax, value, title, cmap in [
        (axes[0], "rna_fraction_detected", "Observed RNA: fraction of nuclei with UMIs", "Blues"),
        (axes[1], "promoter_fraction_accessible", "Observed ATAC: fraction with promoter fragments", "Oranges"),
    ]:
        data = evidence.loc[evidence["gene"].isin(gene_order)].copy()
        vmax = float(np.ceil(np.nanmax(data[value].to_numpy()) * 10) / 10)
        title = f"{title}\n(color and size scaled to 0-{vmax:g})"
        data["x"] = data["gene"].map({gene: index for index, gene in enumerate(gene_order)})
        data["y"] = data[groupby].map({group: index for index, group in enumerate(groups)})
        points = ax.scatter(
            data["x"],
            data["y"],
            s=20 + 260 * data[value].fillna(0) / vmax,
            c=data[value],
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            edgecolors="#555555",
            linewidths=0.3,
        )
        missing = data.loc[data[value].isna()]
        ax.scatter(missing["x"], missing["y"], marker="x", s=18, color="#999999", linewidths=0.8)
        ax.set_xticks(range(len(gene_order)), labels=gene_order, rotation=90, fontsize=8)
        ax.set_yticks(range(len(groups)), labels=groups, fontsize=9)
        ax.set_xlim(-0.7, len(gene_order) - 0.3)
        ax.set_ylim(len(groups) - 0.5, -0.5)
        ax.grid(False)
        ax.set_title(title, fontsize=11)
        plt.colorbar(points, ax=ax, fraction=0.03, pad=0.01)
    fig.tight_layout()
    return fig


def plot_embedding_comparison(container: ad.AnnData, label_key: str, palette: dict) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8))
    for ax, name, title in [
        (axes[0], "rna", "RNA only (PCA)"),
        (axes[1], "atac", "ATAC only (LSI)"),
        (axes[2], "multivi", "Joint (MultiVI)"),
    ]:
        _umap_panel(ax, container.obsm[f"X_umap_{name}"], container.obs[label_key], title, categorical=True, palette=palette)
    for ax in axes[:2]:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
    fig.suptitle("The same nuclei and labels in three separately computed embeddings")
    fig.tight_layout()
    return fig


def plot_cluster_agreement(obs: pd.DataFrame, label_key: str, agreement: pd.DataFrame) -> plt.Figure:
    ari = dict(zip(agreement["comparison"], agreement["adjusted_rand_index"]))
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.6))
    for ax, key, title in [
        (axes[0], "leiden_rna", "RNA-only clusters"),
        (axes[1], "leiden_atac", "ATAC-only clusters"),
    ]:
        table = cluster_contingency(obs, label_key, key)
        sns.heatmap(
            table,
            cmap="Blues",
            vmin=0,
            vmax=1,
            annot=table.map(lambda value: f"{value:.2f}" if value >= 0.05 else ""),
            fmt="",
            annot_kws={"fontsize": 7},
            ax=ax,
            cbar_kws={"label": "Fraction of population"},
        )
        ax.set(xlabel=title, ylabel="")
        ax.set_title(f"{title} (n = {obs[key].nunique()})", fontsize=11)
    axes[0].set_ylabel("Population (joint clusters named by RNA markers)")
    fig.suptitle(
        "How each population is distributed across single-modality clusters"
        f"   |   ARI, RNA-only vs ATAC-only clusters: {ari.get('leiden_rna vs leiden_atac', np.nan):.2f}",
        y=1.02,
    )
    fig.tight_layout()
    return fig


def plot_marker_loci(locus_tables: list[tuple[pd.DataFrame, dict]], vmax: float = 0.6) -> plt.Figure:
    n_loci = len(locus_tables)
    n_groups = len(locus_tables[0][0].index)
    fig, axes = plt.subplots(
        n_loci,
        1,
        figsize=(10.5, (0.24 * n_groups + 1.0) * n_loci),
        layout="constrained",
    )
    axes = np.atleast_1d(axes)
    cmap = plt.get_cmap("Oranges")
    for ax, (table, locus) in zip(axes, locus_tables):
        groups = list(table.index)
        origin = locus["gene_start"]
        starts = (locus["peaks"]["start"].to_numpy() - origin) / 1_000
        widths = np.maximum((locus["peaks"]["end"] - locus["peaks"]["start"]).to_numpy() / 1_000, 0.6)
        for row, group in enumerate(groups):
            for start, width, value in zip(starts, widths, table.loc[group].to_numpy()):
                ax.add_patch(
                    plt.Rectangle((start, row - 0.42), width, 0.84, color=cmap(min(value / vmax, 1.0)), linewidth=0)
                )
        tss_row = len(groups) + 0.1
        tss_end = max((locus["gene_end"] - origin) / 1_000, 0.3)
        ax.plot([0, tss_end], [tss_row, tss_row], color="#173f5f", linewidth=6, solid_capstyle="butt")
        ax.text(tss_end + 0.5, tss_row, f"{locus['gene']} TSS", fontsize=9, color="#173f5f", va="center", fontweight="bold")
        ax.set_xlim((locus["window"][0] - origin) / 1_000, (locus["window"][1] - origin) / 1_000)
        ax.set_ylim(tss_row + 0.8, -0.6)
        ax.set_yticks(range(len(groups)), labels=groups, fontsize=7)
        ax.grid(False)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_xlabel(f"kb from the first annotated {locus['gene']} TSS ({locus['chrom']}:{origin:,})", fontsize=8)
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=plt.Normalize(0, vmax))
    fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.01, label="Fraction of nuclei with a fragment in the peak")
    fig.suptitle("Observed accessibility at peaks within 20 kb of marker-gene TSSs")
    return fig


def write_results(
    container: ad.AnnData,
    thresholds: dict[str, float],
    tables: dict[str, pd.DataFrame],
    output_dir: str | Path,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "qc_thresholds.json").open("w") as handle:
        json.dump(thresholds, handle, indent=2)
    for name, table in tables.items():
        table.to_csv(output_dir / f"{name}.tsv", sep="\t", index=True)
    metadata = container.obs.copy()
    for name in ["rna", "atac", "multivi"]:
        key = f"X_umap_{name}"
        if key in container.obsm:
            metadata[f"UMAP1_{name}"] = container.obsm[key][:, 0]
            metadata[f"UMAP2_{name}"] = container.obsm[key][:, 1]
    metadata.to_csv(output_dir / "cell_metadata.tsv.gz", sep="\t", compression="gzip")
