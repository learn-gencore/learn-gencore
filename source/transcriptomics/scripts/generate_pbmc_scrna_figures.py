from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc


OUTDIR = Path("source/img/transcriptomics/scrna")
OUTDIR.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
sc.settings.datasetdir = "/private/tmp/scanpy-data"
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)

FLAG_PALETTE = {
    "pass": "#1f77b4",
    "high_counts_low_genes": "#d62728",
    "low_genes": "#ff7f0e",
    "high_counts": "#9467bd",
    "high_genes": "#2ca02c",
    "low_counts": "#8c564b",
    "high_mt": "#e377c2",
    "predicted_doublet": "#7f7f7f",
}


def percentile_thresholds(adata, lower_q=2.5, upper_q=97.5):
    gene_lo, gene_hi = np.nanpercentile(
        adata.obs["n_genes_by_counts"],
        [lower_q, upper_q],
    )
    count_lo, count_hi = np.nanpercentile(
        adata.obs["total_counts"],
        [lower_q, upper_q],
    )
    _, mt_hi = np.nanpercentile(
        adata.obs["pct_counts_mt"],
        [lower_q, upper_q],
    )
    return {
        "lower_q": lower_q,
        "upper_q": upper_q,
        "gene_lo": gene_lo,
        "gene_hi": gene_hi,
        "count_lo": count_lo,
        "count_hi": count_hi,
        "mt_hi": mt_hi,
    }


def score_cell_cycle(adata):
    s_genes = [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
        "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU",
        "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
        "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
        "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
        "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
    ]
    g2m_genes = [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
        "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
        "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
        "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
        "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
        "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
        "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
        "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ]
    s_present = [gene for gene in s_genes if gene in adata.var_names]
    g2m_present = [gene for gene in g2m_genes if gene in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_present, g2m_genes=g2m_present)


def encode_qc_flags(row):
    labels = []

    if row["qc_high_counts_low_genes"]:
        labels.append("high_counts_low_genes")
    else:
        if row["qc_low_genes"]:
            labels.append("low_genes")
        if row["qc_high_counts"]:
            labels.append("high_counts")

    if row["qc_high_genes"]:
        labels.append("high_genes")
    if row["qc_low_counts"]:
        labels.append("low_counts")
    if row["qc_high_mt"]:
        labels.append("high_mt")
    if row["qc_predicted_doublet"]:
        labels.append("predicted_doublet")

    return ";".join(labels) if labels else "pass"


def add_qc_flags(adata):
    thresholds = percentile_thresholds(adata)
    adata.uns["qc_thresholds"] = thresholds

    adata.obs["qc_low_genes"] = (
        adata.obs["n_genes_by_counts"] < thresholds["gene_lo"]
    )
    adata.obs["qc_high_genes"] = (
        adata.obs["n_genes_by_counts"] > thresholds["gene_hi"]
    )
    adata.obs["qc_low_counts"] = adata.obs["total_counts"] < thresholds["count_lo"]
    adata.obs["qc_high_counts"] = adata.obs["total_counts"] > thresholds["count_hi"]
    adata.obs["qc_high_mt"] = adata.obs["pct_counts_mt"] > thresholds["mt_hi"]
    adata.obs["qc_high_counts_low_genes"] = (
        adata.obs["qc_high_counts"]
        & adata.obs["qc_low_genes"]
    )

    # PBMC example figures explicitly omit Scrublet doublet detection.
    # Real analyses should run Scrublet and inspect these columns.
    adata.obs["doublet_score"] = 0.0
    adata.obs["predicted_doublet"] = False
    adata.obs["qc_predicted_doublet"] = False

    qc_flag_columns = [
        "qc_low_genes",
        "qc_high_genes",
        "qc_low_counts",
        "qc_high_counts",
        "qc_high_mt",
        "qc_high_counts_low_genes",
        "qc_predicted_doublet",
    ]
    adata.obs["qc_flag"] = adata.obs[qc_flag_columns].any(axis=1)
    adata.obs["qc_flag_count"] = adata.obs[qc_flag_columns].sum(axis=1)
    adata.obs["qc_flags"] = adata.obs[qc_flag_columns].apply(
        encode_qc_flags,
        axis=1,
    )
    adata.obs["primary_qc_flag"] = (
        adata.obs["qc_flags"]
        .str.split(";")
        .str[0]
        .astype("category")
    )


def save_scanpy_plot(filename):
    plt.savefig(OUTDIR / filename, bbox_inches="tight")
    plt.close("all")


def plot_percentile_thresholds(adata):
    t = adata.uns["qc_thresholds"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

    scatter = axes[0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        s=8,
        linewidths=0,
        alpha=0.75,
    )
    axes[0].axvline(t["count_lo"], color="black", linestyle="--", linewidth=1)
    axes[0].axvline(t["count_hi"], color="black", linestyle="--", linewidth=1)
    axes[0].axhline(t["gene_lo"], color="black", linestyle="--", linewidth=1)
    axes[0].axhline(t["gene_hi"], color="black", linestyle="--", linewidth=1)
    axes[0].set_xlabel("total_counts")
    axes[0].set_ylabel("n_genes_by_counts")
    axes[0].set_title("2.5/97.5 percentile guides")
    fig.colorbar(scatter, ax=axes[0], label="pct_counts_mt")

    axes[1].scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        s=8,
        linewidths=0,
        alpha=0.75,
    )
    axes[1].axvline(t["count_lo"], color="black", linestyle="--", linewidth=1)
    axes[1].axvline(t["count_hi"], color="black", linestyle="--", linewidth=1)
    axes[1].axhline(t["mt_hi"], color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("total_counts")
    axes[1].set_ylabel("pct_counts_mt")
    axes[1].set_title("MT upper guide")

    fig.savefig(OUTDIR / "pbmc_percentile_thresholds.png", bbox_inches="tight")
    plt.close(fig)


def plot_cell_cycle_scores(adata):
    sc.pl.violin(
        adata,
        ["S_score", "G2M_score"],
        groupby="phase",
        multi_panel=True,
        jitter=0.4,
        show=False,
    )
    save_scanpy_plot("pbmc_cell_cycle_scores.png")


def plot_pca_elbow(adata, filename, n_pcs=50):
    variance_ratio = adata.uns["pca"]["variance_ratio"][:n_pcs]
    pcs = np.arange(1, len(variance_ratio) + 1)

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(pcs, variance_ratio, marker="o", markersize=3, linewidth=1)
    ax.set_yscale("log")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance ratio, log scale")
    ax.set_title("PCA elbow plot")
    fig.savefig(OUTDIR / filename, bbox_inches="tight")
    plt.close(fig)


def plot_flagged_cells(adata):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)

    primary_flags = adata.obs["primary_qc_flag"].astype(str)
    colors = primary_flags.map(FLAG_PALETTE).fillna("#7f7f7f")
    axes[0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=colors,
        s=8,
        linewidths=0,
        alpha=0.75,
    )
    axes[0].set_xlabel("total_counts")
    axes[0].set_ylabel("n_genes_by_counts")
    axes[0].set_title("Primary QC flag")
    present_flags = [
        flag for flag in FLAG_PALETTE if flag in set(primary_flags)
    ]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="white",
            markerfacecolor=FLAG_PALETTE[flag],
            markersize=5,
            label=flag,
        )
        for flag in present_flags
    ]
    axes[0].legend(handles=handles, frameon=False, fontsize=7, loc="best")

    colors = adata.obs["qc_high_mt"].map({False: "#1f77b4", True: "#ff7f0e"})
    axes[1].scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        c=colors,
        s=8,
        linewidths=0,
        alpha=0.75,
    )
    axes[1].set_xlabel("total_counts")
    axes[1].set_ylabel("pct_counts_mt")
    axes[1].set_title("High MT flag")

    colors = adata.obs["qc_high_counts_low_genes"].map(
        {False: "#1f77b4", True: "#ff7f0e"}
    )
    axes[2].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=colors,
        s=8,
        linewidths=0,
        alpha=0.75,
    )
    axes[2].set_xlabel("total_counts")
    axes[2].set_ylabel("n_genes_by_counts")
    axes[2].set_title("High count, low gene flag")

    fig.savefig(OUTDIR / "pbmc_flagged_cells.png", bbox_inches="tight")
    plt.close(fig)


def run_representation(adata, n_pcs=30, save_diagnostics=False):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

    if save_diagnostics:
        sc.pl.highly_variable_genes(adata, show=False)
        save_scanpy_plot("pbmc_highly_variable_genes.png")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", mask_var="highly_variable")

    if save_diagnostics:
        plot_pca_elbow(adata, "pbmc_pca_elbow.png", n_pcs=50)

        sc.pl.pca(
            adata,
            color=[
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "S_score",
                "G2M_score",
                "phase",
            ],
            ncols=3,
            show=False,
        )
        save_scanpy_plot("pbmc_pca_qc.png")

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0_5")
    sc.tl.umap(adata)


def plot_qc_flags_by_cluster(adata):
    qc_primary_by_cluster = pd.crosstab(
        adata.obs["leiden_0_5"],
        adata.obs["primary_qc_flag"],
        normalize="index",
    )
    ordered_columns = [
        flag for flag in FLAG_PALETTE if flag in qc_primary_by_cluster.columns
    ]
    qc_primary_by_cluster = qc_primary_by_cluster[ordered_columns]

    ax = qc_primary_by_cluster.plot(
        kind="bar",
        stacked=True,
        figsize=(8, 4),
        width=0.85,
        color=[FLAG_PALETTE[flag] for flag in qc_primary_by_cluster.columns],
    )
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("Fraction of cells")
    ax.set_title("Primary QC flag composition by cluster")
    ax.legend(
        title="Primary QC flag",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    ax.figure.tight_layout()
    ax.figure.savefig(OUTDIR / "pbmc_qc_flags_by_cluster.png", bbox_inches="tight")
    plt.close(ax.figure)


def assign_broad_cell_types(adata, marker_genes):
    score_to_cell_type = {}

    for cell_type, genes in marker_genes.items():
        genes_present = [gene for gene in genes if gene in adata.raw.var_names]
        if len(genes_present) == 0:
            continue

        score_name = (
            "score_"
            + cell_type.lower()
            .replace(" ", "_")
            .replace("/", "_")
        )
        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=score_name,
            random_state=0,
            use_raw=True,
        )
        score_to_cell_type[score_name] = cell_type

    score_columns = list(score_to_cell_type)
    cluster_scores = (
        adata.obs
        .groupby("leiden_0_5", observed=True)[score_columns]
        .mean()
    )
    cluster_to_cell_type = (
        cluster_scores
        .idxmax(axis=1)
        .map(score_to_cell_type)
    )

    adata.obs["cell_type"] = (
        adata.obs["leiden_0_5"]
        .map(cluster_to_cell_type)
        .astype("category")
    )


def plot_t_cell_subclustering(adata_qc):
    broad_class = "T cells"
    is_t_cell = adata_qc.obs["cell_type"].astype(str).eq(broad_class)

    if int(is_t_cell.sum()) < 50:
        return

    adata_sub = adata_qc[is_t_cell].copy()
    adata_sub.X = adata_sub.layers["lognorm"].copy()

    sc.pp.highly_variable_genes(
        adata_sub,
        n_top_genes=min(1000, adata_sub.n_vars),
        flavor="seurat",
    )
    sc.pp.scale(adata_sub, max_value=10)
    sc.tl.pca(adata_sub, svd_solver="arpack", mask_var="highly_variable")

    plot_pca_elbow(adata_sub, "pbmc_t_cell_subclustering_elbow.png", n_pcs=50)

    sub_n_pcs = min(20, adata_sub.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=sub_n_pcs)
    sc.tl.leiden(adata_sub, resolution=0.5, key_added="subcluster")
    sc.tl.umap(adata_sub)

    sc.pl.umap(
        adata_sub,
        color=["subcluster", "primary_qc_flag", "phase"],
        ncols=3,
        show=False,
    )
    save_scanpy_plot("pbmc_t_cell_subclustering_umap.png")

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="subcluster",
        method="wilcoxon",
        use_raw=True,
    )
    sc.pl.rank_genes_groups(
        adata_sub,
        n_genes=15,
        sharey=False,
        show=False,
    )
    save_scanpy_plot("pbmc_t_cell_subclustering_rank_genes.png")

    adata_qc.obs["subcluster_label"] = pd.NA
    adata_qc.obs.loc[adata_sub.obs_names, "subcluster_label"] = (
        broad_class + "_" + adata_sub.obs["subcluster"].astype(str)
    )

    sc.pl.umap(
        adata_qc,
        color=["cell_type", "subcluster_label"],
        show=False,
    )
    save_scanpy_plot("pbmc_t_cell_subcluster_labels_full_umap.png")


def main():
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    add_qc_flags(adata)

    sc.pl.violin(
        adata,
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
        ],
        multi_panel=True,
        jitter=0.4,
        show=False,
    )
    plt.savefig(OUTDIR / "pbmc_qc_distributions.png", bbox_inches="tight")
    plt.close("all")
    plot_percentile_thresholds(adata)
    plot_flagged_cells(adata)

    adata_work = adata.copy()
    min_cells = max(3, int(np.ceil(0.001 * adata_work.n_obs)))
    sc.pp.filter_genes(adata_work, min_cells=min_cells)
    adata_work.layers["counts"] = adata_work.X.copy()
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    adata_work.layers["norm"] = adata_work.X.copy()
    sc.pp.log1p(adata_work)
    adata_work.layers["lognorm"] = adata_work.X.copy()
    score_cell_cycle(adata_work)
    plot_cell_cycle_scores(adata_work)
    run_representation(adata_work, n_pcs=30, save_diagnostics=True)

    sc.pl.umap(
        adata_work,
        color=[
            "leiden_0_5",
            "primary_qc_flag",
            "qc_flag_count",
            "n_genes_by_counts",
            "pct_counts_mt",
            "phase",
        ],
        ncols=3,
        show=False,
    )
    plt.savefig(OUTDIR / "pbmc_first_pass_umap_qc.png", bbox_inches="tight")
    plt.close("all")
    plot_qc_flags_by_cluster(adata_work)

    adata.obs["remove_after_qc_review"] = (
        adata.obs["qc_low_genes"]
        | adata.obs["qc_low_counts"]
        | adata.obs["qc_high_counts_low_genes"]
    )
    adata_qc = adata[~adata.obs["remove_after_qc_review"]].copy()
    min_cells = max(3, int(np.ceil(0.001 * adata_qc.n_obs)))
    sc.pp.filter_genes(adata_qc, min_cells=min_cells)
    adata_qc.layers["counts"] = adata_qc.X.copy()
    sc.pp.normalize_total(adata_qc, target_sum=1e4)
    adata_qc.layers["norm"] = adata_qc.X.copy()
    sc.pp.log1p(adata_qc)
    adata_qc.layers["lognorm"] = adata_qc.X.copy()
    adata_qc.raw = adata_qc
    score_cell_cycle(adata_qc)
    run_representation(adata_qc, n_pcs=30)

    sc.pl.umap(
        adata_qc,
        color=[
            "leiden_0_5",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "S_score",
            "G2M_score",
            "phase",
        ],
        ncols=3,
        show=False,
    )
    save_scanpy_plot("pbmc_final_umap_qc.png")
    adata_qc.X = adata_qc.layers["lognorm"].copy()

    marker_genes = {
        "T cells": ["CD3D", "CD3E", "TRAC"],
        "B cells": ["MS4A1", "CD79A", "CD74"],
        "NK cells": ["GNLY", "NKG7", "KLRD1"],
        "Monocytes": ["LYZ", "S100A8", "S100A9"],
        "Dendritic cells": ["FCER1A", "CST3"],
        "Platelets": ["PPBP", "PF4"],
    }
    marker_genes = {
        name: [gene for gene in genes if gene in adata_qc.var_names]
        for name, genes in marker_genes.items()
    }
    marker_genes = {name: genes for name, genes in marker_genes.items() if genes}

    dotplot = sc.pl.dotplot(
        adata_qc,
        marker_genes,
        groupby="leiden_0_5",
        standard_scale="var",
        dendrogram=False,
        use_raw=True,
        show=False,
        return_fig=True,
    )
    dotplot.savefig(OUTDIR / "pbmc_marker_dotplot.png", bbox_inches="tight")
    plt.close("all")

    genes_to_plot = ["CD3D", "MS4A1", "LYZ", "NKG7", "PPBP"]
    genes_to_plot = [gene for gene in genes_to_plot if gene in adata_qc.var_names]
    sc.pl.umap(
        adata_qc,
        color=genes_to_plot,
        vmax="p99",
        use_raw=True,
        ncols=3,
        show=False,
    )
    save_scanpy_plot("pbmc_marker_umaps.png")

    sc.tl.rank_genes_groups(
        adata_qc,
        groupby="leiden_0_5",
        method="wilcoxon",
        use_raw=True,
    )
    sc.pl.rank_genes_groups(
        adata_qc,
        n_genes=15,
        sharey=False,
        show=False,
    )
    save_scanpy_plot("pbmc_rank_genes_groups.png")

    assign_broad_cell_types(adata_qc, marker_genes)
    sc.pl.umap(
        adata_qc,
        color=["leiden_0_5", "cell_type"],
        show=False,
    )
    save_scanpy_plot("pbmc_cell_type_umap.png")

    plot_t_cell_subclustering(adata_qc)


if __name__ == "__main__":
    main()
