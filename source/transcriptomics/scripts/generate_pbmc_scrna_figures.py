from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


OUTDIR = Path("source/img/transcriptomics/scrna")
OUTDIR.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


def mad(x):
    x = np.asarray(x)
    med = np.nanmedian(x)
    return np.nanmedian(np.abs(x - med))


def mad_bounds(x, nmads=3, lower=True, upper=True):
    med = np.nanmedian(x)
    m = mad(x)
    lo = med - nmads * m if lower else -np.inf
    hi = med + nmads * m if upper else np.inf
    return lo, hi


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


def add_qc_flags(adata):
    gene_lo, gene_hi = mad_bounds(adata.obs["n_genes_by_counts"], nmads=3)
    count_lo, count_hi = mad_bounds(adata.obs["total_counts"], nmads=3)
    _, mt_hi = mad_bounds(adata.obs["pct_counts_mt"], nmads=3, lower=False)
    complexity_lo, _ = mad_bounds(
        adata.obs["genes_per_1k_counts"],
        nmads=3,
        lower=True,
        upper=False,
    )

    adata.obs["qc_low_genes"] = adata.obs["n_genes_by_counts"] < gene_lo
    adata.obs["qc_high_genes"] = adata.obs["n_genes_by_counts"] > gene_hi
    adata.obs["qc_low_counts"] = adata.obs["total_counts"] < count_lo
    adata.obs["qc_high_counts"] = adata.obs["total_counts"] > count_hi
    adata.obs["qc_high_mt"] = adata.obs["pct_counts_mt"] > mt_hi
    adata.obs["qc_high_counts_low_complexity"] = (
        adata.obs["qc_high_counts"]
        & (adata.obs["genes_per_1k_counts"] < complexity_lo)
    )

    # Keep the example fast and deterministic. The tutorial text discusses
    # Scrublet doublet detection; these public PBMC figures focus on QC context.
    adata.obs["doublet_score"] = 0.0
    adata.obs["predicted_doublet"] = False
    adata.obs["qc_predicted_doublet"] = False

    qc_flag_columns = [
        "qc_low_genes",
        "qc_high_genes",
        "qc_low_counts",
        "qc_high_counts",
        "qc_high_mt",
        "qc_high_counts_low_complexity",
        "qc_predicted_doublet",
    ]
    adata.obs["qc_flag"] = adata.obs[qc_flag_columns].any(axis=1)
    adata.obs["qc_flag_count"] = adata.obs[qc_flag_columns].sum(axis=1)


def run_representation(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0_5")
    sc.tl.umap(adata)


def main():
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.match(r"^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    adata.obs["genes_per_1k_counts"] = (
        1000
        * adata.obs["n_genes_by_counts"]
        / adata.obs["total_counts"].replace(0, np.nan)
    )
    add_qc_flags(adata)

    sc.pl.violin(
        adata,
        [
            "total_counts",
            "n_genes_by_counts",
            "genes_per_1k_counts",
            "pct_counts_mt",
        ],
        multi_panel=True,
        jitter=0.4,
        show=False,
    )
    plt.savefig(OUTDIR / "pbmc_qc_distributions.png", bbox_inches="tight")
    plt.close("all")

    adata_work = adata.copy()
    min_cells = max(3, int(np.ceil(0.001 * adata_work.n_obs)))
    sc.pp.filter_genes(adata_work, min_cells=min_cells)
    adata_work.layers["counts"] = adata_work.X.copy()
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    adata_work.layers["norm"] = adata_work.X.copy()
    sc.pp.log1p(adata_work)
    adata_work.layers["lognorm"] = adata_work.X.copy()
    score_cell_cycle(adata_work)
    run_representation(adata_work)

    sc.pl.umap(
        adata_work,
        color=[
            "leiden_0_5",
            "qc_flag",
            "n_genes_by_counts",
            "pct_counts_mt",
            "genes_per_1k_counts",
            "phase",
        ],
        ncols=3,
        show=False,
    )
    plt.savefig(OUTDIR / "pbmc_first_pass_umap_qc.png", bbox_inches="tight")
    plt.close("all")

    adata.obs["remove_after_qc_review"] = (
        adata.obs["qc_low_genes"]
        | adata.obs["qc_low_counts"]
        | adata.obs["qc_high_counts_low_complexity"]
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
    run_representation(adata_qc)

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


if __name__ == "__main__":
    main()
