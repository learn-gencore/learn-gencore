#!/usr/bin/env python3
"""Retrain MultiVI with other seeds and report which joint-cluster features recur.

Run after scripts/run_multivi.py. Each seed retrains the same model input
(same nuclei, genes, peaks, and settings); only the training seed changes. The
neighbor graph and Leiden clustering keep seed 17, so differences come from the
model fit.

    python scripts/seed_check.py --seed 23
    python scripts/seed_check.py --seed 42
    python scripts/seed_check.py --seed 101
    python scripts/seed_check.py --combine

`--combine` reloads every saved model, including the main one, and scores each
joint clustering against the main run's labels.
"""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp/learn_gencore_multiome_xdg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/learn_gencore_multiome_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/learn_gencore_multiome_numba")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import anndata as ad
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from multiome_workflow import (
    add_qc_flags,
    best_cluster_share,
    build_embedding_graph,
    build_multivi_input,
    check_observed_libraries,
    choose_qc_thresholds,
    load_multiome_dataset,
    load_multivi,
    score_doublets,
    select_atac_features,
    select_rna_features,
    small_group_masks,
    subset_paired,
    train_multivi,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--combine", action="store_true")
    return parser.parse_args()


def build_input(tutorial_dir: Path, results_dir: Path):
    rna, atac, _ = load_multiome_dataset(tutorial_dir / "data" / "raw", tutorial_dir / "data" / "processed")
    score_doublets(rna, seed=17)
    add_qc_flags(rna.obs, choose_qc_thresholds(rna.obs))
    atac.obs = rna.obs.copy()
    rna_paired, atac_paired = subset_paired(rna, atac)
    genes = select_rna_features(rna_paired, n_top_genes=3_000)
    peaks = select_atac_features(atac_paired, min_cells_per_peak=30, max_peaks=30_000)
    expected = pd.read_csv(results_dir / "model_genes.tsv", sep="\t")["gene"]
    if not pd.Index(genes).equals(pd.Index(expected)):
        raise ValueError("Model genes differ from the main run. Rerun scripts/run_multivi.py first.")
    mdata = build_multivi_input(rna_paired, atac_paired, genes, peaks)
    check_observed_libraries(mdata)
    return mdata


def recurrence_metrics(clusters: pd.Series, reference: pd.DataFrame, groups: dict[str, pd.Series]) -> dict:
    """Summaries of one joint clustering, measured against the main run's labels."""
    obs = reference.assign(cluster=clusters.astype(str).to_numpy())

    # CD14 monocytes (main-run labels): ratio of the highest to the lowest
    # median ATAC fragment count among joint clusters holding >= 50 of them.
    monocytes = obs.loc[obs["population"].eq("CD14 monocytes")]
    sizes = monocytes["cluster"].value_counts()
    large = sizes.index[sizes >= 50]
    medians = monocytes.loc[monocytes["cluster"].isin(large)].groupby("cluster")["atac_fragments"].median()

    row = {
        "n_joint_clusters": int(obs["cluster"].nunique()),
        "ari_vs_main_joint": adjusted_rand_score(obs["leiden_multivi"].astype(str), obs["cluster"]),
        "ari_vs_rna_only": adjusted_rand_score(obs["leiden_rna"].astype(str), obs["cluster"]),
        "ari_vs_atac_only": adjusted_rand_score(obs["leiden_atac"].astype(str), obs["cluster"]),
        "cd14_clusters_ge50": int(len(large)),
        "cd14_depth_ratio_max_min": float(medians.max() / medians.min()) if len(medians) > 1 else 1.0,
    }
    for name, mask in groups.items():
        short = name.split(" ")[0].lower().replace("-like", "")
        row[f"{short}_share_of_best_cluster"] = best_cluster_share(mask, obs["cluster"])["share_of_cluster"]
    return row


def main() -> None:
    args = parse_args()
    tutorial_dir = Path(__file__).resolve().parent.parent
    results_dir = tutorial_dir / "results"
    output_dir = results_dir / "seed_check"
    output_dir.mkdir(parents=True, exist_ok=True)
    mdata = build_input(tutorial_dir, results_dir)

    if not args.combine:
        model = train_multivi(mdata, model_dir=output_dir / f"model_seed_{args.seed}", seed=args.seed)
        epochs = int(len(model.history["reconstruction_loss_train"]))
        pd.DataFrame([{"seed": args.seed, "epochs_trained": epochs}]).to_csv(output_dir / f"seed_{args.seed}.tsv", sep="\t", index=False)
        print(f"Seed {args.seed}: trained {epochs} epochs")
        return

    reference = pd.read_csv(results_dir / "cell_metadata.tsv.gz", sep="\t", index_col=0).loc[mdata.obs_names]
    groups = small_group_masks(
        reference,
        pd.read_csv(results_dir / "dc_markers_by_rna_cluster.tsv", sep="\t", index_col=0),
        pd.read_csv(results_dir / "t_cell_state_markers_by_rna_cluster.tsv", sep="\t", index_col=0),
    )
    model_dirs = {17: results_dir / "models" / "multivi_pbmc_paired"}
    for path in output_dir.glob("model_seed_*"):
        model_dirs[int(re.search(r"(\d+)$", path.name).group(1))] = path

    rows = []
    for seed, model_dir in sorted(model_dirs.items()):
        model = load_multivi(mdata, model_dir)
        container = ad.AnnData(obs=reference[[]].copy())
        container.obsm["X_multivi"] = mdata.obsm["X_multivi"]
        build_embedding_graph(container, "X_multivi", "multivi", resolution=0.8, seed=17)
        history = getattr(model, "history", None) or {}
        rows.append(
            {
                "seed": seed,
                "epochs_trained": int(len(history["reconstruction_loss_train"])) if "reconstruction_loss_train" in history else None,
                **recurrence_metrics(container.obs["leiden_multivi"], reference, groups),
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(results_dir / "seed_check.tsv", sep="\t", index=False)
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
