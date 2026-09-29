# Public PBMC scATAC tutorial through PEAKVI

This folder contains the runnable material for the Learn GenCore scATAC-seq
tutorial. It uses the public 10x Genomics 10k Human PBMCs ATAC v2 dataset and
continues through PEAKVI latent-space construction, UMAP, Leiden clustering,
and a cluster-marker accessibility example.

The workflow preserves raw peak counts, keeps QC and modeling separate, trains
PEAKVI, and audits the latent space for remaining technical structure. The
tutorial uses only public data and calls candidate accessible peaks for the
largest quality-reviewed cluster versus the remaining nuclei.

## Reproduce the tutorial

Create the environment:

```bash
conda env create -f environment.yml
conda activate gencore-scatac-peakvi
```

The shell scripts use Bash. Their command-line dependencies (`curl`, `gzip`,
and `awk`) are included in the Conda environment. Download checks use
`sha256sum`, `shasum`, or the environment's Python interpreter.

Download the public data and calculate fragment-length metrics:

```bash
scripts/download_data.sh
scripts/prepare_fragment_metrics.sh
```

Run the complete analysis:

```bash
python scripts/run_peakvi.py
```

The rendered notebook reloads the model, training history, and cache manifest
created by that command. To train inside the notebook instead, set
`RUN_TRAINING = True` in its PEAKVI cell.

The default executable example retains the 50,000 most variable accessible
regions after QC so it can be rendered on a CPU. Set `--max-peaks 0` in the
notebook implementation or pass a larger value in a research analysis when a
GPU is available. The official scvi-tools documentation notes that a GPU is
effectively required for fast PEAKVI inference on large datasets.

## Repository contents

- `notebooks/scatac_peakvi_tutorial.ipynb`: rendered teaching notebook
- `scripts/download_data.sh`: public 10x downloads
- `scripts/prepare_fragment_metrics.sh`: nucleosome signal and fragment-size QC
- `scripts/peakvi_workflow.py`: reusable analysis and plotting functions
- `scripts/run_peakvi.py`: command-line reproduction of the notebook
- `scripts/build_notebook.py`: regenerates the clean notebook before execution
- `checksums.sha256`: integrity checks for the downloaded 10x files
- `data/raw/`: downloaded inputs, ignored by Git
- `data/processed/`: generated fragment metrics, ignored by Git
- `results/`: trained model and compact results, ignored by Git

The differential-accessibility table is written to
`results/peakvi_da_cluster_0_vs_rest.tsv.gz` in the rendered example. It tests
only the 50,000 peaks retained for the teaching model. The calls describe
cluster-associated accessibility within one donor and do not support
condition-level inference.

## Dataset and citation

Dataset page:
https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-x-2-standard

The dataset is licensed under CC BY 4.0 and was analyzed by 10x Genomics with
Cell Ranger ATAC 2.1.0. The tutorial performs a secondary analysis and records
the software and filtering choices used here.
