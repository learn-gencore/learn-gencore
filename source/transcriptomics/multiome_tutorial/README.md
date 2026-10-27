# Public PBMC multiome tutorial through paired MultiVI

This folder contains the runnable material for the Learn GenCore single-cell
multiome tutorial. It uses the public 10x Genomics PBMC 3k multiome dataset
(RNA and chromatin accessibility measured in the same nuclei). It continues
through joint quality control, single-modality baselines, a paired MultiVI
model, and marker evidence from both measurements.

The workflow keeps the two measurements separate during QC, preserves raw
counts for modeling, and audits the joint representation for technical
structure before interpreting clusters.

## Reproduce the tutorial

Create the environment:

```bash
conda env create -f environment.yml
conda activate gencore-multiome-multivi
```

Download the public data and calculate fragment-length metrics:

```bash
scripts/download_data.sh
scripts/prepare_fragment_metrics.sh
```

Run the complete analysis:

```bash
python scripts/run_multivi.py
```

The rendered notebook reloads the model, training history, and cache manifest
created by that command. To train inside the notebook instead, set
`RUN_TRAINING = True` in its MultiVI cell.

The scripts and notebook set `OMP_NUM_THREADS=1` and import Scanpy before
scvi-tools. On macOS the PyTorch and scikit-learn wheels each bundle an OpenMP
runtime. In the tutorial environment, MultiVI training stalled after Scanpy's
neighbor and UMAP steps unless OpenMP was limited to one thread, and importing
scvi-tools (and with it PyTorch) before Scanpy made numba's parallel QC code
crash. With both settings, the rendered model trained on a CPU.

The default model uses the 3,000 most variable genes, excluding mitochondrial
genes and their nuclear MTRNR2L copies, and the 30,000 peaks with the largest
binary accessibility variance among peaks detected in at least 30 nuclei.
Feature selection is a reported analysis choice because it changes what the
model can represent.

`scripts/seed_check.py` retrains the same model input with three other seeds
(up to about 15 minutes each on one CPU thread) and reports which features of the
joint clustering recur. Run it after `run_multivi.py`.

## Repository contents

- `notebooks/multiome_multivi_tutorial.ipynb`: rendered teaching notebook
- `scripts/download_data.sh`: public 10x downloads and checksum verification
- `scripts/prepare_fragment_metrics.sh`: nucleosome signal and fragment sizes
- `scripts/multiome_workflow.py`: reusable analysis and plotting functions
- `scripts/run_multivi.py`: command-line reproduction of the notebook
- `scripts/build_notebook.py`: regenerates the clean notebook before execution
- `scripts/seed_check.py`: retrains MultiVI with other seeds and reports which
  joint-cluster features recur (run after `run_multivi.py`)
- `checksums.sha256`: integrity checks for the downloaded 10x files
- `claims_ledger.md`: the result behind each statement on the tutorial page
- `data/raw/`: downloaded inputs, ignored by Git
- `data/processed/`: generated fragment metrics, ignored by Git
- `results/`: trained model and compact results, ignored by Git

## Dataset and citation

Dataset page:
https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0

The dataset was analyzed by 10x Genomics with Cell Ranger ARC 2.0.0 against
GRCh38 (`refdata-cellranger-arc-GRCh38-2020-A-2.0.0`). The tutorial performs a
secondary analysis and records the software and filtering choices used here.
