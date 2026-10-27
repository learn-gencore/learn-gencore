# Generated tutorial data

This directory contains small tables created from the raw 10x outputs. Run
`../../scripts/prepare_fragment_metrics.sh` before executing the notebook.

- `cell_barcodes.tsv`: barcodes called as cells by Cell Ranger ARC
- `fragment_metrics.tsv`: per-nucleus nucleosome-free and mononucleosomal
  fragment counts and nucleosome signal
- `fragment_length_histogram.tsv`: fragment-length distribution of called nuclei

Generated files are ignored by Git. The notebook records all thresholds and
recreates its plots and result tables from these inputs.
