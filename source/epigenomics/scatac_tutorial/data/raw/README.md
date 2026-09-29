# Raw public data

This directory is populated by `../../scripts/download_data.sh`. The downloaded
files are ignored by Git because the fragment file is approximately 2.3 GB.

Dataset: 10k Human PBMCs, ATAC v2, Chromium X, Cell Ranger ATAC 2.1.0.

Official dataset page:
https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-x-2-standard

The tutorial uses these Cell Ranger ATAC outputs:

- `filtered_peak_bc_matrix.h5`
- `filtered_tf_bc_matrix.h5`
- `singlecell.csv`
- `peak_annotation.tsv`
- `fragments.tsv.gz` and `fragments.tsv.gz.tbi`
- `summary.csv` and `summary.json`

The 10x dataset is licensed under CC BY 4.0. Cite the dataset page and indicate
that this tutorial performs a secondary analysis of the published 10x output.
