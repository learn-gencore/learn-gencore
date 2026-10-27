# Raw public data

This directory is populated by `../../scripts/download_data.sh`. The downloaded
files are ignored by Git because the fragment file is approximately 540 MB.

Dataset: PBMC from a healthy donor, no cell sorting (3k), Chromium Single Cell
Multiome ATAC + Gene Expression, Cell Ranger ARC 2.0.0, GRCh38 reference
`refdata-cellranger-arc-GRCh38-2020-A-2.0.0`.

Official dataset page:
https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0

The tutorial uses these Cell Ranger ARC outputs, saved without the
`pbmc_unsorted_3k_` prefix:

- `filtered_feature_bc_matrix.h5`: RNA and peak counts for called nuclei
- `per_barcode_metrics.csv`: RNA and ATAC metrics for every barcode
- `atac_peaks.bed` and `atac_peak_annotation.tsv`
- `atac_fragments.tsv.gz` and `atac_fragments.tsv.gz.tbi`
- `summary.csv`

`checksums.sha256` was computed from the files downloaded for the rendered run.
It detects incomplete or changed downloads. It is not a checksum published by
10x Genomics.

Cite the dataset page and indicate that this tutorial performs a secondary
analysis of the published 10x output.
