# Data

`processed/GSE327304_D0_raw_counts.tsv` contains the original integer gene counts for the 12 D0 libraries used in the tutorial. It is built by concatenating the corresponding compressed two-column count files in the NCBI GEO supplementary archive for [GSE327304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327304).

`metadata/GSE327304_D0_sample_manifest.csv` records GEO sample IDs, condition labels, D0 time point, replicate label, and the source archive filename.

To recreate the matrix from the GEO archive, run `python scripts/prepare_counts.py` from the tutorial project root. The archive will be saved under `data/raw/` and is excluded from Git; the compact processed count matrix remains tracked. SRA FASTQ files are not downloaded by this script.
