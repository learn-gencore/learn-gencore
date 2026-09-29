# Bulk RNA-seq analysis tutorial

A reproducible, plot-led bulk RNA-seq notebook using the public Day 0 mouse embryonic stem cell counts in GEO **GSE327304**. The worked comparisons are WT, beta-actin knockout (KO), and nuclear-localization-signal rescue (NLS). The notebook covers sample-level QC, normalization, RLE, PCA, variance and dispersion diagnostics, block and batch design, differential expression, effect sizes, and common pitfalls.

## Data and publication

- Dataset: [NCBI GEO GSE327304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327304)
- Source archive: [GSE327304_RAW.tar](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE327304&format=file)
- Linked manuscript/preprint: [Nuclear beta-actin-dependent chromatin accessibility governs stem cell pluripotency and extracellular matrix gene programs](https://doi.org/10.64898/2026.04.15.718829)

The checked-in `data/processed/GSE327304_D0_raw_counts.tsv` is an unnormalized gene-by-sample count matrix assembled from the 12 D0 sample count files in GEO. `data/metadata/GSE327304_D0_sample_manifest.csv` records sample accession, condition, time point, replicate label, and GEO archive member. The original small GEO archive can be re-downloaded with `python scripts/prepare_counts.py`; downloaded archives are ignored by Git. FASTQ files are not included.

**Sample-count note:** GEO currently lists four D0 libraries per condition, whereas the linked preprint methods appear to describe three biological replicates for RNA-seq. This tutorial uses all 12 D0 count files that GEO currently lists. Treat the output as a reproducible reanalysis until the intended manuscript sample subset and replicate structure are reconciled.

## Setup

Create the environment from the tutorial project root:

```bash
conda env create -f environment.yml
conda activate bulk-rnaseq-tutorial
```

If the processed count matrix is absent, prepare it from the linked GEO archive:

```bash
python scripts/prepare_counts.py
```

Open the notebook from the tutorial project root so its relative paths resolve:

```bash
jupyter lab notebooks/bulk_rnaseq_tutorial.ipynb
```

Run all cells from top to bottom. The DE model is fit with Bioconductor DESeq2 in R (`~ replicate + condition`); Python handles data checks, plots, and notebook presentation. `Rscript` must be available in the active environment. The replicate term is appropriate only if R1-R4 are true matched blocks across conditions; verify that against the experiment record.

## Outputs

The notebook saves figures as PNG files in `results/figures/`, sample and model result tables in `results/tables/`, and keeps the figures/results visible in the executed notebook. The main differential-expression tables are `KO_vs_WT_DESeq2.tsv`, `NLS_vs_KO_DESeq2.tsv`, and `NLS_vs_WT_DESeq2.tsv`.

The GEO archive contains processed gene counts, not FASTQ data. FastQC/MultiQC read quality, strandedness, mapping/quantification rate, and gene-body coverage need the raw reads and a separate workflow such as [nf-core/rnaseq](https://nf-co.re/rnaseq/3.27.0).

## References

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [tximport vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
- [edgeR User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
- [limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
- [variancePartition/DREAM vignette](https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html)

No software license is set in this scaffold. Choose a project license before publishing. The data remain attributed to the GEO submission and linked manuscript.
