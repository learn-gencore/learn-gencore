#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_deseq2.R <count_matrix.tsv> <sample_manifest.csv> <output_directory>")
}
count_path <- args[[1]]
manifest_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

count_tab <- read.delim(count_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_id" %in% colnames(count_tab)) stop("Count matrix must contain a gene_id column")
rownames(count_tab) <- count_tab$gene_id
count_tab$gene_id <- NULL
count_matrix <- as.matrix(count_tab)
storage.mode(count_matrix) <- "integer"

samples <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(samples) <- samples$sample_id
samples$condition <- factor(samples$condition, levels = c("WT", "KO", "NLS"))
samples$condition <- relevel(samples$condition, ref = "WT")
samples$replicate <- factor(samples$replicate, levels = c("R1", "R2", "R3", "R4"))
if (!all(rownames(samples) %in% colnames(count_matrix))) stop("A manifest sample is missing from the count matrix")
count_matrix <- count_matrix[, rownames(samples), drop = FALSE]
if (anyNA(count_matrix) || any(count_matrix < 0L)) stop("Counts must be non-negative integers")

keep <- rowSums(count_matrix >= 10L) >= 4L
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix[keep, , drop = FALSE],
  colData = samples,
  design = ~ replicate + condition
)
dds <- DESeq(dds, quiet = TRUE)

write_table <- function(x, path) {
  tab <- data.frame(gene_id = rownames(x), x, check.names = FALSE)
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_table(counts(dds, normalized = TRUE), file.path(out_dir, "normalized_counts.tsv"))
vst_obj <- vst(dds, blind = FALSE)
write_table(assay(vst_obj), file.path(out_dir, "vst_counts.tsv"))
write.table(
  data.frame(sample_id = names(sizeFactors(dds)), size_factor = unname(sizeFactors(dds))),
  file.path(out_dir, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    gene_id = rownames(dds),
    baseMean = mcols(dds)$baseMean,
    dispGeneEst = mcols(dds)$dispGeneEst,
    dispFit = mcols(dds)$dispFit,
    dispersion = dispersions(dds)
  ),
  file.path(out_dir, "dispersion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

contrasts <- list(
  KO_vs_WT = c("condition", "KO", "WT"),
  NLS_vs_KO = c("condition", "NLS", "KO"),
  NLS_vs_WT = c("condition", "NLS", "WT")
)
summary_rows <- list()
for (nm in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[nm]], alpha = 0.05)
  tab <- as.data.frame(res)
  tab$gene_id <- rownames(tab)
  tab <- tab[, c("gene_id", setdiff(colnames(tab), "gene_id"))]
  tab <- tab[order(tab$pvalue), ]
  write.table(tab, file.path(out_dir, paste0(nm, "_DESeq2.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  summary_rows[[nm]] <- data.frame(
    contrast = nm,
    tested_genes = sum(!is.na(tab$pvalue)),
    FDR_lt_0_05 = sum(!is.na(tab$padj) & tab$padj < 0.05),
    FDR_lt_0_10 = sum(!is.na(tab$padj) & tab$padj < 0.10)
  )
}
write.table(do.call(rbind, summary_rows), file.path(out_dir, "contrast_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("DESeq2", as.character(packageVersion("DESeq2")))), file.path(out_dir, "software_versions.txt"))
message("DESeq2 completed: ", sum(keep), " genes retained; result tables saved to the requested output directory.")
