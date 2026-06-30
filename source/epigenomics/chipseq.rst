ChIP-Seq Analysis: From Aligned BAM Files to Differential Binding & Enrichment
==============================================================================

Post-alignment QC, peak calling, differential binding analysis, and functional annotation of a **ChIP-seq** experiment, starting from already aligned (BAM) reads.

This workshop uses publicly available ENCODE ChIP-seq data for the transcription factor **CTCF** in two human cell lines, **GM12878** and **K562**, each with a matched input/control. Reads are already aligned (BWA, hg38), so the workflow picks up at quality control: signal QC with ``deeptools``, peak calling with ``MACS2``, differential binding between cell lines with ``DiffBind``/``DESeq2``, and annotation/enrichment with ``ChIPseeker`` and ``clusterProfiler``.

Scope is limited to a **single, guided end-to-end pass** comparing CTCF binding between two cell lines. It does not cover read trimming/alignment (assumed already done), motif discovery, or multi-factor analyses. The aim is a reference workflow you can adapt to your own ChIP-seq/CUT&RUN data.

.. note::

   Commands below assume a ``conda``-managed environment on a standalone machine or HPC cluster. SLURM directives are shown where relevant but are ignored if you run commands directly on the command line.

Pipeline overview
-----------------

::

      BAM (aligned reads)
             |
             v
     [1] Post-alignment QC -------- deeptools (bamCoverage, multiBamSummary,
             |                       plotFingerprint, plotCorrelation, plotPCA)
             v
     [2] Peak calling ------------- MACS2 (callpeak, vs matched input)
             |
             v
     [3] Consensus / count matrix - DiffBind (dba.count)
             |
             v
     [4] Differential binding ----- DiffBind + DESeq2 (dba.analyze)
             |
             v
     [5] Annotation ---------------- ChIPseeker (peak-to-gene, genomic context)
             |
             v
     [6] Enrichment analysis ------ clusterProfiler (GO / KEGG / Reactome)

The BAM files are the hub: QC tells you whether the IP worked at all, peak calling turns signal into discrete binding sites, differential analysis compares those sites between conditions, and annotation/enrichment turns genomic coordinates into biology.

Steps
-----

1. Assess IP quality and reproducibility directly from BAMs with **deeptools**.
2. Call peaks per replicate with **MACS2**, against matched input control.
3. Build a consensus peak set and count matrix with **DiffBind**.
4. Run differential binding analysis between GM12878 and K562 with **DiffBind/DESeq2**.
5. Annotate peaks to genes and genomic features with **ChIPseeker**.
6. Run GO/KEGG enrichment on the bound/differential genes with **clusterProfiler**.

Environment
-----------

This workshop spans two different environments:

- **Steps 1-2** (deeptools, MACS2) are command-line tools, run either via
  HPC environment modules or a conda environment -- see "Loading the
  necessary software environment" at the start of each of those steps.
- **Steps 3-6** (DiffBind, DESeq2, ChIPseeker, clusterProfiler) are R/
  Bioconductor analyses, run locally in RStudio on your own machine --
  see `Moving to a local RStudio session`_ at the start of Step 3.
  Package installation for these is embedded directly in each step's R
  code, so there's nothing to set up in advance.

For the command-line tools, versions below are pinned loosely; substitute
what's available to you, but verify flags against your installed version.

.. code:: bash

   # create one environment for the command-line tools
   conda create -n chipseq -c bioconda -c conda-forge \
       deeptools=3.5 macs2=2.2.9.1 samtools=1.20 bedtools=2.31 ucsc-bedclip ucsc-bedgraphtobigwig
   conda activate chipseq

.. note::

   ``deeptools``/``MACS2``/``samtools`` live happily in one conda environment. The R/Bioconductor side is handled separately, locally, in RStudio -- see Step 3.

The tools used across the workshop, with upstream documentation:

-  deeptools — https://deeptools.readthedocs.io/
-  MACS2 — https://github.com/macs3-project/MACS
-  DiffBind — https://bioconductor.org/packages/release/bioc/html/DiffBind.html
-  DESeq2 — https://bioconductor.org/packages/release/bioc/html/DESeq2.html
-  ChIPseeker — https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html
-  clusterProfiler — https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
-  samtools — https://www.htslib.org/
-  bedtools — https://bedtools.readthedocs.io/
-  R — https://cran.r-project.org/
-  RStudio Desktop — https://posit.co/download/rstudio-desktop/

Inputs
------

Before running the pipeline you need:

-  **Coordinate-sorted, indexed BAM files**, one per IP replicate plus a matched input/control, for each condition. This workshop uses 2 cell lines × (2 CTCF replicates + 1 input) = 6 BAM files, aligned to **GRCh38**.
-  A **sample sheet** (``samples.csv``) describing each BAM, its condition, replicate number, and which input it pairs with (format given in `Building the sample sheet <#building-the-sample-sheet>`__).
-  The **hg38 genome** (chromosome sizes file, used by deeptools) and a **TxDb annotation** (used by ChIPseeker), both obtained automatically by the R packages above — no separate download needed.

File formats you will encounter (assumed familiar): BAM, BED, narrowPeak, bigWig.

Getting the data
----------------

This tutorial uses ENCODE CTCF ChIP-seq BAM files for GM12878 and K562
(GRCh38), each with a matched input control. These are ENCODE2-era
experiments, so each biological replicate ships as two technical-replicate
BAM files (same biological replicate, different lanes/runs) that were
merged with ``samtools merge`` into a single BAM per biological replicate
before any downstream analysis, then re-sorted and indexed.

.. list-table:: ENCODE source files used in this tutorial
   :header-rows: 1
   :widths: 12 16 16 28 28

   * - Cell line
     - Experiment
     - Control experiment
     - IP BAM accessions (bio rep -> files, merged)
     - Input BAM accessions (merged)
   * - GM12878
     - `ENCSR000DZN <https://www.encodeproject.org/experiments/ENCSR000DZN/>`__
     - `ENCSR000EYX <https://www.encodeproject.org/experiments/ENCSR000EYX/>`__
     - rep1: ENCFF162QXM, ENCFF747TJH; rep2: ENCFF033WII, ENCFF876QCV
     - ENCFF157YWH, ENCFF658PBP (rep1)
   * - K562
     - `ENCSR000DWE <https://www.encodeproject.org/experiments/ENCSR000DWE/>`__
     - `ENCSR000DWA <https://www.encodeproject.org/experiments/ENCSR000DWA/>`__
     - rep1: ENCFF081HVQ, ENCFF800GVR; rep2: ENCFF494VZW, ENCFF196QVZ
     - ENCFF873PSH, ENCFF897QAS (rep1)

The result is exactly the six BAM files (plus ``.bai`` indexes) used
throughout the rest of this tutorial: ``bam/gm12878_ctcf_rep1.bam``,
``bam/gm12878_ctcf_rep2.bam``, ``bam/gm12878_input.bam``,
``bam/k562_ctcf_rep1.bam``, ``bam/k562_ctcf_rep2.bam``,
``bam/k562_input.bam``.

These pre-merged, pre-indexed BAM/BAI files are available directly from
this shared folder, so you don't need to repeat the download-and-merge
steps yourself:

`Download the merged BAM/BAI files <https://drive.google.com/drive/folders/1_XPBJ1b0VGs6ThYi4dyj2_RkOBXGyIfj?usp=sharing>`__

Download all six ``.bam`` files and their matching ``.bai`` indexes from
that folder into a local ``bam/`` directory before continuing to the QC
steps below.

If you'd rather resolve current accessions programmatically instead (e.g.
for a different factor or cell line), the script below queries the ENCODE
REST API directly:

.. literalinclude:: ../scripts/fetch_encode_bams.sh
   :language: bash
   :caption: scripts/fetch_encode_bams.sh

.. note::

   ENCODE BAMs ship coordinate-sorted but usually unindexed and sometimes without the conventional ``chr`` prefix matching your TxDb. Confirm with ``samtools view -H file.bam | grep ^@SQ`` before proceeding.

The data
--------

The dataset is paired-end CTCF ChIP-seq, GRCh38, two human cell lines:

======================= =========================== ====== ==========
Sample                  Cell line                   Factor Replicates
======================= =========================== ====== ==========
``gm12878_ctcf_rep1/2`` GM12878 (lymphoblastoid)    CTCF   2
``gm12878_input``       GM12878                     input  1
``k562_ctcf_rep1/2``    K562 (myelogenous leukemia) CTCF   2
``k562_input``          K562                        input  1
======================= =========================== ====== ==========

CTCF is a good first ChIP-seq factor to learn on: it gives sharp, abundant, well-characterized peaks (tens of thousands genome-wide), so QC and peak calling behave predictably and there’s a large literature to sanity-check results against.

Building the sample sheet
~~~~~~~~~~~~~~~~~~~~~~~~~

``samples.csv``, used by DiffBind in `Step 3 <#step-3-consensus-peak-set-and-count-matrix-diffbind>`__:

::

   SampleID,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
   GM12878_CTCF_1,GM12878,1,bam/gm12878_ctcf_rep1.bam,GM12878_input,bam/gm12878_input.bam,macs2/GM12878_rep1_peaks.narrowPeak,narrow
   GM12878_CTCF_2,GM12878,2,bam/gm12878_ctcf_rep2.bam,GM12878_input,bam/gm12878_input.bam,macs2/GM12878_rep2_peaks.narrowPeak,narrow
   K562_CTCF_1,K562,1,bam/k562_ctcf_rep1.bam,K562_input,bam/k562_input.bam,macs2/K562_rep1_peaks.narrowPeak,narrow
   K562_CTCF_2,K562,2,bam/k562_ctcf_rep2.bam,K562_input,bam/k562_input.bam,macs2/K562_rep2_peaks.narrowPeak,narrow

The ``Peaks``/``PeakCaller`` columns point at files produced in Step 2 — leave the sheet open and fill them in as you go.

Step 1: Post-alignment QC (deeptools)
-------------------------------------

Loading the necessary software environment (deeptools)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're working on an HPC cluster with environment modules (e.g. the
Jubail cluster at NYUAD), load the relevant module instead of activating
conda:

.. code:: bash

   module purge
   module load all
   module load gencore/3
   module load deeptools

Otherwise, activate the conda (or equivalent) environment containing this
step's tools (see Environment_ above):

.. code:: bash

   conda activate chipseq

Before calling peaks, confirm the IP actually enriched for your factor and
that replicates agree with each other.

If you downloaded the pre-merged BAM/BAI files from the shared folder
above, indexing is already done -- just confirm each ``.bai`` is present
and sits alongside its BAM:

.. code:: bash

   for b in bam/*.bam; do
       [ -f "${b}.bai" ] || echo "MISSING INDEX: ${b}"
   done

If you're working from your own BAMs instead (e.g. freshly merged or
downloaded without indexes), index them first:

.. code:: bash

   for b in bam/*.bam; do samtools index -@ 8 "$b"; done

Signal tracks
~~~~~~~~~~~~~

Convert each BAM to a normalized bigWig for visual inspection in IGV and as input to several deeptools QC plots:

.. code:: bash

   mkdir -p bigwig
   for b in bam/*.bam; do
     name=$(basename "$b" .bam)
     bamCoverage \
         -b "$b" \
         -o bigwig/${name}.bw \
         --binSize 10 \
         --normalizeUsing RPGC \
         --effectiveGenomeSize 2913022398 \
         --extendReads \
         -p 8
   done

Flags:

-  ``--normalizeUsing RPGC`` — normalize to 1x average genome coverage, making samples comparable regardless of sequencing depth.
-  ``--effectiveGenomeSize 2913022398`` — the GRCh38 effective genome size (non-N, mappable bases); see the deeptools table for other genomes/read lengths (https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html).
-  ``--extendReads`` — extend single-end reads to the estimated fragment size; for paired-end BAMs (as here) deeptools uses the actual fragment size from the read pairs instead.

IP strength: plotFingerprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``plotFingerprint`` is the single best deeptools plot for “did the ChIP work?” It plots, per BAM, the cumulative fraction of reads against the cumulative fraction of the genome. A successful IP concentrates reads into a small fraction of the genome (a curve that bows sharply away from the diagonal); input/failed IPs look close to a straight diagonal line.

.. code:: bash

   plotFingerprint \
       -b bam/gm12878_ctcf_rep1.bam bam/gm12878_ctcf_rep2.bam bam/gm12878_input.bam \
          bam/k562_ctcf_rep1.bam bam/k562_ctcf_rep2.bam bam/k562_input.bam \
       --labels GM_rep1 GM_rep2 GM_input K562_rep1 K562_rep2 K562_input \
       --extendReads \
       -p 8 \
       --plotFile fingerprint.png \
       --outRawCounts fingerprint.tab

Expect the two CTCF replicates per cell line to track closely together, both well separated from their matched input.

Replicate reproducibility: multiBamSummary + plotCorrelation / plotPCA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bin the genome and count reads per BAM per bin, then correlate and cluster samples:

.. code:: bash

   multiBamSummary bins \
       --bamfiles bam/gm12878_ctcf_rep1.bam bam/gm12878_ctcf_rep2.bam bam/gm12878_input.bam \
                  bam/k562_ctcf_rep1.bam bam/k562_ctcf_rep2.bam bam/k562_input.bam \
       --labels GM_rep1 GM_rep2 GM_input K562_rep1 K562_rep2 K562_input \
       --binSize 10000 \
       -p 8 \
       -o multiBamSummary.npz

   plotCorrelation \
       -in multiBamSummary.npz \
       --corMethod spearman \
       --whatToPlot heatmap \
       --plotFile correlation_heatmap.png \
       --plotNumbers

   plotPCA \
       -in multiBamSummary.npz \
       -o pca.png

Replicates of the same condition should correlate most strongly with each other and cluster together on the PCA plot; inputs should separate from IP samples on PC1 or PC2. If a replicate is an outlier here, treat its downstream peak calls with caution before you’ve even called peaks.

Strand cross-correlation (optional but informative)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``plotFingerprint`` and correlation plots cover most of what you need, but if you want classic ENCODE-style QC metrics (NSC/RSC), run ``phantompeakqualtools`` or the ``deeptools`` equivalent ``bamPEFragmentSize`` for paired-end fragment size sanity-checking:

.. code:: bash

   bamPEFragmentSize \
       -b bam/gm12878_ctcf_rep1.bam bam/k562_ctcf_rep1.bam \
       --histogram fragment_sizes.png \
       -p 8

Step 2: Peak calling (MACS2)
----------------------------

Loading the necessary software environment (MACS2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're working on an HPC cluster with environment modules (e.g. the
Jubail cluster at NYUAD), load the relevant module instead of activating
conda:

.. code:: bash

   module purge
   module load all
   module load gencore/3
   module load macs2

Otherwise, activate the conda (or equivalent) environment containing this
step's tools (see Environment_ above):

.. code:: bash

   conda activate chipseq

Call peaks per IP replicate against its matched input. CTCF gives sharp narrow peaks, so ``narrowPeak`` mode is appropriate (use ``--broad`` instead for diffuse histone marks like H3K27me3 or H3K36me3).

.. code:: bash

   mkdir -p macs2

   macs2 callpeak \
       -t bam/gm12878_ctcf_rep1.bam \
       -c bam/gm12878_input.bam \
       -f BAMPE \
       -g hs \
       -n GM12878_rep1 \
       --outdir macs2 \
       -q 0.05

   macs2 callpeak \
       -t bam/gm12878_ctcf_rep2.bam \
       -c bam/gm12878_input.bam \
       -f BAMPE -g hs -n GM12878_rep2 --outdir macs2 -q 0.05

   macs2 callpeak \
       -t bam/k562_ctcf_rep1.bam \
       -c bam/k562_input.bam \
       -f BAMPE -g hs -n K562_rep1 --outdir macs2 -q 0.05

   macs2 callpeak \
       -t bam/k562_ctcf_rep2.bam \
       -c bam/k562_input.bam \
       -f BAMPE -g hs -n K562_rep2 --outdir macs2 -q 0.05

Flags:

-  ``-f BAMPE`` — treat the BAM as paired-end and use actual fragment spans rather than estimating fragment length from a single-end model.
-  ``-g hs`` — effective genome size shortcut for human.
-  ``-q 0.05`` — minimum FDR (q-value) cutoff for peak calling; tighten for a more stringent peak set.

Each run produces ``<name>_peaks.narrowPeak`` (the file referenced in ``samples.csv``), ``<name>_summits.bed``, and ``<name>_peaks.xls``. Quick sanity check on peak counts:

.. code:: bash

   wc -l macs2/*_peaks.narrowPeak

CTCF in these cell lines typically yields on the order of 30,000–60,000 peaks per replicate; if you see far fewer, revisit the QC plots above before proceeding.

Step 3: Consensus peak set and count matrix (DiffBind)
------------------------------------------------------

Moving to a local RStudio session
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Steps 3 through 6 are R/Bioconductor analyses, intended to be run on your
own machine in RStudio rather than on the HPC cluster -- it's an
interactive, plot-heavy part of the workflow that's much nicer with a
local GUI than over a remote terminal.

1. **Install R** (if you don't already have it): `<https://cran.r-project.org/>`__
   -- pick the installer for your OS (Windows, macOS, or Linux).
2. **Install RStudio Desktop**: `<https://posit.co/download/rstudio-desktop/>`__
   -- free edition is sufficient.
3. **Copy the files you'll need onto your local machine**, into a single
   project folder (e.g. ``chipseq_workshop/``):

   - ``samples.csv`` (the DiffBind sample sheet from `Building the sample sheet`_)
   - the per-replicate peak calls from Step 2: ``macs2/*_peaks.narrowPeak``
   - the BAM/BAI files from `Getting the data`_: ``bam/*.bam`` and ``bam/*.bam.bai``

   Use ``scp``/``rsync`` from the HPC, or download directly from the
   shared folder if that's where your BAMs live. Update the paths in
   ``samples.csv`` to match wherever you put these files locally.
4. **Open RStudio**, set your working directory to that project folder
   (Session -> Set Working Directory -> Choose Directory, or
   ``setwd("path/to/chipseq_workshop")`` in the console), and run the R
   code blocks below directly in the console or as an R script.

Each R code block below installs any Bioconductor packages it needs the
first time it's run (skipping ones you already have), so the blocks are
self-contained -- you can run Steps 3-6 in one sitting or come back to
any of them independently.

With per-replicate peaks called, build a consensus peak set across all samples and count reads in those regions — the input to differential binding analysis.

.. code:: r

   if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
   }
   BiocManager::install(c("DiffBind", "DESeq2"), update = FALSE, ask = FALSE)

   library(DiffBind)

   samples <- read.csv("samples.csv")
   dbObj <- dba(sampleSheet = samples)

   # overlap rate across samples, sanity check before counting
   dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID)

   # count reads in the consensus peak set
   dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)
   dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID)

``dba.count()`` re-centers and merges overlapping peaks into a consensus set, then counts reads per sample within those consensus intervals — directly analogous to building a counts matrix in RNA-seq, but over genomic intervals instead of genes. Re-running the PCA after counting (rather than just on peak overlap) shows whether IP signal, not just peak presence, separates the two cell lines.

Step 4: Differential binding analysis (DiffBind + DESeq2)
---------------------------------------------------------

Continuing in the same local RStudio project as Step 3 (``dbObj`` carries over).

.. code:: r

   if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
   }
   BiocManager::install(c("DiffBind", "DESeq2"), update = FALSE, ask = FALSE)

   library(DiffBind)

   dbObj <- dba.contrast(dbObj, categories = DBA_CONDITION, minMembers = 2)
   dbObj <- dba.analyze(dbObj, method = DBA_DESEQ2)

   dba.show(dbObj, bContrasts = TRUE)

   res <- dba.report(dbObj, method = DBA_DESEQ2, th = 1)  # th=1 keeps all sites; filter below
   res_df <- as.data.frame(res)

   # significant differential CTCF sites, GM12878 vs K562
   sig <- subset(res_df, FDR < 0.05 & abs(Fold) > 1)
   write.csv(sig, "diffbind_GM12878_vs_K562_sig.csv", row.names = FALSE)

   dba.plotMA(dbObj, method = DBA_DESEQ2)
   dba.plotVolcano(dbObj, method = DBA_DESEQ2)

Under the hood ``dba.analyze(method = DBA_DESEQ2)`` hands the consensus count matrix to DESeq2 exactly as you would for RNA-seq differential expression, modeling read counts per peak as negative binomial and testing condition effects — the difference is the “features” are genomic intervals, not genes. ``Fold`` is the log2 fold change (GM12878 vs K562, by the order of levels in ``Condition``); ``FDR`` is the Benjamini-Hochberg-adjusted p-value.

CTCF binding is largely invariant across cell types for a substantial fraction of sites (the “constitutive” CTCF sites that anchor TAD boundaries), so expect a sizeable shared peak set alongside a smaller set of cell-type- specific differential sites — these differential sites are usually the biologically interesting ones for cell-type-specific 3D genome structure.

Step 5: Peak annotation (ChIPseeker)
------------------------------------

Still in the same local RStudio project (``dbObj``, ``sig`` carry over from Steps 3-4).

Annotate both the full consensus peak set and the differential peak set to nearest genes and genomic features (promoter, exon, intron, intergenic).

.. code:: r

   if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
   }
   BiocManager::install(c(
       "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
       "org.Hs.eg.db", "GenomicRanges"
   ), update = FALSE, ask = FALSE)

   library(ChIPseeker)
   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
   library(org.Hs.eg.db)
   library(GenomicRanges)

   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

   # annotate the full consensus set
   all_peaks <- dba.peakset(dbObj, bRetrieve = TRUE)
   peakAnno <- annotatePeak(
       all_peaks,
       tssRegion = c(-3000, 3000),
       TxDb = txdb,
       annoDb = "org.Hs.eg.db"
   )

   plotAnnoPie(peakAnno)
   plotDistToTSS(peakAnno)

   anno_df <- as.data.frame(peakAnno)
   write.csv(anno_df, "all_peaks_annotated.csv", row.names = FALSE)

   # annotate just the significant differential sites from Step 4
   sig_gr <- makeGRangesFromDataFrame(sig, keep.extra.columns = TRUE)
   sigAnno <- annotatePeak(sig_gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
   sig_anno_df <- as.data.frame(sigAnno)
   write.csv(sig_anno_df, "diff_peaks_annotated.csv", row.names = FALSE)

``plotAnnoPie`` summarizes what fraction of peaks fall in promoters, exons, introns, or intergenic space — for CTCF expect a substantial intergenic and intronic fraction, since CTCF binds far more broadly than just at promoters, unlike e.g. H3K4me3. ``plotDistToTSS`` shows the distance distribution to the nearest transcription start site.

Step 6: Enrichment analysis (clusterProfiler)
---------------------------------------------

Still in the same local RStudio project (``sig_anno_df``, ``anno_df`` carry over from Step 5).

Take the genes nearest the differential CTCF sites (from ``diff_peaks_annotated.csv``) and test for GO/KEGG enrichment relative to all genes near any CTCF peak (the appropriate background for ChIP-seq, rather than the whole genome).

.. code:: r

   if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
   }
   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"), update = FALSE, ask = FALSE)

   library(clusterProfiler)
   library(org.Hs.eg.db)

   diff_genes <- unique(sig_anno_df$geneId)
   universe_genes <- unique(anno_df$geneId)   # all CTCF-bound genes as background

   ego <- enrichGO(
       gene = diff_genes,
       universe = universe_genes,
       OrgDb = org.Hs.eg.db,
       keyType = "ENTREZID",
       ont = "BP",
       pAdjustMethod = "BH",
       pvalueCutoff = 0.05,
       qvalueCutoff = 0.2
   )

   dotplot(ego, showCategory = 20)
   write.csv(as.data.frame(ego), "GO_BP_enrichment_diff_CTCF.csv", row.names = FALSE)

   ekegg <- enrichKEGG(
       gene = diff_genes,
       universe = universe_genes,
       organism = "hsa",
       pvalueCutoff = 0.05
   )

   dotplot(ekegg, showCategory = 20)
   write.csv(as.data.frame(ekegg), "KEGG_enrichment_diff_CTCF.csv", row.names = FALSE)

Using the CTCF-bound gene set as the statistical background (``universe``) rather than the whole genome is the important methodological point here: it controls for the fact that CTCF binding itself is not uniformly distributed across functional gene categories, so testing against the whole genome would inflate enrichment from CTCF’s baseline binding preferences rather than isolating what’s specific to the cell-type-differential sites.

.. note::

   For histone marks or factors with broader regulatory roles, also consider ``ReactomePA::enrichPathway()`` for pathway-level enrichment (installed above), and GSEA-style analysis (``gseGO``/``gseKEGG``) ranking all bound genes by fold change rather than thresholding to a "significant" gene list, which avoids an arbitrary cutoff.

Visualization
-------------

With results generated, bring everything together visually:

-  **IGV** — load the bigWig tracks from Step 1 alongside the ``narrowPeak`` files from Step 2 and the differential site BED file to inspect individual loci: https://igv.org/
-  **deeptools ``plotHeatmap``/``computeMatrix``** — signal heatmaps centered on peak summits or TSSs, useful for a publication-style summary figure:

.. code:: bash

   computeMatrix reference-point \
       -S bigwig/gm12878_ctcf_rep1.bw bigwig/k562_ctcf_rep1.bw \
       -R macs2/GM12878_rep1_summits.bed \
       --referencePoint center -a 2000 -b 2000 \
       -o matrix.gz -p 8

   plotHeatmap -m matrix.gz -o heatmap.png

This concludes the workshop. For questions or comments, treat this as a starting template and adapt the cell lines, factor, and contrasts to your own BAM files.
