ATAC-seq
========

Downstream analysis of bulk ATAC-seq (Assay for Transposase-Accessible
Chromatin) from paired-end Illumina sequencing of **human** samples, run on
locally installed HPC modules. The tutorial covers everything from quality
control through to differential accessibility analysis across multiple
conditions.

ATAC-seq uses the hyperactive Tn5 transposase to simultaneously fragment and
tag accessible (open) chromatin. The distribution of resulting fragments is a
direct readout of chromatin state: sub-nucleosomal fragments (under ~150 bp)
originate from nucleosome-free regions (NFR), whereas progressively longer
fragments reflect mono-, di-, and tri-nucleosomal wrapping. A healthy library
therefore shows a characteristic **nucleosomal ladder** in its fragment-size
distribution, and strong enrichment of signal at transcription start sites
(TSS).

This tutorial follows a four-condition experiment: a wildtype (WT) and three
distinct cell lines (N1, N2, N3), each compared back to WT.

.. note::

   This tutorial begins from **deduplicated, coordinate-sorted BAM files and
   their index files** (``.bai``). It assumes you have already aligned your
   paired-end reads to hg38, marked or removed duplicates, and coordinate
   sorted the result. Alignment, duplicate marking, and sorting are covered in
   the upstream preprocessing material and are not repeated here.

Pipeline overview
-----------------

.. figure:: ../_images/atac_workflow.png
   :alt: ATAC-seq downstream analysis workflow
   :align: center

   The downstream ATAC-seq workflow. The **main path** (green) runs from the
   input BAMs through filtering, peak calling, consensus construction,
   quantification, and differential accessibility. The two QC branches
   (orange) are informational: they validate library quality before committing
   to the analysis but do not alter the BAMs. Peak calling is run twice, on
   individual samples and on per-condition pools; the pooled peaks build the
   consensus set, while the individual peaks populate the DiffBind sample sheet
   (grey dashed auxiliary inputs).

Steps
-----

#. Assess library quality: fragment-size distribution and TSS enrichment.

#. Filter the BAMs: remove mitochondrial reads, apply a mapping-quality
   threshold, apply the Tn5 +4/-5 bp shift, and retain only NFR fragments.

#. Call peaks with MACS2 on each individual sample and on per-condition pooled
   BAMs.

#. Build a consensus (union) peak set from the pooled peaks and remove
   blacklisted regions.

#. Quantify reads per sample in the consensus peaks with ``featureCounts``.

#. Run differential accessibility in R using DiffBind with DESeq2 and
   background normalization, then extract and visualise results per contrast.

Environment
-----------

All commands assume the following modules are loaded. Exact names depend on
your HPC installation; substitute the versions available to you. Note that
peak calling here uses a different module tree (``gencore/2``) from the rest of
the command-line steps (``gencore/3``), so the module environment is switched
where needed.

::

   # BAM filtering, QC, and quantification
   module purge
   module load all gencore/3 deeptools samtools subread bedtools

   # Peak calling
   module purge
   module load all gencore/2 macs2

The differential-accessibility step (Step 6) is run locally in R / RStudio.
Required packages:

::

   install.packages(c("ggplot2", "ggrepel", "dplyr"))
   BiocManager::install(c("DiffBind", "DESeq2", "GenomicRanges"))

Inputs
------

Before running the pipeline you need:

- **Deduplicated, coordinate-sorted BAMs** and their ``.bai`` indexes, one per
  sample.

- **RPKM-normalized bigwig files**, one per sample, for the TSS-enrichment QC
  step. These are typically produced with deeptools ``bamCoverage``.

- A **reference genome build** identifier. This tutorial uses hg38.

- The **ENCODE hg38 blacklist** of problematic regions (downloaded below).

The example data set is 16 samples across four conditions:

.. list-table::
   :header-rows: 1
   :widths: 12 12 38 38

   * - Sample
     - Condition
     - BAM
     - Bigwig
   * - WT-1
     - WT
     - ``..._01_S16_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._01_S16_bamCoverage.bw``
   * - WT-4
     - WT
     - ``..._02_S15_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._02_S15_bamCoverage.bw``
   * - WT-5
     - WT
     - ``..._03_S14_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._03_S14_bamCoverage.bw``
   * - WT-6
     - WT
     - ``..._04_S13_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._04_S13_bamCoverage.bw``
   * - N1-3
     - N1
     - ``..._05_S12_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._05_S12_bamCoverage.bw``
   * - N1-4
     - N1
     - ``..._06_S11_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._06_S11_bamCoverage.bw``
   * - N1-5
     - N1
     - ``..._07_S10_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._07_S10_bamCoverage.bw``
   * - N1-6
     - N1
     - ``..._08_S9_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._08_S9_bamCoverage.bw``
   * - N2-1
     - N2
     - ``..._09_S8_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._09_S8_bamCoverage.bw``
   * - N2-2
     - N2
     - ``..._10_S7_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._10_S7_bamCoverage.bw``
   * - N2-3
     - N2
     - ``..._11_S6_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._11_S6_bamCoverage.bw``
   * - N2-4
     - N2
     - ``..._12_S5_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._12_S5_bamCoverage.bw``
   * - N3-1
     - N3
     - ``..._13_S4_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._13_S4_bamCoverage.bw``
   * - N3-3
     - N3
     - ``..._14_S3_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._14_S3_bamCoverage.bw``
   * - N3-4
     - N3
     - ``..._15_S2_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._15_S2_bamCoverage.bw``
   * - N3-5
     - N3
     - ``..._16_S1_mdup.withrg.csorted.cleaned.aligned.bam``
     - ``..._16_S1_bamCoverage.bw``

Quality control
---------------

Run QC on the raw deduplicated BAMs *before* any filtering, so you have a
baseline picture of library quality. Two metrics are essential for ATAC-seq:
the fragment-size distribution and the TSS enrichment score.

Generating a TSS BED file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TSS enrichment needs a BED file of transcription start sites. Derive one from
the GENCODE v45 annotation for hg38.

::

   # Download GENCODE v45 GTF
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz

   # Extract unique TSS positions, primary chromosomes only
   zcat gencode.v45.annotation.gtf.gz \
     | awk '$3 == "transcript"' \
     | awk 'BEGIN{OFS="\t"} {
         if ($7 == "+") print $1, $4-1, $4, ".", ".", $7;
         else print $1, $5-1, $5, ".", ".", $7
       }' \
     | grep -v "^#" \
     | grep -E "^chr([0-9]+|X|Y)\s" \
     | sort -k1,1 -k2,2n \
     | uniq \
     > hg38_TSS_primary.bed

.. note::

   The primary chromosome sequences in hg38 are identical across all GRCh38
   patch releases (p1 through p14); patches only add or revise alt contigs and
   scaffolds. A TSS file restricted to ``chr1``–``chr22``, ``chrX``, ``chrY``
   is therefore safe to use regardless of which patch your alignment reference
   was built from, as long as you restrict to the primary chromosomes. The
   ``grep -E "^chr([0-9]+|X|Y)\s"`` filter does this; the trailing ``\s``
   matches the tab after the chromosome name (using ``$`` would match
   end-of-line and silently return nothing).

Fragment-size distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use deeptools ``bamPEFragmentSize`` across all samples at once.

::

   mkdir -p QC

   bamPEFragmentSize \
     -b WT-1.bam WT-4.bam WT-5.bam WT-6.bam \
        N1-3.bam N1-4.bam N1-5.bam N1-6.bam \
        N2-1.bam N2-2.bam N2-3.bam N2-4.bam \
        N3-1.bam N3-3.bam N3-4.bam N3-5.bam \
     --samplesLabel WT-1 WT-4 WT-5 WT-6 N1-3 N1-4 N1-5 N1-6 \
                    N2-1 N2-2 N2-3 N2-4 N3-1 N3-3 N3-4 N3-5 \
     --histogram QC/fragment_size_distribution.png \
     --outRawFragmentLengths QC/fragment_sizes_raw.tab \
     -p 28

A healthy ATAC library shows a dominant NFR peak below 150 bp, with smaller
mono- and di-nucleosomal shoulders at roughly 200 bp intervals, and tight
overlap of all sample curves.

TSS enrichment
~~~~~~~~~~~~~~~

Compute a signal matrix centred on every TSS (±2 kb) from the bigwig files,
then plot a profile and a heatmap. TSS enrichment above 7 is acceptable for
human ATAC-seq, and above 10 is good.

::

   computeMatrix reference-point \
     --referencePoint TSS \
     -b 2000 -a 2000 \
     -R hg38_TSS_primary.bed \
     -S WT-1.bw WT-4.bw WT-5.bw WT-6.bw \
        N1-3.bw N1-4.bw N1-5.bw N1-6.bw \
        N2-1.bw N2-2.bw N2-3.bw N2-4.bw \
        N3-1.bw N3-3.bw N3-4.bw N3-5.bw \
     --skipZeros \
     -o QC/TSS_matrix.gz \
     -p 28

   plotProfile \
     -m QC/TSS_matrix.gz \
     -o QC/TSS_enrichment_profile.png \
     --perGroup \
     --plotTitle "TSS Enrichment - All Samples"

   plotHeatmap \
     -m QC/TSS_matrix.gz \
     -o QC/TSS_enrichment_heatmap.png \
     --plotTitle "TSS Enrichment Heatmap" \
     --colorMap RdYlBu_r

What to look for in the QC plots:

- **Fragment size**: a strong sub-150 bp NFR peak, nucleosomal shoulders, and
  tight overlap across samples. A divergent curve flags a depth or quality
  outlier.

- **TSS profile / heatmap**: a sharp, symmetric peak centred exactly on the
  TSS in every sample, consistent peak height within a condition, and a clear
  enrichment band in each heatmap column.

.. warning::

   TSS profiles computed from whole-fragment bigwigs (as here) look broader
   than those from NFR-filtered bigwigs. This is expected at the QC stage and
   is **not** a reason to fail a sample. Profiles sharpen considerably after
   the filtering step below.

BAM filtering
-------------

Several classes of reads in the raw BAM are uninformative or harmful for
ATAC-seq and must be removed before peak calling.

- **Mitochondrial reads (chrM).** The mitochondrial genome is essentially
  naked, highly accessible DNA, and a single cell contains many copies. Without
  removal, chrM reads can dominate an ATAC library (often 30–80 % of reads) and
  dilute nuclear signal.

- **Low-MAPQ reads.** Ambiguously placed multi-mappers are dropped by requiring
  ``MAPQ >= 30``.

- **Non-NFR fragments.** Retaining only sub-nucleosomal fragments (10–150 bp)
  focuses the analysis on genuine open-chromatin signal rather than nucleosome
  positioning.

- **Tn5 insertion offset.** Tn5 inserts as a dimer that leaves a 9 bp
  duplication, so the true insertion point is offset +4 bp on the ``+`` strand
  and -5 bp on the ``-`` strand from the reported alignment start. This shift
  must be corrected for accurate footprinting and motif analysis.

deeptools ``alignmentSieve`` applies all of these in a single pass. The
``--ATACshift`` flag performs the Tn5 correction.

::

   mkdir -p filtered

   alignmentSieve \
     -b WT-1.bam \
     -o filtered/WT-1_filtered.bam \
     --ATACshift \
     --minFragmentLength 10 \
     --maxFragmentLength 150 \
     --minMappingQuality 30 \
     --ignoreDuplicates \
     --blackListFileName hg38_blacklist.bed \
     -p 28

   samtools sort  -@ 28 -o filtered/WT-1_filtered_csorted.bam filtered/WT-1_filtered.bam
   samtools index -@ 28 filtered/WT-1_filtered_csorted.bam
   rm filtered/WT-1_filtered.bam

Download the ENCODE hg38 blacklist (the unified exclusion list recommended by
the Kundaje lab) before running:

::

   wget -L https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
   gunzip ENCFF356LFX.bed.gz
   mv ENCFF356LFX.bed hg38_blacklist.bed

On a cluster, run one job per sample. A generator script that writes the 16
SLURM scripts keeps submission tidy:

::

   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH -A gencore
   #SBATCH --ntasks=1
   #SBATCH --partition=gencore
   #SBATCH --cpus-per-task=28
   #SBATCH --time=12000:00:00
   #SBATCH --mem 420000

   module purge
   module load all gencore/3 deeptools samtools

   cd /path/to/bams
   mkdir -p filtered

   alignmentSieve \
     -b <SAMPLE>.bam \
     -o filtered/<LABEL>_filtered.bam \
     --ATACshift --minFragmentLength 10 --maxFragmentLength 150 \
     --minMappingQuality 30 --ignoreDuplicates \
     --blackListFileName hg38_blacklist.bed -p 28

   samtools sort  -@ 28 -o filtered/<LABEL>_filtered_csorted.bam filtered/<LABEL>_filtered.bam
   samtools index -@ 28 filtered/<LABEL>_filtered_csorted.bam
   rm filtered/<LABEL>_filtered.bam

.. warning::

   Keep the ``#SBATCH`` header minimal. On some cluster configurations adding
   ``--job-name``, ``--output``, or ``--error`` directives causes the job to
   exit silently within a couple of seconds with no logs written, because
   SLURM tries to open the log paths before the script body (which creates
   those directories) has run. If you need custom log paths, create the log
   directory *before* submitting, or use absolute paths that already exist.

Peak calling
------------

Peaks are called with MACS2 in paired-end mode (``-f BAMPE``). Run peak calling
twice: on each individual sample and on per-condition pooled BAMs. The
individual peaks populate the DiffBind sample sheet in Step 6; the pooled peaks
build the consensus set and give more sensitive detection of lower-occupancy
sites.

Individual samples
~~~~~~~~~~~~~~~~~~~

::

   macs2 callpeak \
     -t filtered/WT-1_filtered_csorted.bam \
     -f BAMPE \
     -n WT-1 \
     --outdir peaks/individual/WT-1 \
     --nomodel --shift -100 --extsize 200 \
     --nolambda --keep-dup all \
     -q 0.05 --call-summits \
     -g hs \
     2> peaks/individual/WT-1/WT-1_macs2.log

The parameter choices:

- ``--nomodel --shift -100 --extsize 200`` — the standard ATAC-seq setting.
  Reads are extended around the Tn5 cut site rather than MACS2 building a
  fragment model.

- ``--nolambda`` — disables local background (lambda) estimation. ATAC-seq has
  no input control, and local lambda can bias peak calls, so it is turned off.

- ``--keep-dup all`` — duplicates were already removed upstream; do not let
  MACS2 remove reads again.

- ``--call-summits`` — reports the precise summit within each peak, which is
  needed for downstream motif analysis.

Pooled conditions
~~~~~~~~~~~~~~~~~~

Merge the per-replicate filtered BAMs for each condition, then call peaks on
the pool. Remember to switch module trees between ``samtools`` and ``macs2``.

::

   # Merge replicates (example: WT)
   module purge && module load all gencore/3 samtools

   samtools merge -@ 28 -f filtered/WT_pooled.bam \
       filtered/WT-1_filtered_csorted.bam \
       filtered/WT-4_filtered_csorted.bam \
       filtered/WT-5_filtered_csorted.bam \
       filtered/WT-6_filtered_csorted.bam
   samtools index -@ 28 filtered/WT_pooled.bam

   # Call peaks on the pool
   module purge && module load all gencore/2 macs2

   macs2 callpeak \
     -t filtered/WT_pooled.bam \
     -f BAMPE -n WT_pooled \
     --outdir peaks/pooled/WT \
     --nomodel --shift -100 --extsize 200 \
     --nolambda --keep-dup all -q 0.05 --call-summits \
     -g hs \
     2> peaks/pooled/WT/WT_pooled_macs2.log

Inspecting peak counts
~~~~~~~~~~~~~~~~~~~~~~~~

Before going further, tabulate peak counts per sample. Consistent counts
across replicates are reassuring; a sample with far fewer peaks than its
siblings should be investigated (see `Handling outlier samples`_).

::

   for dir in peaks/individual/*/; do
     sample=$(basename $dir)
     count=$(wc -l < ${dir}/${sample}_peaks.narrowPeak)
     echo "${sample}: ${count} peaks"
   done

Consensus peak set
------------------

A consensus (union) peak set is a single coordinate system shared by all
samples. It is built from the pooled condition peaks, which capture
condition-specific accessible regions that any one sample might miss, and used
as the feature space for quantification and differential testing.

::

   mkdir -p peaks/consensus
   module purge && module load all gencore/3 bedtools

   # Concatenate pooled peaks, sort, merge overlaps, keep primary chromosomes
   cat peaks/pooled/WT/WT_pooled_peaks.narrowPeak \
       peaks/pooled/N1/N1_pooled_peaks.narrowPeak \
       peaks/pooled/N2/N2_pooled_peaks.narrowPeak \
       peaks/pooled/N3/N3_pooled_peaks.narrowPeak \
     | sort -k1,1 -k2,2n \
     | bedtools merge -i stdin \
     | grep -E "^chr([0-9]+|X|Y)\s" \
     > peaks/consensus/consensus_peaks_raw.bed

   # Remove blacklisted regions entirely
   bedtools subtract \
     -a peaks/consensus/consensus_peaks_raw.bed \
     -b hg38_blacklist.bed \
     -A \
     > peaks/consensus/consensus_peaks_final.bed

   wc -l peaks/consensus/consensus_peaks_final.bed

.. note::

   The ``-A`` flag in ``bedtools subtract`` drops any consensus peak that
   overlaps a blacklist region *at all*, rather than trimming only the
   overlapping bases. A peak that straddles a blacklisted region is unreliable
   along its whole length, so removing it entirely is the correct behaviour.

Quantification
--------------

Count the read pairs falling in each consensus peak for every sample. This
produces the peaks × samples count matrix that feeds differential testing.
``featureCounts`` requires the peaks in SAF format.

::

   module purge && module load all gencore/3 subread

   # Convert consensus BED to SAF
   awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
        {print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}' \
     peaks/consensus/consensus_peaks_final.bed \
     > peaks/consensus/consensus_peaks_final.saf

   # Count reads across all samples
   featureCounts \
     -F SAF \
     -a peaks/consensus/consensus_peaks_final.saf \
     -o peaks/consensus/consensus_counts.txt \
     -p --countReadPairs \
     -T 28 -f \
     filtered/WT-1_filtered_csorted.bam \
     filtered/WT-4_filtered_csorted.bam \
     filtered/WT-5_filtered_csorted.bam \
     filtered/WT-6_filtered_csorted.bam \
     filtered/N1-3_filtered_csorted.bam \
     filtered/N1-4_filtered_csorted.bam \
     filtered/N1-5_filtered_csorted.bam \
     filtered/N2-1_filtered_csorted.bam \
     filtered/N2-2_filtered_csorted.bam \
     filtered/N2-3_filtered_csorted.bam \
     filtered/N2-4_filtered_csorted.bam \
     filtered/N3-1_filtered_csorted.bam \
     filtered/N3-3_filtered_csorted.bam \
     filtered/N3-4_filtered_csorted.bam \
     filtered/N3-5_filtered_csorted.bam \
     2> peaks/consensus/featureCounts.log

The per-sample "Successfully assigned" rate in ``featureCounts.log`` is the
equivalent of a FRiP (fraction of reads in peaks) score. Rates of 70 % or above
across all samples indicate the consensus peak set captures the signal well.

Differential accessibility
---------------------------

Why DiffBind rather than plain DESeq2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DESeq2's default normalization (median-of-ratios) assumes that most features
are unchanged between conditions — a reasonable assumption for RNA-seq, where
the bulk of genes are stably expressed. ATAC-seq comparisons between
genetically distinct cell lines can violate this badly: **global chromatin
accessibility may genuinely differ** between conditions, and normalizing on the
peaks themselves can then misattribute a global shift to thousands of
individually "significant" peaks, inflating both fold changes and significance.

DiffBind offers **background normalization** as the remedy. Instead of deriving
size factors from the peaks, it tiles the genome into large bins (default
15 kb) and estimates library-level scaling from read counts in those bins.
Because the bins sample open and closed chromatin impartially, the resulting
factors reflect sequencing-depth differences rather than accessibility shifts.
This is conceptually similar to spike-in normalization and is the recommended
approach whenever conditions may differ in global accessibility.

Sample sheet and counting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

DiffBind is driven by a sample sheet listing, for each sample, its condition,
replicate, BAM, and peak file. The BAMs must be accessible because DiffBind
reads them both to count reads in the consensus peaks and to compute the
background bins.

.. code-block:: r

   library(DiffBind)
   library(GenomicRanges)

   sample_sheet <- data.frame(
     SampleID   = c("WT-1","WT-4","WT-5","WT-6",
                    "N1-3","N1-4","N1-5",
                    "N2-1","N2-2","N2-3","N2-4",
                    "N3-1","N3-3","N3-4","N3-5"),
     Condition  = c("WT","WT","WT","WT",
                    "N1","N1","N1",
                    "N2","N2","N2","N2",
                    "N3","N3","N3","N3"),
     Replicate  = c(1,2,3,4, 1,2,3, 1,2,3,4, 1,2,3,4),
     bamReads   = paste0("bams/",
                    c("WT-1","WT-4","WT-5","WT-6",
                      "N1-3","N1-4","N1-5",
                      "N2-1","N2-2","N2-3","N2-4",
                      "N3-1","N3-3","N3-4","N3-5"),
                    "_filtered_csorted.bam"),
     Peaks      = paste0("np/",
                    c("WT-1","WT-4","WT-5","WT-6",
                      "N1-3","N1-4","N1-5",
                      "N2-1","N2-2","N2-3","N2-4",
                      "N3-1","N3-3","N3-4","N3-5"),
                    "_peaks.narrowPeak"),
     PeakCaller = "narrow",
     stringsAsFactors = FALSE
   )

   # Load the consensus peaks as the counting interval set
   consensus_bed <- read.table("consensus_peaks_final.bed",
                               header = FALSE, sep = "\t",
                               col.names = c("Chr","Start","End"))
   consensus_gr <- GRanges(seqnames = consensus_bed$Chr,
                           ranges   = IRanges(consensus_bed$Start,
                                              consensus_bed$End))

   dba_obj <- dba(sampleSheet = sample_sheet)
   dba_obj <- dba.count(dba_obj,
                        peaks                 = consensus_gr,
                        score                 = DBA_SCORE_READS,
                        fragmentSize          = 0,
                        bUseSummarizeOverlaps = TRUE)

Background normalization, contrasts, and analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: r

   # The key step: size factors from 15 kb background bins
   dba_obj <- dba.normalize(dba_obj,
                            method     = DBA_DESEQ2,
                            normalize  = DBA_NORM_BACKGROUND,
                            background = 15000)

   # WT as the reference level for all contrasts
   dba_obj <- dba.contrast(dba_obj,
                           reorderMeta = list(Condition = "WT"))

   dba_obj <- dba.analyze(dba_obj,
                          method     = DBA_DESEQ2,
                          bBlacklist = FALSE,   # already filtered upstream
                          bGreylist  = FALSE)

Extracting and saving results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Retrieve the full result table for each contrast (N1 vs WT, N2 vs WT, N3 vs
WT), then apply significance thresholds. It is good practice to save both a
**standard** set (FDR < 0.05, |log2FC| > 1) and a **strict** set that adds a
minimum-concentration filter to suppress low-count peaks.

.. code-block:: r

   comparisons <- list(
     list(name = "N1_vs_WT", contrast = 1),
     list(name = "N2_vs_WT", contrast = 2),
     list(name = "N3_vs_WT", contrast = 3)
   )

   for (comp in comparisons) {
     res <- dba.report(dba_obj,
                       contrast    = comp$contrast,
                       method      = DBA_DESEQ2,
                       th          = 1,       # return all peaks
                       bUsePval    = FALSE,
                       fold        = 0,
                       bNormalized = TRUE)
     df <- as.data.frame(res)

     # Standard significance
     sig <- subset(df, FDR < 0.05 & abs(Fold) > 1)
     write.csv(sig, paste0("results/", comp$name, "_DA_significant.csv"),
               row.names = FALSE)

     # Strict: additionally require Conc >= log2(20)
     sig_strict <- subset(sig, Conc >= log2(20))
     write.csv(sig_strict, paste0("results/strict/", comp$name,
               "_DA_significant.csv"), row.names = FALSE)
   }

Sample-level diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~

A PCA and a correlation heatmap of the normalized counts are the quickest check
that replicates cluster by condition before trusting any contrast.

.. code-block:: r

   dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
   dba.plotHeatmap(dba_obj)

Well-separated condition clusters on PC1/PC2 and tight within-condition
grouping indicate that the chromatin landscapes differ substantially between
conditions and that the replicates are reproducible.

Handling outlier samples
------------------------

Always reconcile peak counts and library depth across samples before
finalising. A replicate with a dramatically lower peak count than its group is
usually undersequenced rather than biologically distinct, and should be
excluded.

::

   # Compare library depth of a suspected outlier against its group
   samtools flagstat filtered/N1-6_filtered_csorted.bam

.. warning::

   In the example data set, sample **N1-6** carried only ~12.6 M mapped reads
   after filtering, versus ~110–118 M for the other N1 replicates — roughly a
   ten-fold deficit. This produced only ~145 k called peaks against ~330 k for
   N1-3, N1-4, and N1-5. N1-6 was excluded from all downstream analysis, and
   the N1 pooled peak call was re-run using only N1-3, N1-4, and N1-5.

After excluding N1-6, the final layout used for differential analysis is:

.. list-table::
   :header-rows: 1
   :widths: 20 55 10

   * - Condition
     - Replicates
     - n
   * - WT
     - WT-1, WT-4, WT-5, WT-6
     - 4
   * - N1
     - N1-3, N1-4, N1-5
     - 3
   * - N2
     - N2-1, N2-2, N2-3, N2-4
     - 4
   * - N3
     - N3-1, N3-3, N3-4, N3-5
     - 4

.. note::

   Three replicates per condition is the minimum for a well-powered
   differential test with DESeq2 / DiffBind, so dropping N1-6 still leaves N1
   adequately powered. If excluding an outlier would drop a condition below
   three replicates, consider resequencing rather than proceeding.
