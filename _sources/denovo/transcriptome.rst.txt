.. _denovo-transcriptome:

###########################################
De Novo Transcriptome Assembly & Annotation
###########################################

De novo assembly, quality assessment, functional annotation, and expression
quantification of a **transcriptome** from paired-end short reads (Illumina),
run on locally installed HPC modules.

This workshop uses publicly available RNA-Seq data from the **non-model
prokaryote** *Staphylococcus aureus* grown under three carbon sources
(glucose, fructose, pyruvate; two replicates each). No reference genome is
assumed: the transcriptome is built directly from the reads with Trinity,
assessed, annotated against homology databases, and then used as a reference
for abundance estimation.

Scope is limited to a **single, guided end-to-end pass** on this dataset.
It does not attempt to compare assemblers, tune parameters per species, or
cover every downstream analysis. The aim is a reference workflow you can
adapt to your own data.

.. note::

   The analysis runs on the NYUAD High Performance Computing (HPC) cluster,
   but every step works on any standalone machine (server, laptop, desktop)
   provided the software stack below is installed. We recommend ``conda`` for
   installing and maintaining bioinformatics software.


Pipeline overview
=================

.. figure:: /_images/denovo-workflow.png
   :alt: De novo transcriptome assembly and annotation pipeline
   :align: center

   The workflow. Reads are assembled once with Trinity into a reference
   transcriptome; that assembly is then assessed (read mapping, Nx/ExN50,
   BUSCO, QUAST), annotated (TransDecoder ORFs + Trinotate homology search),
   and finally quantified per sample (RSEM) for differential expression
   (DESeq2). The assembly is the hub: every downstream step consumes it.

.. note::

   Drop your own workflow image at ``_images/denovo-workflow.png`` (or change
   the path in the ``figure`` directive). If you have ``sphinx-design`` or
   Graphviz enabled, this is also a natural place for a generated diagram.


Steps
=====

#. Assemble the reads into a reference transcriptome with **Trinity**.
#. Assess the assembly: read mapping rate, contig **Nx / ExN50**, **BUSCO**
   completeness, and **QUAST** mis-assembly metrics.
#. Predict ORFs with **TransDecoder** and annotate with **Trinotate**
   (BLAST/DIAMOND vs. SwissProt, PFAM via HMMER, EggNOG, TmHMM, SignalP).
#. Quantify expression per sample with **RSEM**, build expression matrices,
   and run differential expression with **DESeq2**.
#. Visualize alignments and annotations in **IGV** and **TrinotateWeb**.


Environment
===========

All commands assume the following modules are loaded. Exact names and versions
depend on your HPC installation; substitute the versions available to you.

.. code-block:: bash

   module purge
   module load all
   module load gencore/2
   module load trinity/2.15.1
   module load salmon/0.14.2
   module load numpy/1.26.4
   module load bowtie2/2.4.5
   module load samtools/1.10
   module load hisat2/2.2.1

.. note::

   Different stages load different stacks (BUSCO, QUAST, and Trinotate each
   need their own environment). Those are given inline in the relevant
   sections below rather than all up front, because they conflict if loaded
   together. Parameters can change between tool versions, so verify flags
   against the version you actually have.

The tools used across the workshop, with upstream documentation:

- Trinity — `<https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_
- Trinotate — `<https://github.com/Trinotate/Trinotate/wiki>`_
- Kallisto — `<https://github.com/pachterlab/kallisto>`_
- Salmon — `<https://github.com/COMBINE-lab/salmon>`_
- Augustus — `<https://bioinf.uni-greifswald.de/augustus/>`_
- rnaQUAST — `<https://github.com/ablab/rnaquast>`_
- BUSCO — `<https://busco.ezlab.org/>`_
- SAMtools — `<https://www.htslib.org/>`_
- IGV — `<https://www.igv.org/>`_


Inputs
======

Before running the pipeline you need:

- **Paired-end RNA-Seq reads** in FASTQ format. This dataset has 6 samples
  (12 files), 3 conditions × 2 replicates. Public accessions:
  ``SRR28281136``–``SRR28281141`` from the SRA
  (`<https://www.ncbi.nlm.nih.gov/sra>`_).

- A **samples / conditions file** (``cond.txt``), tab-delimited, mapping each
  replicate to its read pair (format given in
  :ref:`denovo-running-trinity`).

- A **pre-generated assembly** (``trinity_assembly.Trinity.fasta``) shipped
  with the data. Assembly is CPU/memory/time intensive, so for the workshop
  it is provided ready-made; the commands to regenerate it are shown anyway.

- A **pre-initialized Trinotate database** (``myTrinotate.sqlite``), also
  shipped with the data so you do not have to download and build the homology
  databases during the session.

All other files in the pipeline are produced by the commands below.

File formats you will encounter (assumed familiar): FASTA, FASTQ, SAM, BAM.


Getting the data
================

We use the NYUAD HPC cluster. Connect, move to your scratch space, and copy
the workshop data.

Connecting from macOS / Linux
-----------------------------

#. Open **Terminal** and run ``ssh NetID@jubail.abudhabi.nyu.edu``. Enter your
   NYU password when prompted.
#. Move to your scratch directory: ``cd $SCRATCH``.
#. Create and enter a working directory:
   ``mkdir de_novo_transcriptome && cd de_novo_transcriptome``.
#. Copy the data: ``cp -r /scratch/gencore/nd48/workshop_de_novo_trans/data .``

Connecting from Windows
-----------------------

#. Open **PuTTY** and set **Host name** = ``jubail.abudhabi.nyu.edu``,
   **Port** = ``22``, then click **Open**.
#. Enter your NetID and password when prompted.
#. Move to your scratch directory: ``cd $SCRATCH``.
#. Create and enter a working directory:
   ``mkdir de_novo_transcriptome && cd de_novo_transcriptome``.
#. Copy the data: ``cp -r /scratch/gencore/nd48/workshop_de_novo_trans/data .``


The data
========

The datasets are publicly available sequencing reads in FASTQ format,
downloadable from the SRA under accessions ``SRR28281136``, ``SRR28281137``,
``SRR28281138``, ``SRR28281139``, ``SRR28281140``, and ``SRR28281141``.

There are 6 samples organized as paired-end sequences (12 files), spanning
3 conditions with 2 replicates each. The organism is *Staphylococcus aureus*,
and the conditions represent bacterial growth under three carbon sources:
glucose, fructose, and pyruvate.


The scripts
===========

The data folder ships with several shell scripts:

- ``assembly.sh`` — assembles the transcriptome with Trinity.
- ``qc.sh`` — assesses the assembly using read alignment, BUSCO, and QUAST.
- ``annotate.sh`` — runs TransDecoder and Trinotate for annotation.
- ``quantify.sh`` — quantifies expression and runs DGE (RSEM + DESeq2).
- ``de_novo.sh`` — combines all of the above into a single submission.

.. note::

   To run the complete analysis as one job, submit ``de_novo.sh``. The
   individual scripts exist so each stage can be run and inspected on its own.


Assembly with Trinity
======================

Most sequencing analyses begin with quality checking and trimming to remove
low-quality bases and adapter contamination. That step is not covered here;
see the quality-checking section of the Variant Detection & Annotation
workshop
(`<https://github.com/nizardrou/Variant-Detection-and-Annotation-workshop>`_).

We use **Trinity** for assembly. There are alternatives, but Trinity produces
reliable assemblies across organism sizes and ships a rich set of helper
scripts for post-processing. The golden rule still applies: no single tool is
best for every dataset, so for real projects run several assemblers and
compare.

.. note::

   The data folder already contains a pre-built assembly,
   ``trinity_assembly.Trinity.fasta``. Assembly takes 2+ hours and is
   resource-heavy, so the workshop skips the run itself and uses the
   pre-built FASTA for all downstream steps. The Trinity commands are shown
   (and are commented out in ``de_novo.sh``) so you can reproduce them later.

SLURM arguments
---------------

The first lines of the script are SLURM directives telling the scheduler where
and how to run the job. They are ignored if you are not submitting with
``sbatch``, so you can leave them in place when running locally.

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=28
   #SBATCH --time=04:00:00
   #SBATCH --mem 100000

De novo assembly is resource intensive, so we request 100 GB of memory and 28
CPUs. Larger datasets and larger transcriptomes need more. A useful rule of
thumb is at least **1 GB of memory per 1 million reads**; for very large read
sets, run digital normalization before assembling
(`<https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization>`_).

Loading the software stack
--------------------------

These ``module load`` lines make the tools available. They are specific to the
NYUAD HPC; adapt them to your environment.

.. code-block:: bash

   module purge
   module load all
   module load gencore/2
   module load trinity/2.15.1
   module load salmon/0.14.2
   module load numpy/1.26.4
   module load bowtie2/2.4.5
   module load samtools/1.10
   module load hisat2/2.2.1

.. _denovo-running-trinity:

Running the assembly
--------------------

The command below runs the Trinity assembly: input reads are FASTQ
(``--seqType fq``), memory capped at 100 GB, 28 CPUs, output prefixed
``trinity_assembly``, and samples described in the tab-delimited ``cond.txt``.

.. code-block:: bash

   Trinity \
       --seqType fq \
       --samples_file cond.txt \
       --CPU 28 \
       --output trinity_assembly \
       --max_memory 100G

``cond.txt`` is tab-delimited, one row per replicate, in the column order:

.. code-block:: text

   CONDITION  CONDITION_REPLICATE_NUMBER  READ1_FASTQ_FILE_NAME  READ2_FASTQ_FILE_NAME

For this dataset it looks like:

.. code-block:: text

   fructose   fructose_rep1   fru_1_1.fastq   fru_1_2.fastq
   fructose   fructose_rep2   fru_2_1.fastq   fru_2_2.fastq
   glucose    glucose_rep1    glu_1_1.fastq   glu_1_2.fastq
   glucose    glucose_rep2    glu_2_1.fastq   glu_2_2.fastq
   pyruvate   pyruvate_rep1   pyr_1_1.fastq   pyr_1_2.fastq
   pyruvate   pyruvate_rep2   pyr_2_1.fastq   pyr_2_2.fastq

Specifying samples this way pays off later: the same ``cond.txt`` is reused
for per-sample abundance estimation after the assembly is built.

Submitting the assembly
-----------------------

Submit the script to SLURM:

.. code-block:: bash

   sbatch assembly.sh

Or run it directly on the command line:

.. code-block:: bash

   chmod 755 assembly.sh
   ./assembly.sh


Assessing the assembly quality
==============================

Before annotating or quantifying, confirm the assembly is sound. We use three
complementary angles: how well the reads map back, length/expression
statistics, and gene-content completeness. The QC steps live in ``qc.sh``.

Assess using read mapping
-------------------------

Ideally the reads that produced the assembly should map back to it at a high
rate. A low rate can signal insufficient coverage, contamination, poor quality,
a complex transcriptome the assembler struggled with, or partially spliced
genes pointing to an extraction-protocol issue.

Use a splice-aware aligner — here **HISAT2** with **SAMtools**.

#. Index the transcriptome:

   .. code-block:: bash

      hisat2-build -p 28 trinity_assembly.Trinity.fasta trinity_assembly

#. Align the reads:

   .. code-block:: bash

      hisat2 \
          -x trinity_assembly \
          -p 28 \
          -1 fru_1_1.fastq,fru_2_1.fastq,glu_1_1.fastq,glu_2_1.fastq,pyr_1_1.fastq,pyr_2_1.fastq \
          -2 fru_1_2.fastq,fru_2_2.fastq,glu_1_2.fastq,glu_2_2.fastq,pyr_1_2.fastq,pyr_2_2.fastq \
          -S aln.sam

#. Convert to BAM, coordinate-sort, index, and gather alignment stats:

   .. code-block:: bash

      samtools view -@ 28 -b -o aln.bam aln.sam
      samtools sort -@ 28 -o aln.sorted.bam aln.bam
      samtools index -@ 28 aln.sorted.bam
      samtools flagstat aln.sorted.bam > stats.txt

Look at ``stats.txt`` for the overall alignment rate. Aim for **>80%**; this
assembly reaches ~93%, meaning most reads are represented in the assembly.

.. code-block:: text

   61353459 reads; of these:
     61353459 (100.00%) were paired; of these:
       10998445 (17.93%) aligned concordantly 0 times
       36304307 (59.17%) aligned concordantly exactly 1 time
       14050707 (22.90%) aligned concordantly >1 times
       ----
       10998445 pairs aligned concordantly 0 times; of these:
         1097373 (9.98%) aligned discordantly 1 time
       ----
       9901072 pairs aligned 0 times concordantly or discordantly; of these:
         19802144 mates make up the pairs; of these:
           8727873 (44.08%) aligned 0 times
           7094790 (35.83%) aligned exactly 1 time
           3979481 (20.10%) aligned >1 times
   92.89% overall alignment rate

Dedicated alignment-QC tools (e.g. Qualimap) give a more detailed breakdown,
but this is enough for an initial check.

Contig Nx and ExN50 statistics
------------------------------

Trinity ships helper scripts for length statistics. **N50** is the contig
length at which 50% of the total assembly length is contained in contigs of
that length or longer — a contiguity indicator, where higher *may* mean a
better assembly. It says nothing about correctness, so never use it alone.

.. code-block:: bash

   TrinityStats.pl trinity_assembly.Trinity.fasta > assembly_stats.txt

.. code-block:: text

   ################################
   ## Counts of transcripts, etc.
   ################################
   Total trinity 'genes':	3040
   Total trinity transcripts:	3794
   Percent GC: 34.52

   ########################################
   Stats based on ALL transcript contigs:
   ########################################

   	Contig N10: 14232
   	Contig N20: 10849
   	Contig N30: 8233
   	Contig N40: 5784
   	Contig N50: 4271

   	Median contig length: 310
   	Average contig: 1122.70
   	Total assembled bases: 4259519

   #####################################################
   ## Stats based on ONLY LONGEST ISOFORM per 'GENE':
   #####################################################

   	Contig N10: 14232
   	Contig N20: 10599
   	Contig N30: 7935
   	Contig N40: 5575
   	Contig N50: 4175

   	Median contig length: 302
   	Average contig: 1097.71
   	Total assembled bases: 3337029

We read the gene/transcript counts and total size alongside N50. At ~3K genes
we are in the expected ballpark for this organism. The ~4.2 Mb of assembled
bases is larger than expected, which may reflect background genomic DNA
contamination or captured plasmids.

A more transcriptome-appropriate statistic is the **Ex90N50** and **Ex90 gene
count**: N50 computed over only the most highly expressed genes that account
for 90% of normalized expression. See the Trinity wiki
(`<https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats>`_).

Assessing completeness with BUSCO
---------------------------------

**BUSCO** (Benchmarking Universal Single-Copy Orthologs) assesses an assembly
by searching for conserved single-copy core genes. It works on eukaryotes and
prokaryotes; higher scores mean a more complete assembly. Unlike the read- and
length-based checks above, BUSCO evaluates assembly *content*.

Load its environment:

.. code-block:: bash

   module purge
   module load all
   module load gencore/2
   export PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/bin:$PATH"
   export AUGUSTUS_CONFIG_PATH="/scratch/gencore/.eb/2.0/software/augustus/3.4.0/config/"
   module load busco/5.2.0

Run BUSCO in transcriptome mode with automatic prokaryote lineage selection:

.. code-block:: bash

   busco \
       -i trinity_assembly.Trinity.fasta \
       --out BUSCO \
       -m tran \
       --auto-lineage-prok \
       -c 28

Results are written to the ``BUSCO`` folder.

Detecting mis-assemblies with QUAST
-----------------------------------

**QUAST** (Quality Assessment Tool for Genome Assemblies) adds mis-assembly
detection on top of the length statistics. Mis-assemblies are contigs stitched
together that should not be — often around repeats or low-complexity regions,
or where coverage is too low. They are found by mapping reads back and looking
for regions of zero or near-zero coverage: an assembled region with no read
support should not exist. QUAST can also wrap ab initio gene finding and BUSCO,
so it makes sense to run them together.

QUAST takes a single read pair, not a comma-separated list, so concatenate all
read 1s and all read 2s first:

.. code-block:: bash

   cat *_1.fastq > read1.fastq
   cat *_2.fastq > read2.fastq

Load the QUAST environment:

.. code-block:: bash

   module purge
   module load gencore
   module load Miniconda3/4.7.10
   source activate /scratch/gencore/conda3/envs/quast_env
   export PATH="/scratch/gencore/software/genemarks_t:$PATH"

Run QUAST:

.. code-block:: bash

   quast.py \
       -t 28 \
       -f \
       -b \
       -o QUAST \
       --rna-finding \
       --bam aln.sorted.bam \
       -1 read1.fastq \
       -2 read2.fastq \
       trinity_assembly.Trinity.fasta

Flags:

- ``-t 28`` — number of CPUs.
- ``-f`` — enable gene finding (GeneMark; the prokaryote default we want).
- ``-b`` — run BUSCO.
- ``-o QUAST`` — output folder.
- ``--rna-finding`` — detect ribosomal RNA genes.
- ``--bam aln.sorted.bam`` — the sorted BAM from the read-mapping step, used
  for structural-variation detection and coverage histograms.
- ``-1 read1.fastq -2 read2.fastq`` — the concatenated reads, aligned to the
  assembly to detect mis-assemblies and coverage.
- ``trinity_assembly.Trinity.fasta`` — the assembly.

Output is written to the ``QUAST`` folder.


Annotation with Trinotate
==========================

With QC passed, annotate the assembly. Broadly there are two approaches:

#. **Ab initio gene finding** — predicting genes from sequence alone
   (AUGUSTUS, GeneMark, Glimmer); this is what BUSCO and QUAST used above.
#. **Homology-based searching** — what **Trinotate** uses, comparing against
   known databases (SwissProt/UniProt, PFAM, Gene Ontology, EggNOG).

Trinotate does **not** predict gene locations; it annotates the genes you give
it, in FASTA format. Here the assembly already is that gene set, so we supply
it directly. The annotation steps are in ``annotate.sh``.

.. note::

   As of March 2024, Trinotate is no longer under active development.

Load the Trinotate environment and make a safety copy of the assembly:

.. code-block:: bash

   module purge
   module load gencore
   module load Miniconda3/4.7.10
   source activate /scratch/gencore/conda3/envs/Trinotate/
   export PATH="/scratch/gencore/software/Trinotate-Trinotate-v4.0.2:$PATH"
   export TRINOTATE_DATA_DIR=/scratch/gencore/trinotate_addons/databases/4.0.2

   cp trinity_assembly.Trinity.fasta transcripts.fasta

Predict ORFs with TransDecoder
------------------------------

Identify the longest open reading frames in the transcripts:

.. code-block:: bash

   TransDecoder.LongOrfs -t transcripts.fasta

These ORFs are inferred from sequence alone. To add supporting evidence,
search the longest ORFs against SwissProt with BLASTP and against PFAM with
HMMER before the final prediction:

.. code-block:: bash

   blastp \
       -query transcripts.fasta.transdecoder_dir/longest_orfs.pep \
       -db /scratch/gencore/trinotate_addons/databases/4.0.2/uniprot_sprot.pep \
       -max_target_seqs 1 \
       -outfmt 6 \
       -evalue 1e-5 \
       -num_threads 28 > blastp.outfmt6

   hmmscan \
       --cpu 28 \
       --domtblout pfam.domtblout \
       /scratch/gencore/trinotate_addons/databases/4.0.2/Pfam-A.hmm \
       transcripts.fasta.transdecoder_dir/longest_orfs.pep

Re-run prediction, retaining ORFs with homology hits:

.. code-block:: bash

   TransDecoder.Predict \
       -t transcripts.fasta \
       --retain_pfam_hits pfam.domtblout \
       --retain_blastp_hits blastp.outfmt6

This produces ``transcripts.fasta.transdecoder.{bed,cds,gff3,pep}``.

Initialize and load the Trinotate database
-------------------------------------------

The data folder already includes a built ``myTrinotate.sqlite``. You only need
to create it once; the command below is shown for reference (it downloads and
builds PFAM, UniProt/SwissProt, EggNOG, etc.) and can be reused across
projects.

.. code-block:: bash

   Trinotate --create \
       --db myTrinotate.sqlite \
       --trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2 \
       --use_diamond

Import the transcripts and their predicted protein translations:

.. code-block:: bash

   Trinotate \
       --db myTrinotate.sqlite \
       --init \
       --gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
       --transcript_fasta transcripts.fasta \
       --transdecoder_pep transcripts.fasta.transdecoder.pep

Run Trinotate
-------------

.. code-block:: bash

   Trinotate --db myTrinotate.sqlite --CPU 28 \
       --transcript_fasta transcripts.fasta \
       --transdecoder_pep transcripts.fasta.transdecoder.pep \
       --trinotate_data_dir /scratch/gencore/trinotate_addons/databases/4.0.2/ \
       --run "swissprot_blastp swissprot_blastx pfam tmhmmv2 EggnogMapper" \
       --use_diamond

Flags:

- ``--CPU 28`` — use 28 CPUs.
- ``--transcript_fasta transcripts.fasta`` — the assembly as input.
- ``--transdecoder_pep ...`` — the TransDecoder ORFs.
- ``--trinotate_data_dir ...`` — the database directory used at init.
- ``--run "..."`` — the analysis modules to run.
- ``--use_diamond`` — use DIAMOND instead of NCBI BLAST, which is much faster.

Results load into the SQLite database automatically.

SignalP for Gram-positive bacteria
-----------------------------------

``SignalP6`` predicts signal peptides and cleavage sites. Trinotate runs it in
*eukaryotic* mode by default, which is wrong for our Gram-positive bacterium.
Run it independently in ``gram+`` mode, then load the results:

.. code-block:: bash

   /scratch/gencore/trinotate_addons/signalp-4.1/signalp \
       -f short \
       -t gram+ \
       -n signalp.out \
       transcripts.fasta.transdecoder.pep

   Trinotate --db myTrinotate.sqlite --LOAD_signalp signalp.out

Generate the annotation report
-------------------------------

Finally, export a tab-delimited report of all annotations.

.. warning::

   In the Trinotate version used here, report generation produces an **empty**
   report (this may be specific to our installation). The workaround is to
   generate the report with an earlier Trinotate version. Switch environments,
   then run the report command:

   .. code-block:: bash

      conda deactivate
      module purge
      module load gencore/1
      module load gencore_rnaseq/1.0
      module load gencore_annotation
      module load gencore_dev
      module load gencore_trinity/1.0

      Trinotate myTrinotate.sqlite report > trinotate_annotation_report.xls
      import_transcript_names.pl myTrinotate.sqlite trinotate_annotation_report.xls

   This is exactly why a ``conda``-managed stack is valuable: switching
   versions is a module swap, not a reinstall.

Open ``trinotate_annotation_report.xls`` in a spreadsheet to inspect the
annotations.


Expression quantification
=========================

The assembly is built, assessed, and annotated. The final stage quantifies
expression per sample and compares conditions — essentially a reference-based
RNA-Seq analysis using our own assembly as the reference. Trinity ships helper
scripts that integrate cleanly with the outputs so far, so we stay within
Trinity; ``quantify.sh`` holds these steps.

Trinity supports several differential-expression packages:

- edgeR — `<http://bioconductor.org/packages/release/bioc/html/edgeR.html>`_
- DESeq2 — `<http://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_
- limma/voom — `<http://bioconductor.org/packages/release/bioc/html/limma.html>`_
- ROTS — `<http://www.btk.fi/research/research-groups/elo/software/rots/>`_

.. warning::

   This dataset has only **two** replicates per condition, not three. With
   duplicates you can compute fold changes but cannot perform robust
   statistical testing of differential expression. Treat the DGE output here
   as illustrative.

Switch back to the Trinity stack (now also loading RSEM):

.. code-block:: bash

   module purge
   module load all
   module load gencore/2
   module load trinity/2.15.1
   module load salmon/0.14.2
   module load numpy/1.26.4
   module load bowtie2/2.4.5
   module load samtools/1.10
   module load hisat2/2.2.1
   module load rsem/1.3.3

Estimate abundance per sample
-----------------------------

A counts matrix records how much of each transcript each sample expresses.
There are alignment-based methods (RSEM, HISAT2, Bowtie2) and alignment-free,
k-mer methods (Kallisto, Salmon); a benchmark is here
(`<https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4869-5>`_).
We use **RSEM** with **Bowtie2** via Trinity's wrapper:

.. code-block:: bash

   align_and_estimate_abundance.pl \
       --transcripts trinity_assembly.Trinity.fasta \
       --seqType fq \
       --samples_file cond.txt \
       --est_method RSEM \
       --thread_count 28 \
       --gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
       --prep_reference \
       --aln_method bowtie2

``--prep_reference`` indexes the assembly with Bowtie2 first; skip it if the
assembly is already indexed. Output is per sample — one folder per replicate.
For example:

.. code-block:: text

   ls -1 fructose_rep1/
   bowtie2.bam
   bowtie2.bam.for_rsem.bam
   bowtie2.bam.ok
   RSEM.genes.results
   RSEM.isoforms.results
   RSEM.isoforms.results.ok
   RSEM.stat

The two files of interest are ``RSEM.genes.results`` and
``RSEM.isoforms.results`` — quantification at gene and transcript level.

Build expression matrices
--------------------------

Combine the per-replicate results into experiment-wide matrices with
``abundance_estimates_to_matrix.pl``:

.. code-block:: bash

   abundance_estimates_to_matrix.pl \
       --est_method RSEM \
       --gene_trans_map trinity_assembly.Trinity.fasta.gene_trans_map \
       --name_sample_by_basedir \
       --out_prefix RSEM.out \
       fructose_rep1/RSEM.isoforms.results \
       fructose_rep2/RSEM.isoforms.results \
       glucose_rep1/RSEM.isoforms.results \
       glucose_rep2/RSEM.isoforms.results \
       pyruvate_rep1/RSEM.isoforms.results \
       pyruvate_rep2/RSEM.isoforms.results

.. note::

   Supply the **isoform** results, not the gene results — this step computes
   both gene- and isoform-level matrices in one pass, and supplying the gene
   files instead raises an error.

This produces:

.. code-block:: text

   RSEM.out.gene.counts.matrix
   RSEM.out.gene.TMM.EXPR.matrix
   RSEM.out.gene.TPM.not_cross_norm
   RSEM.out.gene.TPM.not_cross_norm.runTMM.R
   RSEM.out.gene.TPM.not_cross_norm.TMM_info.txt
   RSEM.out.isoform.counts.matrix
   RSEM.out.isoform.TMM.EXPR.matrix
   RSEM.out.isoform.TPM.not_cross_norm
   RSEM.out.isoform.TPM.not_cross_norm.runTMM.R
   RSEM.out.isoform.TPM.not_cross_norm.TMM_info.txt

Differential expression with DESeq2
-----------------------------------

Run DGE with Trinity's ``run_DE_analysis.pl``, using DESeq2:

.. code-block:: bash

   module load bioconductor-deseq2/1.32.0

   run_DE_analysis.pl \
       --matrix RSEM.out.gene.counts.matrix \
       --method DESeq2 \
       --samples_file dge.txt \
       --output deseq2 \
       --contrasts comp.txt

Results are written to the ``deseq2`` folder.

.. note::

   As an alternative to the command line, you can use the NYUAD NASQAR2
   instance (`<http://nasqar2.abudhabi.nyu.edu/>`_) and its DESeq2 app
   (`<http://nasqar2.abudhabi.nyu.edu/deseq2shiny/>`_). For more downstream
   options within Trinity (e.g. extracting and clustering differentially
   expressed transcripts), see
   `<https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression>`_.


Visualization
=============

With the results generated, you can bring everything together visually:

- **TrinotateWeb** presents the assembly, annotations, and expression in one
  interface — `<https://github.com/Trinotate/Trinotate/wiki/TrinotateWeb>`_.
- **IGV** loads the sorted BAM and assembly to inspect read support for
  individual transcripts.

This concludes the workshop. For questions or comments, raise an issue on the
workshop repository.
