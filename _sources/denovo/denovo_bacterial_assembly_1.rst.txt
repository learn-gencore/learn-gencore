.. _denovo-bacterial-assembly:

=====================================================
De novo Bacterial Genome Assembly
=====================================================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

In this tutorial we will go from raw paired-end Illumina short reads to an
annotated, quality-assessed bacterial genome assembly. We will cover:

#. Read quality assessment with **FastQC**
#. Quality and adapter trimming with **fastp**
#. Re-assessment of trimmed reads with **FastQC**
#. *De novo* assembly with **SPAdes**
#. Assembly quality assessment with **QUAST** and **BUSCO**
#. Genome annotation with **Prokka**
#. Extra value-added steps: species/contamination check (Kraken2),
   MLST typing, coverage/depth assessment, and a MultiQC summary report

This tutorial follows the same structure and spirit as the NYU GenCore
`de novo transcriptome assembly tutorial <https://learn.gencore.bio.nyu.edu/denovo/transcriptome.html>`_,
adapted here for bacterial genome assembly.

.. note::

   Commands below are shown both for users of our HPC (via ``module load``)
   and for users running locally (via ``conda``/``mamba``, ``Docker``, or
   pip/standalone binaries). Pick whichever column applies to you.

Prerequisites
=============

* Basic familiarity with the Linux command line
* An HPC account with Slurm/PBS access, **or** a local machine/laptop with
  conda/mamba and at least 8 GB RAM and ~10 GB free disk space
* The example dataset (see :ref:`data-download` below)

.. _data-download:

Example Dataset
================

For this tutorial we use a real, publicly available paired-end Illumina
MiSeq bacterial whole-genome sequencing dataset: a methicillin-resistant
*Staphylococcus aureus* (MRSA) isolate from Hikichi *et al.* 2019,
*"Complete Genome Sequences of Eight Methicillin-Resistant Staphylococcus
aureus Strains Isolated from Patients in Japan"* (`Microbiology Resource
Announcements <https://doi.org/10.1128/mra.01212-19>`_).

This dataset is small (genome ~2.9 Mb), well characterized, has a closed
reference genome available for benchmarking, and is already used in a
peer-reviewed `Galaxy Training Network tutorial
<https://training.galaxyproject.org/training-material/topics/assembly/tutorials/mrsa-illumina/tutorial.html>`_,
so results in this tutorial can be cross-checked against a known-good
analysis.

+--------------------+---------------------------------------------+
| Run accession      | DRR187559 (DDBJ/SRA/ENA, BioProject DRP004939)|
+--------------------+---------------------------------------------+
| Organism           | *Staphylococcus aureus* (MRSA)               |
+--------------------+---------------------------------------------+
| Platform           | Illumina MiSeq, 2 x 300 bp paired-end        |
+--------------------+---------------------------------------------+
| Approx. genome size| 2,914,567 bp                                 |
+--------------------+---------------------------------------------+
| Reads              | ``DRR187559_1.fastq.gz``, ``DRR187559_2.fastq.gz`` |
+--------------------+---------------------------------------------+

.. admonition:: Shared drive download
   :class: tip

   For convenience, these reads are mirrored on our group's shared
   Google Drive folder so you don't need to query SRA/ENA yourself:

   * Folder: https://drive.google.com/drive/folders/1cSxLTe7R_YDWq2DCJnuiRFTi-jAHtzuJ?usp=sharing
   * Contains: ``DRR187559_1.fastq.gz`` (forward/R1) and
     ``DRR187559_2.fastq.gz`` (reverse/R2)

   Download the two files from the folder above (click each file ->
   Download), or, if you'd rather use the command line, install
   `gdown <https://github.com/wkentaro/gdown>`_ and pull the whole folder:

   .. code-block:: bash

      mkdir -p ~/denovo_tutorial/raw_data && cd ~/denovo_tutorial/raw_data

      pip install --user gdown
      gdown --folder "https://drive.google.com/drive/folders/1cSxLTe7R_YDWq2DCJnuiRFTi-jAHtzuJ" -O .

   This should leave you with ``DRR187559_1.fastq.gz`` and
   ``DRR187559_2.fastq.gz`` in ``~/denovo_tutorial/raw_data``.

   .. note::

      ``gdown`` may not be available/loadable as a module on every HPC.
      If it isn't, just download the two files manually through the Drive
      web UI on your laptop and ``scp``/``rsync`` them to the cluster, e.g.:

      .. code-block:: bash

         rsync -avP DRR187559_1.fastq.gz DRR187559_2.fastq.gz \
             your_username@hpc.example.edu:~/denovo_tutorial/raw_data/

Downloading the raw data yourself (alternative)
------------------------------------------------

If the shared-drive links above are unavailable, the same reads can be
obtained directly from public archives.

**Option A — SRA Toolkit (works for any HPC/local install):**

.. code-block:: bash

   # HPC
   module purge && module load all gencore/3 sratoolkit
   prefetch DRR187559
   fasterq-dump --split-files DRR187559
   gzip DRR187559_1.fastq DRR187559_2.fastq

   # Local (conda)
   conda create -n sra -c bioconda -c conda-forge sra-tools -y
   conda activate sra
   prefetch DRR187559
   fasterq-dump --split-files DRR187559
   gzip DRR187559_*.fastq

**Option B — Direct from ENA (no tool install needed):**

Browse to https://www.ebi.ac.uk/ena/browser/view/DRR187559, open the
"Download" tab, and copy the FASTQ FTP links, e.g.:

.. code-block:: bash

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR187/059/DRR187559/DRR187559_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR187/059/DRR187559/DRR187559_2.fastq.gz

.. note::

   Double check the exact path on the ENA browser page before using the
   ``wget`` command above — ENA's FTP directory numbering depends on the
   exact accession digits and can differ from the example shown here.

Reference genome (for QUAST/annotation comparison, optional)
--------------------------------------------------------------

Any closed, annotated *S. aureus* RefSeq assembly works well as a
reference for QUAST comparison purposes (it does not need to be the exact
same strain). For example, *S. aureus* NCTC 8325
(``GCF_000013425.1``) can be fetched via NCBI datasets:

.. code-block:: bash

   datasets download genome accession GCF_000013425.1 --include genome,gff3
   unzip ncbi_dataset.zip

Working with your own data
----------------------------

Everything in this tutorial works directly on your own paired-end
``*_R1.fastq.gz`` / ``*_R2.fastq.gz`` (or ``*_1.fastq.gz`` / ``*_2.fastq.gz``)
files — simply substitute your filenames wherever ``DRR187559_1/2`` appears.

Step 1: Quality Assessment with FastQC
=========================================

Before any processing, we inspect the raw reads to check base quality,
adapter content, GC content, and duplication levels.

.. code-block:: bash

   mkdir -p ~/denovo_tutorial/qc/fastqc_raw
   cd ~/denovo_tutorial

**On the HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 fastqc
   fastqc -o qc/fastqc_raw -t 4 raw_data/DRR187559_1.fastq.gz raw_data/DRR187559_2.fastq.gz

**Locally (conda):**

.. code-block:: bash

   conda create -n qc -c bioconda -c conda-forge fastqc -y
   conda activate qc
   fastqc -o qc/fastqc_raw -t 4 raw_data/DRR187559_1.fastq.gz raw_data/DRR187559_2.fastq.gz

**Locally (Docker, no install needed):**

.. code-block:: bash

   docker run --rm -v $PWD:/data -w /data biocontainers/fastqc:v0.11.9_cv8 \
       fastqc -o qc/fastqc_raw -t 4 raw_data/DRR187559_1.fastq.gz raw_data/DRR187559_2.fastq.gz

Open ``qc/fastqc_raw/DRR187559_1_fastqc.html`` in a browser and review:

* **Per base sequence quality** — should be mostly green (Phred > 28);
  Illumina runs commonly show some quality drop-off toward the 3' end.
* **Adapter content** — flags if adapters need trimming.
* **Sequence duplication levels** and **overrepresented sequences** —
  high values can indicate PCR duplication or contamination.

Step 2: Quality and Adapter Trimming with fastp
==================================================

We use ``fastp`` to remove adapters, trim low-quality bases from read
ends, and discard reads that become too short to be useful for assembly.

.. code-block:: bash

   mkdir -p ~/denovo_tutorial/trimmed
   cd ~/denovo_tutorial

**On the HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 fastp
   fastp \
       -i raw_data/DRR187559_1.fastq.gz \
       -I raw_data/DRR187559_2.fastq.gz \
       -o trimmed/DRR187559_1.trimmed.fastq.gz \
       -O trimmed/DRR187559_2.trimmed.fastq.gz \
       --detect_adapter_for_pe \
       --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
       --length_required 30 \
       --thread 4 \
       --json trimmed/DRR187559.fastp.json \
       --html trimmed/DRR187559.fastp.html

**Locally (conda):**

.. code-block:: bash

   conda create -n trim -c bioconda -c conda-forge fastp -y
   conda activate trim
   fastp \
       -i raw_data/DRR187559_1.fastq.gz -I raw_data/DRR187559_2.fastq.gz \
       -o trimmed/DRR187559_1.trimmed.fastq.gz -O trimmed/DRR187559_2.trimmed.fastq.gz \
       --detect_adapter_for_pe \
       --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
       --length_required 30 --thread 4 \
       --json trimmed/DRR187559.fastp.json --html trimmed/DRR187559.fastp.html

Parameter notes:

* ``--detect_adapter_for_pe`` — auto-detects adapters by overlap analysis
  of read pairs (no adapter FASTA needed).
* ``--cut_front``/``--cut_tail`` with a sliding window — trims low-quality
  ends (mirrors the approach commonly used with Trimmomatic).
* ``--length_required 30`` — discards reads shorter than 30 bp post-trimming,
  since very short reads complicate assembly graphs.

Step 3: Re-Assessment of Trimmed Reads
=========================================

Re-run FastQC on the trimmed reads to confirm improvement, and inspect the
fastp HTML report directly (it already contains before/after comparison
panels).

.. code-block:: bash

   mkdir -p qc/fastqc_trimmed

**HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 fastqc
   fastqc -o qc/fastqc_trimmed -t 4 \
       trimmed/DRR187559_1.trimmed.fastq.gz trimmed/DRR187559_2.trimmed.fastq.gz

**Locally:**

.. code-block:: bash

   conda activate qc
   fastqc -o qc/fastqc_trimmed -t 4 \
       trimmed/DRR187559_1.trimmed.fastq.gz trimmed/DRR187559_2.trimmed.fastq.gz

Compare ``qc/fastqc_raw`` vs ``qc/fastqc_trimmed`` (or simply open
``trimmed/DRR187559.fastp.html``) and confirm:

* Per-base quality scores improved or stayed high across the full read length
* Adapter content dropped to ~0%
* Read length distribution shifted down slightly (expected, due to trimming)
* GC content is unchanged (trimming should not bias GC%)

Step 4: De Novo Assembly with SPAdes
========================================

With clean reads in hand, we assemble the genome using **SPAdes**, a
widely used assembler for bacterial/small genomes.

.. code-block:: bash

   mkdir -p ~/denovo_tutorial/assembly

**On the HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 spades
   spades.py \
       -1 trimmed/DRR187559_1.trimmed.fastq.gz \
       -2 trimmed/DRR187559_2.trimmed.fastq.gz \
       --isolate \
       -o assembly/DRR187559_spades \
       -t 8 -m 32

**Locally (conda):**

.. code-block:: bash

   conda create -n spades -c bioconda -c conda-forge spades -y
   conda activate spades
   spades.py \
       -1 trimmed/DRR187559_1.trimmed.fastq.gz \
       -2 trimmed/DRR187559_2.trimmed.fastq.gz \
       --isolate \
       -o assembly/DRR187559_spades \
       -t 8 -m 16

Notes:

* ``--isolate`` mode is recommended by the SPAdes authors specifically
  for high-coverage isolate bacterial genomes (it disables some
  metagenome/low-coverage-oriented heuristics and tends to produce a
  cleaner assembly than the default mode).
* Adjust ``-t`` (threads) and ``-m`` (memory in GB) to match your
  allocation/machine.
* The key output is ``assembly/DRR187559_spades/contigs.fasta``
  (and ``scaffolds.fasta``).

On an HPC with Slurm, submit as a batch job, e.g.:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=spades_assembly
   #SBATCH --cpus-per-task=8
   #SBATCH --mem=32G
   #SBATCH --time=04:00:00
   #SBATCH --output=spades_%j.log

   module purge && module load all gencore/3 spades
   spades.py -1 trimmed/DRR187559_1.trimmed.fastq.gz \
             -2 trimmed/DRR187559_2.trimmed.fastq.gz \
             --isolate -o assembly/DRR187559_spades -t 8 -m 32

Step 5: Assembly Quality Assessment
=======================================

QUAST
-----

**QUAST** reports assembly statistics (N50, number of contigs, total
length, GC content, misassemblies if a reference is given, etc.).

**HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 quast
   quast.py assembly/DRR187559_spades/contigs.fasta \
       -o qc/quast_DRR187559 \
       -t 4 \
       -r reference/GCF_000013425.1_genomic.fna \
       -g reference/GCF_000013425.1_genomic.gff

**Locally:**

.. code-block:: bash

   conda create -n quast -c bioconda -c conda-forge quast -y
   conda activate quast
   quast.py assembly/DRR187559_spades/contigs.fasta -o qc/quast_DRR187559 -t 4

.. note::

   The ``-r``/``-g`` reference arguments are optional but recommended when
   a reasonably close reference genome is available — they enable
   misassembly detection and genome-fraction metrics. Omit them to run
   QUAST reference-free.

Key metrics to check in ``qc/quast_DRR187559/report.html``:

* **Number of contigs** — fewer, larger contigs generally indicate a
  better assembly
* **N50** — the contig length at which 50% of the assembly is contained
  in contigs of that size or larger; higher is better
* **Total length** — should be close to the expected genome size
  (~2.9 Mb for this MRSA isolate)
* **GC content** — should match the expected species GC% (~32-33% for
  *S. aureus*)
* **L50, largest contig, # misassemblies** (if reference supplied)

BUSCO
-----

**BUSCO** assesses assembly completeness by searching for a curated set
of near-universal single-copy orthologs expected in your lineage
(e.g. ``bacteria_odb10`` or a more specific lineage such as
``bacillales_odb10``).

**HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 busco
   busco \
       -i assembly/DRR187559_spades/contigs.fasta \
       -o busco_DRR187559 \
       --out_path qc \
       -l bacteria_odb10 \
       -m genome \
       -c 4

**Locally (conda):**

.. code-block:: bash

   conda create -n busco -c bioconda -c conda-forge busco -y
   conda activate busco
   busco -i assembly/DRR187559_spades/contigs.fasta \
       -o busco_DRR187559 --out_path qc -l bacteria_odb10 -m genome -c 4

Check ``qc/busco_DRR187559/short_summary*.txt`` for the proportion of
**Complete (Single + Duplicated)**, **Fragmented**, and **Missing** BUSCOs.
A good bacterial assembly typically shows >95% complete BUSCOs.

Step 6: Genome Annotation
=============================

We annotate the assembled contigs with **Prokka**, a fast, widely used
bacterial annotation pipeline (predicts CDS, rRNA, tRNA, and more, and
outputs GFF3/GenBank/FASTA protein files).

**HPC:**

.. code-block:: bash

   module purge && module load all gencore/3 prokka
   prokka \
       --outdir annotation/DRR187559_prokka \
       --prefix DRR187559 \
       --genus Staphylococcus --species aureus \
       --cpus 8 \
       assembly/DRR187559_spades/contigs.fasta

**Locally (conda):**

.. code-block:: bash

   conda create -n prokka -c bioconda -c conda-forge prokka -y
   conda activate prokka
   prokka --outdir annotation/DRR187559_prokka --prefix DRR187559 \
       --genus Staphylococcus --species aureus --cpus 8 \
       assembly/DRR187559_spades/contigs.fasta

Key outputs in ``annotation/DRR187559_prokka/``:

* ``DRR187559.gff`` — annotation in GFF3 format (genes, CDS, RNA features)
* ``DRR187559.faa`` — predicted protein sequences
* ``DRR187559.ffn`` — nucleotide sequences of predicted genes
* ``DRR187559.txt`` — summary statistics (number of CDS, rRNAs, tRNAs, etc.)

Step 7: Additional Value-Added Steps
========================================

The following steps are optional but commonly add real value to a
bacterial assembly workflow:

7a. Species/Contamination Check with Kraken2
------------------------------------------------

Confirms the reads/assembly actually belong to the expected species and
flags contamination (common in bacterial WGS from clinical/environmental
samples).

.. code-block:: bash

   module purge && module load all gencore/3 kraken2   # HPC
   # or: conda create -n kraken2 -c bioconda -c conda-forge kraken2 -y

   kraken2 --db /path/to/kraken2_standard_db \
       --paired trimmed/DRR187559_1.trimmed.fastq.gz trimmed/DRR187559_2.trimmed.fastq.gz \
       --report qc/kraken2_DRR187559.report \
       --output qc/kraken2_DRR187559.out \
       --threads 4

7b. MLST Typing
------------------

Multi-locus sequence typing assigns a sequence type (ST), useful for
strain tracking and outbreak investigations.

.. code-block:: bash

   module purge && module load all gencore/3 mlst   # HPC
   # or: conda create -n mlst -c bioconda -c conda-forge mlst -y

   mlst assembly/DRR187559_spades/contigs.fasta > qc/DRR187559.mlst.tsv

7c. Coverage / Depth Assessment
-----------------------------------

Map trimmed reads back to the assembly to confirm even coverage across
contigs (low/uneven coverage in regions can flag misassemblies).

.. code-block:: bash

   module purge && module load all gencore/3 bwa samtools   # HPC

   bwa index assembly/DRR187559_spades/contigs.fasta
   bwa mem -t 8 assembly/DRR187559_spades/contigs.fasta \
       trimmed/DRR187559_1.trimmed.fastq.gz trimmed/DRR187559_2.trimmed.fastq.gz \
       | samtools sort -@ 4 -o qc/DRR187559.sorted.bam
   samtools index qc/DRR187559.sorted.bam
   samtools depth -a qc/DRR187559.sorted.bam > qc/DRR187559.depth.txt

7d. Antimicrobial Resistance (AMR) Gene Screening
------------------------------------------------------

Since this is an MRSA isolate, screening for resistance genes adds
biological context (this also nicely follows on from the GTN tutorial
referenced above).

.. code-block:: bash

   module purge && module load all gencore/3 abricate   # HPC
   # or: conda create -n abricate -c bioconda -c conda-forge abricate -y

   abricate --db resfinder assembly/DRR187559_spades/contigs.fasta > qc/DRR187559.amr.tsv

7e. Consolidated MultiQC Report
-----------------------------------

Aggregate FastQC, fastp, QUAST, and BUSCO outputs into a single summary
HTML report.

.. code-block:: bash

   module purge && module load all gencore/3 multiqc   # HPC
   # or: conda create -n multiqc -c bioconda -c conda-forge multiqc -y

   multiqc qc/ trimmed/ -o qc/multiqc_report

Summary
=========

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Step
     - Tool
     - Purpose
   * - 1. Raw QC
     - FastQC
     - Assess raw read quality, adapter content, duplication
   * - 2. Trimming
     - fastp
     - Remove adapters, trim low-quality bases, filter short reads
   * - 3. Trimmed QC
     - FastQC / fastp report
     - Confirm trimming improved read quality
   * - 4. Assembly
     - SPAdes
     - Reconstruct the genome from short reads
   * - 5. Assembly QC
     - QUAST, BUSCO
     - Assess contiguity (N50, # contigs) and completeness (% BUSCOs)
   * - 6. Annotation
     - Prokka
     - Predict genes, CDS, rRNA, tRNA
   * - 7. Extras
     - Kraken2, MLST, BWA/samtools, ABRicate, MultiQC
     - Contamination check, strain typing, coverage, AMR genes, summary report

Further Reading
==================

* SPAdes manual: https://github.com/ablab/spades
* QUAST manual: http://quast.sourceforge.net/
* BUSCO user guide: https://busco.ezlab.org/
* Prokka: https://github.com/tseemann/prokka
* fastp: https://github.com/OpenGene/fastp
* Galaxy Training — MRSA assembly tutorial (data source for this tutorial):
  https://training.galaxyproject.org/training-material/topics/assembly/tutorials/mrsa-illumina/tutorial.html
