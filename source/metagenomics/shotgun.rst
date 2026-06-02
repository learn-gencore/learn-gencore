Series of Metagenomics Workshop 
================================

Shotgun
----------------------
 
During this course, participants will learn about the concepts and steps involved in analyzing
the metagenomic data obtained from short read high throughput sequencing.
 
Participants are expected to have some basic knowledge of metagenome/microbiome as well as
some biological knowledge on the subject matter.
 
Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance
Computing environment, it is possible to run this workshop on any other system provided that the
necessary software is installed, and the data is uploaded to that environment. It will also be
necessary to change the input and output directories/files to accommodate such an environment.
 

 
Connecting to the HPC
----------------------
 
Using a Mac Machine
~~~~~~~~~~~~~~~~~~~
 
1. Open the **Terminal** app and type:
 
   .. code-block:: bash
 
      ssh NetID@jubail.abudhabi.nyu.edu
 
2. Enter your NYU password when prompted and hit enter.
3. Once logged in, navigate to your personal **SCRATCH** directory or any other subdirectory
   where you want to run this workshop.
 
Using a Windows Machine
~~~~~~~~~~~~~~~~~~~~~~~
 
1. Open the **PuTTY** app, and fill out the fields as follows:
 
   - **Host name:** ``jubail.abudhabi.nyu.edu``
   - **Port:** ``22``
   - Click **Open**.
 
2. Enter your NetID and password when prompted.
3. Once logged in, navigate to your personal **SCRATCH** directory or any other subdirectory
   where you want to run this workshop.
 
----
 
MODULE 1
--------
 
The Dataset
~~~~~~~~~~~
 
The shotgun metagenomic sequencing data used for this hands-on training was obtained from a
recently published study in *Microbiome* by Jeilu et al. 2025
(`https://doi.org/10.1186/s40168-025-02276-7 <https://doi.org/10.1186/s40168-025-02276-7>`_).
The authors deposited the data in the NCBI under the BioProject accession number **PRJNA1228129**.
 
Briefly, the study was undertaken for metagenomic profiling of airborne microbial communities
from aircraft filters and face masks. This study compared the following four groups:
 
- **Aircraft filter** — aircraft cabin filter removed from a commercial airplane after 8039 h of usage
- **Travel** — face masks worn by air travelers
- **Hospital** — face masks worn by healthcare professionals in a hospital
- **Control** — unworn
 
For the sake of time and resources, we will be using the following two groups in this workshop
to analyze microbial profiling:
 
- Aircraft filter
- Control
 
We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop; however,
you can run all of the analysis on any stand-alone machine (server, personal laptop/desktop, etc.)
provided that the necessary software packages are pre-installed.
 
Setting Up the Working Directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
.. code-block:: bash
 
   mkdir -p /scratch/$USER/metagenomics/
   mkdir -p /scratch/$USER/metagenomics/raw_fastqs/
   cd /scratch/$USER/metagenomics/
   cp -r /scratch/gencore/metagenomics/* .
 
Quality Control and Trimming (QC/QT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
We will check the quality of reads, then perform quality trimming to remove low-quality reads,
adapter sequences, and host (human) sequences. This step produces clean, high-quality reads free
from host contamination for downstream analysis.
 
We will use **KneadData**, part of the Biobakery suite
(`https://github.com/biobakery/kneaddata <https://github.com/biobakery/kneaddata>`_).
 
Load the KneadData module:
 
.. code-block:: bash
 
   module purge
   module load all
   module load gencore/3
   module load kneaddata/0.12.4
 
Run ``kneaddata --help`` to explore available options, then run the tool:
 
.. code-block:: bash
 
   kneaddata \
     -i1 raw_fastqs/AF1_1.fastq \
     -i2 raw_fastqs/AF1_2.fastq \
     --reference-db /scratch/Reference_Genomes/Public/Metagenomic/kneaddata_hg39/ \
     --reference-db /scratch/Reference_Genomes/Public/Metagenomic/kneaddata_phix/ \
     --trimmomatic-options="SLIDINGWINDOW:4:20" --trimmomatic-options="MINLEN:50" \
     --run-trim-repetitive \
     --output analysis/kneaddata/ \
     --output-prefix AF1 \
     --remove-intermediate-output \
     --fastqc fastqc \
     --threads 14
 
Run FastQC on the cleaned output:
 
.. code-block:: bash
 
   fastqc analysis/kneaddata/AF1_paired_1.fastq -t 14 -o analysis/kneaddata/fastqc
   fastqc analysis/kneaddata/AF1_paired_2.fastq -t 14 -o analysis/kneaddata/fastqc
 
.. note::
 
   The results of this run have been shared with you. Compare the results you obtain by running
   QC/QT with the provided ready-made results. Feel free to ask questions or discuss if you
   obtain very different results between the two sets.
 
----
 
MODULE 2
--------
 
Taxonomic Classification with Kraken2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
Load the required modules and run Kraken2:
 
.. code-block:: bash
 
   module purge && \
   module load all && \
   module load gencore/3 && \
   module load kraken2/2.17.1
 
   k2 classify \
     --db /scratch/Reference_Genomes/Public/Metagenomic/kraken2/bacteria \
     --threads 14 \
     --confidence 0.1 \
     --output analysis/kraken2/AF1.kraken \
     --unclassified-out analysis/kraken2/AF1_unbact#.fastq \
     --report analysis/kraken2/AF1_profile.txt \
     --report-zero-counts \
     --paired analysis/kneaddata/AF1_paired_1.fastq analysis/kneaddata/AF1_paired_2.fastq
 
Abundance Estimation with Bracken
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
.. code-block:: bash
 
   module purge && \
   module load all gencore/3 && \
   module load bracken/3.1
 
   bracken -d /scratch/Reference_Genomes/Public/Metagenomic/kraken2/bacteria \
     -i analysis/kraken2/AF1_profile.txt \
     -o analysis/bracken/AF1_S.txt \
     -r 100 -l S -t 10 > analysis/bracken/AF1_log.txt
 
Diversity Analysis
------------------
 
Load modules:
 
.. code-block:: bash
 
   module load all gencore/3
   module load R/4.3.1
   module load biopython/1.85
 
Alpha Diversity
~~~~~~~~~~~~~~~~
 
.. code-block:: bash
 
   python \
     /scratch/Reference_Genomes/Public/Metagenomic/KrakenTools/DiversityTools/alpha_diversity.py \
     -f analysis/bracken/AF1_S.txt \
     -a Sh \
     > analysis/alpha_diversity/AF1_Sh.txt
 
Beta Diversity
~~~~~~~~~~~~~~~
 
.. code-block:: bash
 
   python \
     /scratch/Reference_Genomes/Public/Metagenomic/KrakenTools/DiversityTools/beta_diversity.py \
     -i analysis/kraken2/*profile.txt --type kreport2 -l S \
     > analysis/beta_diversity/merged_beta_div.txt
 
Diversity Plots
~~~~~~~~~~~~~~~~
 
Set up the plots directory and merge diversity outputs:
 
.. code-block:: bash
 
   mkdir -p diversity_plots/
   cd diversity_plots/
   cp analysis/kraken2/*_profile.txt .
   cp analysis/bracken/*_Sh.txt .
 
   python /scratch/Reference_Genomes/Public/Metagenomic/KrakenTools/DiversityTools/beta_diversity.py \
     -i *_profile.txt --type kreport2 -l S > merged_beta_div.txt
 
   (echo -e "Sample\tDiversity"; \
     awk 'FNR==1{f=FILENAME; sub(/^.*\//,"",f); sub("_Sh.txt","",f)} \
          /Shannon/{print f"\t"$NF}' *_Sh.txt) > merged_alpha_div.txt
 
   rm *_profile.txt
   rm *_Sh.txt
 
Add group annotations to the alpha diversity table:
 
.. code-block:: bash
 
   Rscript -e '
     df <- read.table("merged_alpha_div.txt", header=TRUE, stringsAsFactors=FALSE)
     df$Group <- ifelse(grepl("^AF", df$Sample), "Filter",
                   ifelse(grepl("^C", df$Sample), "Control", "Other"))
     write.table(df, "merged_alpha_div_with_group.txt", sep="\t", quote=FALSE, row.names=FALSE)
   '
 
Generate a Shannon diversity boxplot with Wilcoxon test:
 
.. code-block:: bash
 
   Rscript -e '
     df <- read.table("merged_alpha_div_with_group.txt", header=TRUE, stringsAsFactors=FALSE)
     ymax <- max(df$Diversity) * 1.15
     png("diversity_boxplot_stats.png", width=800, height=600, res=120)
     boxplot(Diversity ~ Group, data=df, col=c("steelblue","orange"),
             ylab="Shannon Diversity", main="Shannon Diversity by Group",
             ylim=c(min(df$Diversity), ymax))
     stripchart(Diversity ~ Group, data=df, vertical=TRUE, method="jitter",
                pch=16, col="black", add=TRUE)
     w <- wilcox.test(Diversity ~ Group, data=df)
     p_txt <- paste0("Wilcoxon p = ", signif(w$p.value, 3))
     text(x=1.5, y=ymax, labels=p_txt)
     dev.off()
   '
 
Clean the beta diversity matrix:
 
.. code-block:: bash
 
   Rscript -e '
     args <- commandArgs(trailingOnly=TRUE)
     infile <- args[1]; outfile <- args[2]
     lines <- readLines(infile)
     samp_lines <- grep("^#", lines, value=TRUE)
     samps <- sub("^#\\d+\\s+","", samp_lines)
     samps <- sub("_profile.txt.*","", samps)
     start <- grep("^x\\s", lines)
     mat <- read.table(text=lines[start:length(lines)], header=TRUE,
                       stringsAsFactors=FALSE, check.names=FALSE)
     mat <- mat[, -1]
     mat[mat=="x.xxx"] <- NA
     mat <- as.matrix(sapply(mat, as.numeric))
     rownames(mat) <- samps; colnames(mat) <- samps
     mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
     diag(mat) <- 0
     write.table(mat, file=outfile, sep="\t", quote=FALSE, col.names=NA)
   ' merged_beta_div.txt cleaned_beta_matrix.txt
 
Generate a PCoA plot of beta diversity:
 
.. code-block:: bash
 
   Rscript -e '
     library(vegan)
     mat <- as.matrix(read.table("cleaned_beta_matrix.txt", header=TRUE,
                                  check.names=FALSE, row.names=1))
     d <- as.dist(mat)
     p <- cmdscale(d, k=2, eig=TRUE)
     scores <- as.data.frame(p$points)
     scores$Sample <- rownames(scores)
     scores$Group <- ifelse(grepl("^AF", scores$Sample), "AirFilter",
                      ifelse(grepl("^C", scores$Sample), "Control", "Other"))
     ad <- adonis2(d ~ Group, data=scores)
     pval <- ad$"Pr(>F)"[1]
     var_expl <- round(p$eig / sum(p$eig) * 100, 1)
     png("beta_pcoa.png", width=800, height=600, res=120)
     par(mar=c(5,5,4,2))
     cols <- c(AirFilter="steelblue", Control="orange", Other="gray")
     plot(scores$V1, scores$V2,
          xlab=paste0("PCoA1 (", var_expl[1], "%)"),
          ylab=paste0("PCoA2 (", var_expl[2], "%)"),
          pch=16, col=cols[scores$Group],
          main="PCoA of beta diversity")
     legend("topright", legend=unique(scores$Group),
            col=cols[unique(scores$Group)], pch=16, bty="n")
     text(x=min(scores$V1), y=max(scores$V2),
          labels=paste0("PERMANOVA p = ", signif(pval, 3)), adj=c(0,1))
     dev.off()
   '