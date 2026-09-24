FastQ format
============

For more information on the FastQ format, see the `MAQ FastQ
documentation <http://maq.sourceforge.net/fastq.shtml>`_.

FastQ is the standard format for storing the raw reads that come off a
sequencer. It extends the FastA format by pairing each base with a
quality score that describes how confident the sequencer was in that
basecall.

Format
------

Each read in a FastQ file is described by four lines:

1. A line beginning with ``@`` that holds the sequence identifier and an
   optional description.
2. The raw sequence (the bases).
3. A line beginning with ``+``, optionally followed by the same
   identifier again.
4. The quality scores, encoded as ASCII characters, one per base in
   line 2.

Example
-------

A single FastQ read looks like this:

.. code-block:: text

   @HWI-ST911:111:C0N4WACXX:5:1101:2249:2216 1:Y:18:TTAGGC
   TTAGGCAGGACAGCTCAGGGCATGAAGTTGTTAATTCAGGACAGGGCATGT
   +
   CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHIJJJ

The identifier line packs a lot of information. Breaking the first
example apart:

* ``HWI-ST911`` — the unique instrument name
* ``111`` — the run id
* ``C0N4WACXX`` — the flowcell id
* ``5`` — the flowcell lane
* ``1101`` — the tile number within the lane
* ``2249`` — the x-coordinate of the cluster within the tile
* ``2216`` — the y-coordinate of the cluster within the tile
* ``1`` — the member of a pair (1 or 2; paired-end reads only)
* ``Y`` — whether the read failed the filter (Y = filtered out, N =
  passed)
* ``18`` — control bits (0 when none are on)
* ``TTAGGC`` — the index (barcode) sequence

What software use FastQ?
------------------------

FastQ is the entry point for almost every downstream tool. A few
examples:

* Aligners such as `Bowtie 2
  <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ and `TopHat2
  <https://ccb.jhu.edu/software/tophat/index.shtml>`_.
* Assemblers such as `Velvet
  <https://github.com/dzerbino/velvet>`_ and `SPAdes
  <https://github.com/ablab/spades>`_.
* Quality control and trimming tools such as `Trimmomatic
  <http://www.usadellab.org/cms/?page=trimmomatic>`_ and `FastQC
  <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_.

How are these files generated?
------------------------------

FastQ files are produced by the sequencer's basecalling software, which
converts the raw signal measured for each cluster into a sequence of
bases and their associated quality scores.
