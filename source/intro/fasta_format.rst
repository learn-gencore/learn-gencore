FastA format
============

For more information on the FastA format, see the `NCBI BLAST
documentation
<https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp>`_.

FastA is one of the simplest and most widely used file formats in
bioinformatics. It is a text-based format for representing nucleotide or
amino acid sequences, in which each residue is represented by a single
letter.

Format
------

A FastA file is made up of one or more records. Each record has two
parts:

* A single **header line** that begins with a greater-than symbol
  (``>``), followed by an identifier and an optional description.
* One or more **sequence lines** containing the nucleotide or amino acid
  sequence.

Example
-------

Below is the beginning of a FastA record for chromosome 1:

.. code-block:: text

   >Chr1
   CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAAT
   CTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTT
   CTCTGGTTGAAAATCATTGTGTATATAATGATAATTTTATCGTTTTTATGTAATTGCTTA
   TTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGT

The header line names the sequence (here, ``Chr1``), and every line
after it, up to the next ``>``, is the sequence itself.

Software that use FastA format
------------------------------

Many bioinformatics tools require input in FastA format. Aligners, for
example, need the reference genome supplied as a FastA file, and
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_ searches a query
sequence (in FastA format) against a database of known sequences.

How are these files generated?
------------------------------

FastA files usually come from public sequence databases (for example a
reference genome or a set of gene models downloaded from Ensembl or
NCBI), or are produced as the output of an assembly program.
