BED format
==========

For more information on the BED format, see the `UCSC Genome Browser FAQ
<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_.

BED (Browser Extensible Data) is a simple, flexible, tab-delimited format
for describing genomic regions. Each line defines a feature by its
location on a chromosome, and optionally carries extra information such
as a name, a score, and display instructions.

Required fields
---------------

Every BED line must have at least the first three fields:

* **chrom** — the name of the chromosome or scaffold (for example
  ``chr1``).
* **chromStart** — the start position of the feature, counting from 0.
* **chromEnd** — the end position of the feature. The feature spans the
  bases from ``chromStart`` up to but not including ``chromEnd``.

A minimal BED file looks like this:

.. code-block:: text

   chr7	127471196	127472363
   chr7	127472363	127473530
   chr7	127473530	127474697

Optional fields
---------------

Nine further fields may follow the first three, in this fixed order:

* **name** — a label for the feature.
* **score** — a score between 0 and 1000, used for shading in a browser.
* **strand** — ``+`` or ``-``.
* **thickStart** — the position at which the feature starts being drawn
  thickly.
* **thickEnd** — the position at which the feature stops being drawn
  thickly.
* **itemRgb** — an ``R,G,B`` colour value for the feature.
* **blockCount** — the number of blocks (exons) in the feature.
* **blockSizes** — a comma-separated list of block sizes.
* **blockStarts** — a comma-separated list of block start positions,
  relative to ``chromStart``.

A BED file using several of the optional fields:

.. code-block:: text

   chr7	127471196	127472363	Pos1	0	+	127471196	127472363	255,0,0
   chr7	127472363	127473530	Pos2	0	+	127472363	127473530	255,0,0
   chr7	127473530	127474697	Neg1	0	-	127473530	127474697	0,0,255

Track lines
-----------

A BED file can begin with a ``track`` line that tells a genome browser
how to display the features that follow:

.. code-block:: text

   track name="ItemRGBDemo" description="Item RGB demonstration" itemRgb="On"
   chr7	127471196	127472363	Pos1	0	+	127471196	127472363	255,0,0
   chr7	127472363	127473530	Pos2	0	+	127472363	127473530	255,0,0

BedGraph
--------

A closely related format, **bedGraph**, is used to display continuous
data (such as coverage) along the genome. Each line gives a region and a
single value:

.. code-block:: text

   track type=bedGraph name="BedGraph Format" description="BedGraph format"
   chr19	49302000	49302300	-1.0
   chr19	49302300	49302600	-0.75
   chr19	49302600	49302900	-0.50

Software that use BED format
----------------------------

BED is understood by nearly every genome browser and interval tool. A
very common one is `BEDTools
<https://bedtools.readthedocs.io/en/latest/>`_, a toolkit for comparing,
merging, and manipulating genomic intervals.
