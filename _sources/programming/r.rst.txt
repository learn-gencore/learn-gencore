R / RShiny
==========

This tutorial is a hands-on introduction to data analysis in R for
people who are new to the language. It follows the same workflow as
the :doc:`python` tutorial — load, filter, slice, transform, merge,
save — using the same bulk RNA-seq dataset (gene expression across
eight tissues in human and mouse).

The material is adapted from the *Applied data analysis for
beginners* module of the NYUAD Core Bioinformatics Workshops
(*Start Analytics 2025*). Source data and notebooks are at
`NYUAD-Core-Bioinformatics/StartAnalytics_2025
<https://github.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/>`_.

No prior R experience is assumed.

What you need to know first
---------------------------

A few terms that come up over and over:

* **Package** — a chunk of pre-written code you can use. We use:

  * ``dplyr`` — verbs for data frames (``filter``, ``select``,
    ``mutate``, ``left_join``) and the pipe ``%>%``.
  * ``tidyr`` — reshape data between wide and long form
    (``pivot_longer``).
  * ``tibble`` — helpers like ``rownames_to_column``.
  * ``reshape2`` — older reshaper, provides ``melt``.
  * ``ggplot2`` — the standard plotting library.
  * ``gridExtra`` — arrange several ggplot panels in a grid
    (``grid.arrange``).
  * ``pheatmap`` — pretty heatmaps with clustering.

* **data.frame** — a table with named rows and named columns. The
  R equivalent of a ``pandas`` DataFrame. After we read a CSV
  with ``row.names = 1`` the first column of the file becomes the
  row labels.

* **matrix** — a 2-D grid of values that are *all the same type*
  (usually all numeric). Most of our math (sums, normalization,
  correlation) is done on a matrix.

* **The pipe** ``%>%`` — takes whatever is on the left and feeds
  it as the first argument to the function on the right.
  ``x %>% f(y)`` is the same as ``f(x, y)``. It lets you read
  transformations left-to-right instead of inside-out.

* **Row and column ops** — ``rowSums(x)`` sums each row,
  ``colSums(x)`` sums each column. In Python you would write
  ``df.sum(axis=1)`` and ``df.sum(axis=0)``; in R the function
  name itself tells you the direction.

* **The transpose trick** — ``t(t(x) / colSums(x))`` divides each
  *column* of ``x`` by its column total. Without the two ``t()``
  calls, R would try to divide row-wise and recycle the values
  the wrong way. (More on this in the normalization step.)

Workflow overview
-----------------

1. **Load** the data into a data frame.
2. **Filter** — inspect samples and genes; spot anything empty.
3. **Slice** the table to focus on samples or genes of interest.
4. **Transform** the numbers (normalize, log-transform).
5. **Merge** human and mouse data and check correlation.
6. **Save** results to disk.

Data shape
~~~~~~~~~~

* **rows** — genes (e.g. ``STAG2``, ``GOSR2``).
* **columns** — samples (e.g. ``human_brain_4``).
* **values** — raw sequencing read counts.

Before any analysis, ask yourself: *Are there empty samples? Do
all genes have a value? Which features are irrelevant? Are there
outliers or missing data?* Garbage in, garbage out.

Environment setup
-----------------

Install the packages once from the R console:

.. code:: r

    install.packages(c("dplyr", "tidyr", "tibble", "reshape2",
                       "ggplot2", "gridExtra", "pheatmap"))

Load them at the top of each script or notebook session:

.. code:: r

    library(dplyr)
    library(gridExtra)
    library(reshape2)
    library(ggplot2)
    library(tidyr)
    library(tibble)
    library(pheatmap)

If you are working inside a Jupyter notebook with the R kernel,
this option keeps the printed tables short so they fit on screen:

.. code:: r

    options(repr.matrix.max.rows = 10, repr.matrix.max.cols = 6)

Inputs
------

The example uses three gzipped CSVs (``read.csv`` handles ``.gz``
automatically — no need to unzip):

Download the three files from
`NYUAD-Core-Bioinformatics/StartAnalytics_2025 / notebooks / data
<https://github.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/tree/main/notebooks/data>`_
and place them in a ``data/`` folder next to your script or
notebook:

* ``data/human_expression.csv.gz`` — raw read counts, 14,744 genes
  × 8 human tissue samples.
* ``data/mouse_expression.csv.gz`` — raw read counts, 14,744 genes
  × 8 mouse tissue samples.
* ``data/human_mouse_gene_map.csv.gz`` — a lookup table pairing
  each human gene symbol (``STAG2``) with its mouse equivalent
  (``Stag2``).

Or grab all three from the terminal:

.. code:: bash

    mkdir -p data
    BASE=https://raw.githubusercontent.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/main/notebooks/data
    curl -L -o data/human_expression.csv.gz       $BASE/human_expression.csv.gz
    curl -L -o data/mouse_expression.csv.gz       $BASE/mouse_expression.csv.gz
    curl -L -o data/human_mouse_gene_map.csv.gz   $BASE/human_mouse_gene_map.csv.gz

1. Loading data
---------------

``read.csv`` reads a CSV file into a data frame. The argument
``row.names = 1`` says "use the first column of the file as the
row labels" — so gene names become the row names rather than a
plain data column.

.. code:: r

    mouse    <- read.csv('data/mouse_expression.csv.gz', row.names = 1)
    human    <- read.csv('data/human_expression.csv.gz', row.names = 1)
    gene_map <- read.csv('data/human_mouse_gene_map.csv.gz', row.names = 1)

A few quick ways to look at what just got loaded:

.. code:: r

    dim(human)       # 14744    8     (rows, cols)
    head(human)      # first 6 rows
    colnames(human)  # sample names
    rownames(human)[1:5]   # first 5 gene names

In a notebook, just typing ``human`` and running the cell prints
the table. In a plain script use ``print(human)``.

.. note::

   ``View(human)`` (capital V) opens the table in a
   spreadsheet-style viewer — but this only works inside
   **RStudio**. It will not do anything useful from a plain R
   console, a terminal session, or most Jupyter setups.

2. Filtering
------------

Filtering means dropping rows or columns that carry no useful
information. At this stage we just *inspect* — the actual drop
happens later if needed.

Empty samples
~~~~~~~~~~~~~

``colSums`` adds up the values column-by-column, giving total
reads per sample:

.. code:: r

    colSums(human)

In this dataset every sample has millions of reads, so no empty
columns to drop.

Genes with no expression
~~~~~~~~~~~~~~~~~~~~~~~~

``rowSums`` adds up across each row — the total reads per gene
across all samples:

.. code:: r

    gene_count_across_tissues <- rowSums(human)
    head(gene_count_across_tissues)
    summary(gene_count_across_tissues)

``summary()`` is R's quick descriptive-stats function (min, 1st
quartile, median, mean, 3rd quartile, max). It is the closest
analogue to pandas' ``.describe()``.

How many genes have zero expression in every tissue?
``gene_count_across_tissues == 0`` returns a logical vector
(``TRUE``/``FALSE``), and ``summary()`` on a logical vector
prints a count of each value:

.. code:: r

    summary(gene_count_across_tissues == 0)
    # Mode   FALSE    TRUE
    # logical 14070     674

So 674 genes are flat zero across all tissues; 14,070 are
expressed in at least one tissue.

3. Slicing
----------

Slicing means picking a subset — particular rows, particular
columns, or both. The base-R bracket form is ``x[rows, cols]``;
leave a side blank to keep "all".

Pick columns whose name contains a word
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``grep("brain", colnames(human))`` returns the *positions*
(integers) of column names that match. Wrap it in
``colnames(human)[...]`` to get the names themselves back:

.. code:: r

    brain_cols <- colnames(human)[grep("brain", colnames(human))]
    brain_cols                  # "human_brain_4"

Pick rows by gene name
~~~~~~~~~~~~~~~~~~~~~~

Same idea applied to ``rownames``. The *DLX* genes are a family
of homeobox transcription factors active in the brain.

.. code:: r

    dlx_rows <- rownames(human)[grep("DLX", rownames(human))]
    dlx_rows

Extract the subset of data
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: r

    human[dlx_rows, brain_cols]

The order is always ``[rows, cols]``. Forgetting the comma is a
common mistake — ``human[dlx_rows]`` would try to interpret
``dlx_rows`` as column names and break.

4. Transform data
-----------------

Raw counts are hard to interpret directly for two reasons:

1. A handful of highly-expressed genes dominate any plot or
   statistic (outliers).
2. Different samples were sequenced to different depths, so a
   higher count can reflect deeper sequencing rather than higher
   biological expression.

We apply two standard transformations:

* **Normalize by library size** (counts per million, CPM) — scale
  each sample so column totals match, making samples comparable.
* **Log-transform** (``log1p(x)`` = ``log(1 + x)``) — compresses
  the long tail of highly-expressed genes. ``log1p(0)`` is
  ``0``, so zero counts stay defined.

Wide vs long format
~~~~~~~~~~~~~~~~~~~

Before the first plot, one concept that trips up everyone new to
``ggplot2``: the difference between **wide** and **long** format.

Our expression table is in **wide** format — one column per
sample, one row per gene. It is convenient to read and to do
matrix math on:

============  ===============  ===============  ===============
Gene          human_brain_4    human_heart_2    human_liver_2
============  ===============  ===============  ===============
STAG2                    2513             1734             1262
GOSR2                      40              187               87
DLX6                      142                0                0
============  ===============  ===============  ===============

The same data in **long** format has *one row per observation* —
i.e. one row per (gene, sample) pair — with the sample name
moved into its own column:

============  =================  =============
Gene          Sample             Expression
============  =================  =============
STAG2         human_brain_4               2513
STAG2         human_heart_2               1734
STAG2         human_liver_2               1262
GOSR2         human_brain_4                 40
GOSR2         human_heart_2                187
GOSR2         human_liver_2                 87
DLX6          human_brain_4                142
DLX6          human_heart_2                  0
DLX6          human_liver_2                  0
============  =================  =============

Why do we need long format? Because ``ggplot2``'s ``aes()`` maps
*columns* of data to visual properties. To draw "one box per
sample, coloured by sample," ggplot needs a single ``Sample``
column to point ``x =`` and ``fill =`` at. A wide table has no
such column — the sample identity is hidden in the *names* of
the columns.

The function that converts wide → long is ``tidyr::pivot_longer``:

.. code:: r

    long <- as.data.frame(human) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene,                # collapse everything except Gene
                   names_to = "Sample",  # column names go here
                   values_to = "Expression")  # values go here

* ``-Gene`` means "all columns except ``Gene``" — those are the
  ones that get stacked.
* ``names_to`` is the new column that will hold the old column
  names.
* ``values_to`` is the new column that will hold the old values.

You will see this exact three-line pattern many times below —
whenever a matrix needs to become a plot.

A quick tour of ggplot2
~~~~~~~~~~~~~~~~~~~~~~~

With long format sorted out, here is the ggplot mental model.
``ggplot2`` is based on the *grammar of graphics*: every plot is
built up by adding **layers** with the ``+`` operator. A minimal
plot needs three things:

1. **Data** — a data frame in long format (one row per
   observation).
2. **Aesthetics** — ``aes(...)`` maps columns of the data to
   visual properties (x, y, colour, fill, size, ...).
3. **Geometry** — a ``geom_*()`` layer that says *how* to draw
   (``geom_point``, ``geom_boxplot``, ``geom_violin``,
   ``geom_tile``, ...).

A complete plot reads like a sentence:

.. code:: r

    ggplot(data, aes(x = sample, y = expression)) +   # data + mapping
        geom_boxplot() +                               # how to draw
        labs(title = "...", x = "...", y = "...") +    # labels
        theme_bw()                                     # overall look

Two practical notes:

* **Long format is required.** Each row must be one observation.
  A wide table with one column per sample will not plot directly
  — reshape it first with ``pivot_longer``.
* **The plus signs matter.** ``+`` belongs at the *end* of the
  previous line; if you put it at the start of the next line, R
  thinks the previous line was complete and runs it on its own.

Boxplot of raw counts
~~~~~~~~~~~~~~~~~~~~~

``ggplot2`` needs long format, so we pivot the wide matrix first.
``rownames_to_column("Gene")`` moves the row labels into a
regular column so ``pivot_longer`` can leave it alone.

.. code:: r

    human_long <- as.data.frame(human) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

    boxplot_human <- ggplot(human_long,
                            aes(x = Sample, y = Expression, fill = Sample)) +
      geom_boxplot(outlier.shape = 21, outlier.colour = "black",
                   outlier.fill = "white", color = "black") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Box Plot of expression by samples",
           x = "Tissue sample", y = "Expression") +
      theme(legend.position = "none")

    boxplot_human

The raw plot is unreadable: a few enormous outlier genes squash
the rest of the distribution onto the baseline.

Log1p transform
~~~~~~~~~~~~~~~

``log1p()`` applied to a data frame returns a new data frame
where every value has been replaced by ``log(1 + value)``. The
original is not modified.

.. code:: r

    log1p_human <- log1p(human)

    log1p_human_long <- as.data.frame(log1p_human) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

    sample_colors <- c(
        "lightblue", "lightgreen", "lightpink", "lightsalmon",
        "lightyellow", "powderblue", "lavender", "thistle"
    )

    ggplot(log1p_human_long,
           aes(x = Sample, y = Expression, fill = Sample)) +
      geom_boxplot(width = 0.5,
                   outlier.shape = 21, outlier.colour = "black",
                   outlier.fill = "white", color = "black") +
      scale_fill_manual(values = sample_colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Box Plot of expression by samples",
           x = "Tissue Sample", y = "log1p Expression") +
      theme(legend.position = "none")

Now the distributions are visible — but samples still differ in
overall level because of library size. The next step fixes that.

Normalize by library size (CPM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normalization is one line:

.. code:: r

    normed_human <- t(t(human) / colSums(human)) * 1e6
    normed_human[dlx_rows, ]

Read it from the inside out:

1. ``colSums(human)`` — total reads per sample, a length-8
   vector.
2. ``t(human)`` — *transpose* the table so samples are rows and
   genes are columns.
3. ``t(human) / colSums(human)`` — now R divides each row of the
   transposed table by the matching element of the vector, which
   is what we want (per-sample division).
4. ``t(...)`` — transpose back so genes are rows again.
5. ``* 1e6`` — rescale to counts per million.

Without the two transposes, R would recycle the column totals
across rows in the wrong direction and produce nonsense.

Compare all three views side by side. ``grid.arrange`` from
``gridExtra`` lays out separate ggplot objects in a grid.

.. code:: r

    options(repr.plot.width = 18, repr.plot.height = 6)  # notebook only

    # raw counts
    boxplot_human <- ggplot(human_long,
                            aes(x = Sample, y = Expression, fill = Sample)) +
      geom_boxplot(width = 0.5, outlier.shape = 21,
                   outlier.colour = "black", outlier.fill = "white",
                   color = "black") +
      scale_fill_manual(values = sample_colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      labs(title = "Raw expression", x = "Tissue sample", y = "Expression")

    # normalized
    normed_human_long <- as.data.frame(normed_human) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

    boxplot_normed_human <- ggplot(normed_human_long,
                                   aes(x = Sample, y = Expression, fill = Sample)) +
      geom_boxplot(width = 0.5, outlier.shape = 21,
                   outlier.colour = "black", outlier.fill = "white",
                   color = "black") +
      scale_fill_manual(values = sample_colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      labs(title = "Normed", x = "Tissue sample", y = "Expression")

    # log1p of normalized
    log1p_normed_human <- log1p(normed_human)
    log1p_normed_human_long <- as.data.frame(log1p_normed_human) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

    boxplot_log1p_normed_human <- ggplot(log1p_normed_human_long,
                                         aes(x = Sample, y = Expression, fill = Sample)) +
      geom_boxplot(width = 0.5, outlier.shape = 21,
                   outlier.colour = "black", outlier.fill = "white",
                   color = "black") +
      scale_fill_manual(values = sample_colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      labs(title = "l1p Normed", x = "Tissue sample", y = "Expression")

    grid.arrange(boxplot_human, boxplot_normed_human,
                 boxplot_log1p_normed_human, ncol = 3)

Violin plots show the full distribution shape. ``geom_violin`` +
``geom_boxplot`` on top of each other is a useful combination —
the violin gives density, the inner box gives the quartiles.

.. code:: r

    violinplot_human <- ggplot(human_long,
                               aes(x = Sample, y = Expression, fill = Sample)) +
      geom_violin() +
      geom_boxplot(width = 0.5, outlier.shape = 21,
                   outlier.colour = "black", outlier.fill = "white",
                   color = "black") +
      scale_fill_manual(values = sample_colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      labs(title = "Raw expression (violin)",
           x = "Tissue sample", y = "Expression")

    # build violinplot_normed_human and violinplot_log1p_normed_human
    # the same way, swapping in normed_human_long / log1p_normed_human_long

    options(repr.plot.width = 18, repr.plot.height = 12)
    grid.arrange(boxplot_human, boxplot_normed_human,
                 boxplot_log1p_normed_human,
                 violinplot_human, violinplot_normed_human,
                 violinplot_log1p_normed_human,
                 ncol = 3, nrow = 2)

Fold change before and after transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Fold change** asks: how much higher (or lower) is gene X in
sample B compared to a reference sample A? The simple formula is
``(B - A) / A``: ``1`` means twice as high, ``-0.5`` means half,
``0`` means equal.

Computing it on raw counts gives implausibly large values driven
by sequencing depth. Computing it on log1p-normalized values
gives interpretable numbers.

The pattern below is from the notebook: ``sapply`` walks over
each non-reference column, computes the fold-change against the
reference, and returns the result as a matrix.

.. code:: r

    ref <- "human_brain_4"
    a   <- normed_human[dlx_rows, ref]
    col_subset <- colnames(normed_human)[colnames(normed_human) != ref]

    normed_human_fold_changes <- sapply(col_subset, function(c) {
        b <- normed_human[dlx_rows, c]
        (b - a) / a
    })

* ``sapply(X, FUN)`` — apply ``FUN`` to each element of ``X``
  and try to simplify the result to a matrix or vector.
* ``function(c) { ... }`` — an anonymous function. ``c`` here is
  one column name from ``col_subset`` at a time.

The result will contain ``NaN`` (from ``0 / 0``) and ``Inf``
(from dividing by zero). Patch those up:

.. code:: r

    # 0 / 0 -> 0
    normed_human_fold_changes[is.na(normed_human_fold_changes)] <- 0

    # division by zero -> fall back to the raw value
    inf_mask <- is.infinite(normed_human_fold_changes)
    normed_human_fold_changes[inf_mask] <-
        human[dlx_rows, col_subset][inf_mask]

    normed_human_fold_changes <-
        as.data.frame(normed_human_fold_changes, row.names = dlx_rows)

Do the same calculation for raw counts and for log1p-normalized
counts:

.. code:: r

    # raw
    a <- human[dlx_rows, ref]
    col_subset <- colnames(human)[colnames(human) != ref]
    human_fold_changes <- sapply(col_subset, function(c) {
        b <- human[dlx_rows, c]
        (b - a) / a
    })
    human_fold_changes[is.na(human_fold_changes)] <- 0
    inf_mask <- is.infinite(human_fold_changes)
    human_fold_changes[inf_mask] <- human[dlx_rows, col_subset][inf_mask]

    # log1p of normalized
    l1p_normed_human <- log1p(normed_human)
    a <- l1p_normed_human[dlx_rows, ref]
    col_subset <- colnames(l1p_normed_human)[colnames(l1p_normed_human) != ref]
    l1p_normed_human_fold_changes <- sapply(col_subset, function(c) {
        b <- l1p_normed_human[dlx_rows, c]
        (b - a) / a
    })
    l1p_normed_human_fold_changes[is.na(l1p_normed_human_fold_changes)] <- 0
    inf_mask <- is.infinite(l1p_normed_human_fold_changes)
    l1p_normed_human_fold_changes[inf_mask] <-
        l1p_normed_human[dlx_rows, col_subset][inf_mask]

Heatmap of fold changes
~~~~~~~~~~~~~~~~~~~~~~~

Build a heatmap with ``geom_tile`` (one rectangle per cell) and
``geom_text`` (the number inside it). ``scale_fill_gradient2``
makes blue → white → red with white at zero.

.. code:: r

    human_fold_changes_long <- as.data.frame(human_fold_changes) %>%
      rownames_to_column("Gene") %>%
      mutate(Gene = factor(Gene, levels = dlx_rows)) %>%  # preserve order
      pivot_longer(-Gene, names_to = "Sample", values_to = "FoldChange")

    plot1 <- ggplot(human_fold_changes_long,
                    aes(x = Sample, y = Gene, fill = FoldChange)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white",
                           high = "red", midpoint = 0) +
      geom_text(aes(label = sprintf("%.2f", FoldChange)), size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_discrete(limits = rev) +
      labs(title = "Raw FC", x = "Sample", y = "Gene",
           fill = "Fold Change") +
      theme(legend.position = "none")

Repeat for ``normed_human_fold_changes`` and
``l1p_normed_human_fold_changes`` (call them ``plot2`` and
``plot3``) and lay them side by side:

.. code:: r

    grid.arrange(plot1, plot2, plot3, ncol = 3)

Raw fold-changes contain values like ``19`` driven entirely by
library size; the log1p-normalized panel shows interpretable
values in roughly ``[-1, 1]``.

5. Merging data
---------------

To compare expression across species, we need to align genes:
``STAG2`` in human and ``Stag2`` in mouse are the same gene. The
orthologue map provides that link.

First inspect the mouse side. Note that ``grep("DLX", toupper(...))``
matches *DLX* regardless of capitalization (the mouse symbols are
``Dlx4`` etc., not ``DLX4``):

.. code:: r

    mouse_brain_cols <- colnames(mouse)[grep("brain", colnames(mouse))]
    mouse_dlx_rows   <- rownames(mouse)[grep("DLX", toupper(rownames(mouse)))]

    mouse[mouse_dlx_rows, ]

Normalize the mouse table the same way:

.. code:: r

    normed_mouse <- log1p(t(t(mouse) / colSums(mouse)) * 1e6)
    normed_mouse[mouse_dlx_rows, ]

Filter the gene map to just the DLX orthologue pairs and look up
those rows in ``normed_mouse``:

.. code:: r

    gene_map[gene_map$mouse %in% mouse_dlx_rows, ]
    normed_mouse[gene_map[gene_map$mouse %in% mouse_dlx_rows, ]$mouse, ]

To join the gene map onto the expression tables we need the
gene-name to be a regular column (not a row name), because
``left_join`` matches by column. Build "expression + key
column" data frames for each species:

.. code:: r

    mouse_with_genes <- as.data.frame(normed_mouse)
    mouse_with_genes$mouse <- rownames(mouse_with_genes)

    normed_human_with_genes <- as.data.frame(log1p(normed_human))
    normed_human_with_genes$human <- rownames(normed_human_with_genes)

Join each species onto the gene map. ``left_join`` keeps every
row from the left table and pulls in matching rows from the right
table — unmatched rows on the right are dropped.

.. code:: r

    star_normed_mouse <- gene_map %>%
        left_join(mouse_with_genes, by = "mouse")
    rownames(star_normed_mouse) <- rownames(gene_map)

    star_normed_human <- gene_map %>%
        left_join(normed_human_with_genes, by = "human")
    rownames(star_normed_human) <- rownames(gene_map)

Now both tables have the same row order (one row per ortholog
pair). Join them and drop the duplicated key columns. After the
second join, the ``human`` column appears twice — ``left_join``
disambiguates with ``.x`` / ``.y`` suffixes, which we then drop:

.. code:: r

    merged_df <- star_normed_mouse %>%
        left_join(star_normed_human, by = "mouse") %>%
        select(-human.x, -mouse, -human.y)
    rownames(merged_df) <- rownames(star_normed_mouse)

Sample correlation
~~~~~~~~~~~~~~~~~~

A correlation heatmap is a quick sanity check: matched tissues
(human brain ↔ mouse brain) should correlate more strongly than
mismatched tissues, even across species.

``cor()`` computes pairwise Pearson correlation. The argument
``use = "pairwise.complete.obs"`` tells it to drop ``NA``\ s
per pair instead of returning ``NA`` for the whole matrix if any
``NA`` is present.

.. code:: r

    cor_matrix <- cor(merged_df,
                      use = "pairwise.complete.obs",
                      method = "pearson")
    cor_matrix

Two ways to visualize it.

**Option A — pheatmap** (with built-in clustering):

.. code:: r

    options(repr.plot.width = 12, repr.plot.height = 12)  # notebook only

    pheatmap(cor_matrix,
             display_numbers = TRUE,
             number_format = "%.2f",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_method = "complete",
             fontsize_number = 5,
             angle_col = 90,
             main = "Sample Correlation Heatmap")

The clustered version reorders rows and columns by similarity,
making tissue clusters jump out. Setting
``cluster_rows = FALSE, cluster_cols = FALSE`` keeps the original
order if you prefer.

**Option B — ggplot tile heatmap.** This needs the matrix
reshaped to long form. ``reshape2::melt`` is the classic way:

.. code:: r

    cor_long <- melt(cor_matrix)
    head(cor_long)
    # Var1            Var2     value
    # mouse_adipose_1 mouse_adipose_1 1.000
    # ...

    ggplot(cor_long, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Correlation Heatmap",
           x = "", y = "", fill = "Correlation") +
      coord_fixed() +
      scale_y_discrete(limits = rev)

If you prefer to stay inside the tidyverse, ``pivot_longer``
gives the same result. You must first turn the matrix into a
data frame so it has named columns:

.. code:: r

    cor_df <- as.data.frame(cor_matrix) %>%
        rownames_to_column(var = "Var1")

    long_data <- pivot_longer(cor_df,
                              cols = -Var1,
                              names_to = "Var2",
                              values_to = "value")

    ggplot(long_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Correlation Heatmap",
           x = "", y = "", fill = "Correlation") +
      coord_fixed() +
      scale_y_discrete(limits = rev)

Cross-species correlation only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A focused version: correlate the eight mouse samples against the
eight human samples (an 8 × 8 matrix instead of 16 × 16). Pick
the species-specific columns with ``grep`` and pattern-match on
the ``^mouse_`` / ``^human_`` prefixes.

.. code:: r

    mouse_columns <- grep("^mouse_", colnames(merged_df), value = TRUE)
    human_columns <- grep("^human_", colnames(merged_df), value = TRUE)

    mouse_data <- merged_df[, mouse_columns]
    human_data <- merged_df[, human_columns]

    cor_matrix <- cor(mouse_data, human_data,
                      use = "pairwise.complete.obs",
                      method = "pearson")
    cor_matrix

    cor_long <- melt(cor_matrix)
    ggplot(cor_long, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Mouse vs Human Correlation",
           x = "Human", y = "Mouse", fill = "Correlation") +
      coord_fixed() +
      scale_y_discrete(limits = rev)

6. Saving data
--------------

Write the merged table to disk so the next step in your pipeline
can read it back. ``write.csv`` is the base-R writer:

.. code:: r

    write.csv(merged_df, "./data/human_mouse_merged_expression.csv")

Useful options:

* ``row.names = FALSE`` — do not write the row names as a
  leading unnamed column.
* End the filename with ``.gz`` and combine with ``gzfile()``,
  or use ``readr::write_csv`` which handles ``.gz`` directly.

7. Results
-----------

The expected output files generated throughout this tutorial are
available in the course repository:

`Results directory <https://github.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/tree/main/results>`_


Next steps
----------

The merged table produced here is the standard input for the
downstream tutorials in this resource:

* :doc:`../transcriptomics/bulk_rna_seq` — differential expression.
* :doc:`../transcriptomics/over_representation_analysis` —
  GO / KEGG over-representation on a DE gene list.
* :doc:`../transcriptomics/gene_set_enrichment_analysis` — GSEA on
  a ranked gene list.
