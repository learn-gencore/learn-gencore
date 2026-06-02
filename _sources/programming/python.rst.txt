Python
======

This tutorial is a hands-on introduction to data analysis in Python
for people who are new to the language. It walks through a typical
workflow on a real biology dataset: bulk RNA-seq gene expression
measured across eight tissues in human and mouse. By the end you
will have loaded a table from disk, cleaned it up, made some plots,
combined two datasets, and saved the result back out.

The material is adapted from the *Applied data analysis for
beginners* module of the NYUAD Core Bioinformatics Workshops
(*Start Analytics 2025*). Source data and notebooks are at
`NYUAD-Core-Bioinformatics/StartAnalytics_2025
<https://github.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/>`_.

No prior Python experience is assumed. If you already know Python,
skim the "What you need to know first" boxes and read the code.

What you need to know first
---------------------------

A few terms that come up over and over:

* **Library** (also called a "package" or "module") — a chunk of
  pre-written code you can use. We use four:

  * ``pandas`` — tables of data (think Excel sheets in code).
  * ``numpy`` — fast math on numbers, arrays, and matrices.
  * ``matplotlib`` — the standard plotting library.
  * ``seaborn`` — nicer-looking plots built on top of matplotlib.

* **DataFrame** — a table with named rows and named columns. This
  is what ``pandas`` gives you. Each column can hold a different
  type (numbers, text, dates).

* **Index** — the labels on the rows of a DataFrame. In our case
  the index is the gene name (e.g. ``STAG2``); the column names
  are the sample names (e.g. ``human_brain_4``).

* **Axis** — ``axis=0`` means "go down columns" (operate on rows).
  ``axis=1`` means "go across rows" (operate on columns). It is
  the opposite of what most people guess on the first try.

Workflow overview
-----------------

A typical analysis on a gene-expression table follows six steps:

1. **Load** the data into a DataFrame.
2. **Filter** out empty samples and unexpressed genes.
3. **Slice** the table to focus on samples or genes of interest.
4. **Transform** the numbers (normalize, log-transform) so plots
   and comparisons are meaningful.
5. **Merge** datasets (e.g. human + mouse) and check that the
   biology makes sense via a correlation heatmap.
6. **Export** the cleaned table to CSV for downstream work.

Data shape
~~~~~~~~~~

The data table looks like this:

* **rows** — genes (one per row, e.g. ``STAG2``, ``GOSR2``, ...).
* **columns** — samples / experiments
  (e.g. ``human_brain_4``, ``mouse_liver_5``).
* **values** — raw "read counts": how many sequencing reads from
  that sample mapped to that gene. A higher number means the gene
  is more strongly expressed in that sample.

Before any analysis, ask yourself: *Are there empty samples? Do all
genes have a value? Which features are irrelevant? Are there
outliers or missing data?* Garbage in, garbage out.

Environment setup
-----------------

Install the libraries once from a terminal:

.. code:: bash

    pip install numpy pandas matplotlib seaborn

Then, at the top of your Python script or notebook, import them.
The ``as np`` / ``as pd`` part gives each library a short nickname
so you do not have to type the full name every time.

.. code:: python

    import numpy as np                 # numerical math
    import pandas as pd                # tables (DataFrames)
    import matplotlib.pyplot as plt    # plotting
    import seaborn as sns              # nicer plots

    import warnings
    warnings.filterwarnings("ignore")  # hide non-critical warnings

Inputs
------

The example uses three CSV files (provided as gzipped CSVs in the
workshop repo). Python's ``pandas`` reads ``.csv.gz`` natively — no
need to unzip them first.

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
  each human gene symbol (uppercase, e.g. ``STAG2``) with its
  mouse equivalent (capitalized, e.g. ``Stag2``).

Or grab all three from the terminal:

.. code:: bash

    mkdir -p data
    BASE=https://raw.githubusercontent.com/NYUAD-Core-Bioinformatics/StartAnalytics_2025/main/notebooks/data
    curl -L -o data/human_expression.csv.gz       $BASE/human_expression.csv.gz
    curl -L -o data/mouse_expression.csv.gz       $BASE/mouse_expression.csv.gz
    curl -L -o data/human_mouse_gene_map.csv.gz   $BASE/human_mouse_gene_map.csv.gz

1. Loading data
---------------

``pd.read_csv`` reads a CSV file into a DataFrame. The
``index_col=0`` argument says "use the first column of the file as
the row labels" — so gene names become the index rather than a
plain data column.

.. code:: python

    mouse    = pd.read_csv('data/mouse_expression.csv.gz', index_col=0)
    human    = pd.read_csv('data/human_expression.csv.gz', index_col=0)
    gene_map = pd.read_csv('data/human_mouse_gene_map.csv.gz', index_col=0)

A few quick ways to look at what just got loaded:

.. code:: python

    human.shape       # (rows, cols) tuple — here (14744, 8)
    human.head()      # first 5 rows
    human.columns     # the column (sample) names
    human.index       # the row (gene) names

Typing the name of a DataFrame on its own line in a notebook
prints a preview of it. In a regular script use ``print(human)``.

2. Filtering
------------

Filtering means dropping rows or columns that carry no useful
information.

Empty samples
~~~~~~~~~~~~~

``.sum()`` adds up the values. By default it sums each column,
giving the total number of reads per sample. Any sample with a
total of 0 sequenced badly and should be removed.

.. code:: python

    human.sum()        # one total per column

In this dataset every sample has millions of reads, so nothing to
drop here.

Genes with no expression
~~~~~~~~~~~~~~~~~~~~~~~~

To check the *other* direction — genes that are zero in every
sample — pass ``axis=1`` so the sum runs across each row instead
of each column:

.. code:: python

    human.sum(axis=1).describe()   # summary stats: min, max, mean, etc.

``.describe()`` prints quick summary statistics for whatever it is
called on. Useful for sanity-checking distributions.

There are 674 genes with zero counts in every tissue. We can see
them with a **boolean filter** — ``human.sum(1) == 0`` produces a
Series of True/False, and using it inside ``human[...]`` keeps
only the rows where it is True:

.. code:: python

    empty_genes = human[human.sum(1) == 0]
    empty_genes.shape      # (674, 8)

Drop those rows and keep only genes expressed in at least one
tissue. ``.gt(0)`` is the same as ``> 0`` in this context:

.. code:: python

    human = human[human.sum(1).gt(0)]
    human.shape            # (14070, 8)

.. note::

   ``human.sum(1).gt(0)`` and ``~human.sum(1).eq(0)`` (read: "not
   equal to zero") do the same thing. ``~`` is the "not" operator
   for boolean Series.

3. Slicing
----------

Slicing means picking a piece of the table — particular rows,
particular columns, or both. Square brackets ``[]`` are the main
tool.

Pick columns whose name contains a word
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The expression ``[c for c in human.columns if "brain" in c]`` is a
**list comprehension**. Read it left-to-right: "make a list of
``c`` for each ``c`` in ``human.columns``, but only when the text
``brain`` appears in ``c``."

.. code:: python

    brain_cols = [c for c in human.columns if "brain" in c]
    brain_cols       # ['human_brain_4']
    human[brain_cols]

Pick rows by gene name
~~~~~~~~~~~~~~~~~~~~~~

Same trick on the row labels (``human.index``). The *DLX* genes
are a family of homeobox transcription factors active in the
brain.

.. code:: python

    dlx_rows = [g for g in human.index if "DLX" in g]
    human.loc[dlx_rows]

``human.loc[...]`` selects rows by their *label*. Plain
``human[...]`` selects columns. They look similar but they target
different axes, which is a common point of confusion.

Pick rows and columns at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``.loc`` accepts both row and column lists:

.. code:: python

    subset = human.loc[dlx_rows, brain_cols]
    n_rows, n_cols = subset.shape
    print(f"subset has: {n_rows} rows, and {n_cols} cols")
    subset

The ``f"..."`` thing is an **f-string**: Python plugs the value of
``n_rows`` and ``n_cols`` into the string for you.

.. note::

   ``.loc[rows, cols]`` selects by label.
   ``.iloc[rows, cols]`` selects by integer position
   (``.iloc[0:5]`` = first five rows). Mixing them is a common
   source of off-by-one bugs — pick one and stick with it per
   call.

4. Transforming data
--------------------

Raw counts are hard to interpret directly for two reasons:

1. A handful of highly-expressed genes are so much bigger than
   everything else that they dominate any plot or statistic
   (outliers).
2. Different samples were sequenced to different depths, so a gene
   can look "higher" in one sample purely because that sample was
   sequenced more deeply, not because of biology.

We fix both with two standard transformations:

* **Normalize by library size** (counts per million, CPM) — scale
  each sample so all samples sum to the same total, making them
  comparable.
* **Log-transform** (``log1p(x)`` = ``log(1 + x)``) — compresses
  the long tail of highly-expressed genes. We use ``log1p``
  instead of plain ``log`` so that zero counts stay defined
  (``log(0)`` is undefined; ``log1p(0) = 0``).

A quick tour of matplotlib and seaborn
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before the first plot, a short mental model.

**matplotlib** is the foundational plotting library. Two ways to
use it:

* **Pandas/seaborn shortcuts** — many objects know how to plot
  themselves. ``human.plot.box(...)`` or
  ``sns.violinplot(data=human)`` build a matplotlib figure for
  you and return the ``Axes`` (one panel) so you can customize
  it.
* **Direct matplotlib** — ``fig, axes = plt.subplots(...)``
  creates an explicit figure and one or more ``Axes`` panels;
  you then pass ``ax=...`` to the shortcut to draw into a
  specific panel.

Useful pieces you will see below:

* ``plt.title("...")`` / ``plt.xlabel(...)`` / ``plt.ylabel(...)``
  — labels on the current panel.
* ``plt.xticks(rotation=90)`` — tilt tick labels so long sample
  names fit.
* ``plt.tight_layout()`` — auto-fix overlapping titles/labels
  when you have multiple panels.
* ``plt.show()`` — display the figure (only needed in plain
  scripts; notebooks show it automatically).

**seaborn** is built on top of matplotlib and adds two things we
use:

* Nicer defaults and palettes (``palette="pastel"``,
  ``cmap="coolwarm"``).
* Higher-level plot types that work directly on DataFrames:
  ``sns.boxplot``, ``sns.violinplot``, ``sns.heatmap``,
  ``sns.clustermap``.

One practical note: when you mix pandas's ``.plot.box()`` with
seaborn calls and ``plt.subplots``, all three are drawing onto
the *same* matplotlib figure underneath. That is why a
``plt.title(...)`` call right after a ``.plot.box(...)`` line
still affects the right plot.

Boxplot of raw counts
~~~~~~~~~~~~~~~~~~~~~

A boxplot summarizes each sample's distribution: the box is the
middle 50% of values, the whiskers extend outward, and outliers
appear as individual dots.

.. code:: python

    colors = [
        "lightblue", "lightgreen", "lightpink", "lightsalmon",
        "lightyellow", "powderblue", "lavender", "thistle",
    ]

    ax = human.plot.box(figsize=(8, 6), patch_artist=True)
    for patch, color in zip(ax.artists, colors):
        patch.set_facecolor(color)
    plt.title("Box plot of expression by samples", fontsize=14)
    plt.xlabel("Tissue sample", fontsize=12)
    plt.ylabel("Expression", fontsize=12)
    plt.xticks(rotation=90)   # tilt sample names so they fit
    plt.show()

The raw plot is unreadable: a few enormous outlier genes squash
the rest of the distribution onto the baseline.

Log1p transform
~~~~~~~~~~~~~~~

``np.log1p(human)`` returns a new DataFrame where every value has
been replaced with ``log(1 + value)``. The original ``human``
table is *not* modified.

.. code:: python

    ax = np.log1p(human).plot.box(figsize=(8, 6), patch_artist=True)
    for i, patch in enumerate(ax.patches[:len(human.columns)]):
        patch.set_facecolor(colors[i])
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)
    plt.title("Box plot of log1p expression by samples", fontsize=14)
    plt.xlabel("Tissue sample", fontsize=12)
    plt.ylabel("log1p Expression", fontsize=12)
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

Now the distributions are visible. But samples still differ in
overall level because of library size — that is what the next
step fixes.

Normalize by library size (CPM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Divide every value in a column by that column's total (so each
column now sums to 1), then multiply by one million for readable
numbers:

.. code:: python

    normed_human = human.div(human.sum(axis=0), axis=1) * 1e6
    normed_human.loc[dlx_rows]

* ``human.sum(axis=0)`` gives one total per column.
* ``human.div(..., axis=1)`` divides each row of ``human`` by
  those totals, column-by-column.
* ``* 1e6`` rescales to counts per million.

Visually compare all three views — raw, normalized, log-normalized
— side by side. ``plt.subplots(1, 3)`` makes one figure with
three panels.

.. code:: python

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    titles = ['Raw expression', 'Normed', 'l1p Normed']
    for i, (df, ax) in enumerate(
        zip([human, normed_human, np.log1p(normed_human)], axes)
    ):
        boxes = df.boxplot(
            ax=ax, patch_artist=True,
            return_type="dict", boxprops=dict(color="black"),
        )
        for patch, color in zip(boxes['boxes'], colors):
            patch.set_facecolor(color)
        ax.set_title(titles[i])
        if i == 0:
            ax.set_ylabel("Expression")
        ax.set_xticklabels(df.columns, rotation=90)
    plt.tight_layout()
    plt.show()

Violin plots show the full distribution shape, not just the
quartiles:

.. code:: python

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for i, (df, ax) in enumerate(
        zip([human, normed_human, np.log1p(normed_human)], axes)
    ):
        sns.violinplot(data=df, ax=ax, palette="pastel", inner="quartile")
        ax.set_title(titles[i])
        ax.set_xticklabels(df.columns, rotation=90)
    plt.tight_layout()
    plt.show()

Fold change before and after transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Fold change** asks: how much higher (or lower) is gene X in
sample B compared to sample A? The simple formula is
``(B - A) / A``: a value of ``1`` means twice as high, ``-0.5``
means half as high, ``0`` means equal.

Computing fold change on raw counts gives implausibly large values
driven by sequencing depth. Computing it on log1p-normalized
values gives interpretable numbers.

A small helper function. Defining a function avoids repeating the
same lines three times:

.. code:: python

    def fold_change(df, ref, rows):
        a = df.loc[rows, ref]                       # reference column
        cols = df.columns[df.columns != ref]        # everything else
        fc = {c: (df.loc[rows, c] - a) / a for c in cols}
        fc = pd.DataFrame(fc).fillna(0)             # NaN -> 0
        fc = fc.mask(np.isinf(fc), df.loc[rows, cols])  # inf -> raw val
        return fc

    fc_raw        = fold_change(human,                 "human_brain_4", dlx_rows)
    fc_normed     = fold_change(normed_human,          "human_brain_4", dlx_rows)
    fc_l1p_normed = fold_change(np.log1p(normed_human),"human_brain_4", dlx_rows)

``.fillna(0)`` replaces missing values (NaN, from ``0 / 0``) with
zero. ``.mask(condition, other)`` replaces values where
``condition`` is True with values from ``other`` — here we
replace infinities (from dividing by zero) with the raw value, as
a fallback.

Heatmap comparison — colour-coded so the difference is obvious:

.. code:: python

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    sns.heatmap(fc_raw,        annot=True, fmt=".2f",
                cmap="coolwarm", ax=axes[0], cbar=False)
    axes[0].set_title("Raw expression FC")
    sns.heatmap(fc_normed,     annot=True, fmt=".2f",
                cmap="coolwarm", ax=axes[1], cbar=False)
    axes[1].set_title("Normed FC")
    sns.heatmap(fc_l1p_normed, annot=True, fmt=".2f",
                cmap="coolwarm", ax=axes[2], cbar=True)
    axes[2].set_title("l1p Normed FC")
    plt.tight_layout()
    plt.show()

Raw fold-changes contain values like ``19`` or ``-1`` driven
entirely by library size; the log1p-normalized panel shows
interpretable values in roughly ``[-1, 1]``.

5. Merging data
---------------

To compare expression across species, we need to align genes:
``STAG2`` in human is the same gene as ``Stag2`` in mouse. The
orthologue map provides that link.

.. code:: python

    gene_map.head()
    #                          human          mouse
    # Stag2_STAG2              STAG2          Stag2
    # Stag1_STAG1              STAG1          Stag1
    # ...

Process the mouse table the same way (normalize + log1p), then
merge each species' table onto the gene map. ``merge`` is
``pandas``'s version of a SQL join — it lines up rows by a key
column.

.. code:: python

    normed_mouse = np.log1p(mouse.div(mouse.sum(axis=0), axis=1) * 1e6)

    star_normed_mouse = normed_mouse.merge(
        gene_map, left_index=True, right_on='mouse',
    )
    star_normed_human = np.log1p(normed_human).merge(
        gene_map, left_index=True, right_on='human',
    )

* ``left_index=True`` — use the left table's row index as the
  join key.
* ``right_on='mouse'`` — match against the ``mouse`` column on
  the right table.

After both merges, the two tables share the same row labels
(``Stag2_STAG2``, etc.). Join them and drop the duplicate
``human``/``mouse`` columns introduced by the second merge
(``suffixes=('', '_drop')`` tags the duplicates so we can find
and remove them):

.. code:: python

    merged_df = star_normed_mouse.merge(
        star_normed_human,
        right_index=True, left_index=True,
        suffixes=('', '_drop'),
    )
    merged_df = merged_df.loc[
        :, ~merged_df.columns.str.endswith('_drop')
    ]
    merged_df.shape    # (14744, 18)

Sample correlation
~~~~~~~~~~~~~~~~~~

A correlation heatmap is a quick sanity check: matched tissues
(human brain ↔ mouse brain) should correlate more strongly than
mismatched tissues, even across species.

``.corr()`` computes the pairwise Pearson correlation between
every pair of columns.

.. code:: python

    correlation_matrix = merged_df.loc[
        :, ~merged_df.columns.isin(['human', 'mouse'])
    ].corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix,
                annot=True, cmap='coolwarm',
                fmt=".2f", linewidths=0.5)
    plt.title('Correlation Heatmap')
    plt.show()

    sns.clustermap(correlation_matrix,
                   annot=True, cmap='coolwarm',
                   fmt=".2f", linewidths=0.5)
    plt.show()

``sns.clustermap`` reorders the rows and columns so similar
samples sit next to each other, making tissue clusters jump out
even when the original ordering was arbitrary.

6. Exporting data
-----------------

Save the cleaned, transformed, merged table to disk so the next
step in your pipeline (differential expression, enrichment
analysis, ...) can read it back.

.. code:: python

    merged_df.to_csv('./data/human_mouse_merged_expression.csv')

Useful options:

* ``index=False`` — do not write the row labels (use this if you
  do not want a leading unnamed column when you re-read the file).
* Pass a ``.csv.gz`` filename and ``pandas`` writes it gzipped for
  free.

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
