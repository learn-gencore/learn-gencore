# Next-Generation Sequencing Analysis Resources

Resources for mastering NGS analysis, maintained by the Bioinformatics team at the NYU Center for Genomics and Systems Biology in [New York](http://gencore.bio.nyu.edu) and [Abu Dhabi](https://cgsb.abudhabi.nyu.edu).

Each module provides a step-by-step guide to using standard analysis pipelines with state-of-the-art bioinformatics tools. Modules include sample datasets and scripts that can be accessed on NYU's HPC facility.

**Published site:** https://learn.gencore.bio.nyu.edu/

## Topics covered

- **Transcriptomics:** bulk RNA-seq, single-cell RNA-seq, multiome, gene set enrichment analysis
- **Variant detection:** whole-exome sequencing (WES), whole-genome sequencing (WGS)
- **Epigenomics:** ATAC-seq, ChIP-seq, scATAC
- **De novo assembly:** genome, transcriptome, bacterial assembly
- **Metagenomics:** shotgun, amplicon (16S)
- **Primary analysis:** demultiplexing (Pheniqs)
- **Visualization:** IGV
- **Programming:** Python, R, shell, SLURM

## Contributing

We use [Sphinx](https://www.sphinx-doc.org/) with reStructuredText (`.rst`). While Markdown is simpler, Sphinx provides better structure, supports cross-referencing, and scales more effectively for larger and more technical documentation.

### Authentication setup

We use SSH for authentication between your laptop and GitHub. Set up SSH keys and verify the connection is working:
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account

### Clone the repository

```bash
git clone git@github.com:learn-gencore/learn-gencore.git
cd learn-gencore/
```

### Create a new branch

```bash
git checkout main
git pull origin main
git checkout -b <branch-name>
git branch  # verify your new branch is created and active
```

### Make changes

Navigate to the relevant directory under `source/` and edit or add `.rst` files.

```bash
cd source/<path-to-your-topic>
```

When you are done, stage and commit:

```bash
git status
git add .
git status  # verify modified files
git commit -m "Add a meaningful commit message"
```

### Push and create a pull request

```bash
git branch   # verify the branch name
git push origin <branch-name>
```

Then go to https://github.com/learn-gencore/learn-gencore/. You will see a prompt like "Compare & pull request" for your recently pushed branch. Click it, add a clear title and description, and click "Create pull request".

### Pull request approval

- Any team member can approve pull requests. There is no need to wait for a specific person's approval.
- Users cannot approve their own pull requests (e.g., if User A creates a PR, User B must approve it).
- Please review pull requests carefully before approving.

**Approval procedure:**

1. Open the pull request.
2. Click the **Files changed** tab and review the changes.
3. Click **Submit review** and select **Approve**.
4. Select **Squash and merge** from the merge options dropdown.
5. Click **Confirm squash and merge**.

### Clean up after merge

```bash
git checkout main
git pull origin main
git branch -D <branch-name>
git branch  # verify deletion
```

For every new change, follow the same process from the beginning. Do not reuse old branches for new changes.

## GitHub notifications

To receive notifications for new pull requests, enable watch settings in GitHub:

**Watch** > **Custom** > **Pull Requests** > **Apply**

## RST formatting reference

```rst
Title
=====

Subheading
----------

Bold -> **word**
Italics -> *word*

Bullets:
* Item 1
* Item 2

Hyperlink: `Text <https://example.com>`_
```

Heading underlines must be at least as long as the heading text.

For more details: https://sphinx-intro-tutorial.readthedocs.io/en/latest/rst_intro.html
