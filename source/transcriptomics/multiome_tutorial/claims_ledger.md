# Claims ledger for `multiome.rst`

Each statement of a result on the tutorial page is listed with the computation
it rests on, where that computation lives, what it establishes, and what it does
not. Numbers come from the rendered run of `scripts/run_multivi.py` (seed 17,
377 epochs). The `results/` files are not versioned; the executed notebook
reproduces every table and the runner regenerates them with
`python scripts/run_multivi.py --reuse-model`.

Code references: `wf` = `scripts/multiome_workflow.py`, `run` =
`scripts/run_multivi.py`, `seed` = `scripts/seed_check.py`.

## Dataset and inputs

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| 3,009 called nuclei; 81,156 peaks; ARC 2.0.0; GRCh38 2020-A | `data/raw/summary.csv` ("Estimated number of cells", "ATAC Number of peaks", "Pipeline version"); fragment-file header `reference_version=2020-A` | Values reported by 10x for this output | — |
| `gex_umis_count` is zero for every barcode; exonic + intronic UMIs equal matrix totals | Check in session: 0 called barcodes with nonzero `gex_umis_count`; exonic+intronic equal `total_counts` for 3,009/3,009 nuclei. Loader comment `wf:load_multiome_dataset` | The substitution is exact for called nuclei | Why ARC leaves the column at zero |
| Fragment file about 540 MB | `ls -l data/raw/atac_fragments.tsv.gz` = 537,308,082 bytes | Size of the downloaded file | — |
| About 1.9 s per epoch on one CPU thread | Progress bar in the run log: 1.86 s/it | Speed on the rendering machine (Apple arm64, macOS) | Speed on other hardware |
| Training stalled without `OMP_NUM_THREADS=1`; importing scvi-tools (and with it PyTorch) before Scanpy crashed numba QC | Reproduced in session: synthetic MultiVI training hung at 0% CPU after Scanpy neighbors/UMAP and finished with `OMP_NUM_THREADS=1`; `calculate_qc_metrics` exited 139 when `scvi` was imported first and succeeded when Scanpy was first | Behaviour in this environment on macOS arm64 | Behaviour on Linux; the exact mechanism (three OpenMP runtimes present: conda, scikit-learn, and torch wheels) |

## 1. Cell calling

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| 14 called nuclei have fewer than 100 RNA UMIs; lowest has 29 ATAC fragments | `per_barcode_metrics.csv`, `is_cell == 1`, exonic+intronic UMIs; min `atac_fragments` | Counts at the low end of calling | Whether these are real nuclei |
| Most called nuclei sit above background on both axes | Figure 01 (`wf:plot_joint_cell_calling`) | Visual separation | A cell-calling accuracy; nuclei below ~100 UMIs or fragments overlap the background range |
| Large nucleosome-free component, mono-nucleosomal component, weak shoulder above 300 bp; ~10 bp oscillations | Figure 02 from `data/processed/fragment_length_histogram.tsv` | Shape of the pooled distribution | Per-nucleus quality |

## 2. QC thresholds

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| Threshold values in the table | `results/qc_thresholds.json`; rule `wf:QC_QUANTILES`, `wf:choose_qc_thresholds` | The data-derived values used | That they are biologically correct cutoffs |
| Rules reuse scRNA (2.5/97.5) and scATAC (1/99/99.5) tutorials; upper depth 99.5 in both modalities is this page's choice | `single_cell_rna_seq.rst` (`lower_q = 2.5`, `upper_q = 97.5`); `epigenomics/scatac_tutorial/scripts/peakvi_workflow.py` `choose_qc_thresholds` (0.01, 0.995, 0.99) | Consistency with existing pages | — |
| Median intronic fraction 0.43; ARC counts intronic UMIs | `results/joint_qc_obs.h5ad` median `intronic_fraction` = 0.431; exonic+intronic = matrix totals | Intronic UMIs are in the counts | — |
| 77 called nuclei below 500 UMIs | Per-barcode metrics, called nuclei | Size of the low-depth tail | — |
| `tss_fraction` is not TSS enrichment | Definition in `wf:load_multiome_dataset` (TSS fragments / fragments) | Naming | — |

## 3. Joint QC groups

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| 2,785 / 102 / 103 / 19 nuclei; 67 / 2 / 5 / 0 predicted doublets; 74 doublets total; 2,718 in model | `results/retention.tsv`, `results/qc_flag_counts.tsv` | Group sizes under these thresholds | Group sizes under other thresholds |
| ATAC-only medians: 320 UMIs, 32% mito, intronic 0.17 (vs 0.43 all called), FRiP 0.72 (vs 0.67), TSS 0.59 (vs 0.56), 4,431 fragments (vs 10,153); comparison group is "pass both" (2,785 nuclei, doublets included) | Medians by `qc_category` from `results/joint_qc_obs.h5ad`; notebook section 5 shows the ATAC-only medians | Descriptive contrast with the passing nuclei | Cause of the low RNA. "Suggests nuclei whose own RNA was lost or diluted" is an interpretation |
| RNA-only medians: 2,891 UMIs, 345 fragments | Same source | Descriptive | — |

## 4. Features and baselines

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| 3,000 genes, MT and MTRNR2L excluded; MTRNR2L copies track the mitochondrial fraction | `wf:select_rna_features`; `results/model_genes.tsv`; `results/depth_group_gene_mito_link.tsv` | Feature rule; the MTRNR2L–mitochondrial association | That isolation quality is the cause ("likely" on the page). That exclusion is optimal. Supporting history: a first run that kept MT genes split CD14 monocytes by MT transcripts (Wilcoxon top genes all MT-*), and a second run that kept MTRNR2L12 still listed it among model genes |
| 64,305 peaks detected in ≥30 of 2,718 nuclei; 30,000 kept | `results/run_manifest.json` (`n_peaks_detected_min_cells`); `wf:select_atac_features`; printed in notebook section 6 | Feature counts | — |
| First LSI component r = 0.99 with log10 fragments; dropped | `results/lsi_depth_correlations.tsv` (component 1: 0.992) | Depth association of LSI 1 | — |
| Resolution 0.8 is a fixed choice, set before looking at the results and not tuned | `run:parse_args` default `--resolution 0.8`; no resolution sweep was run for the main clusterings (only `wf:resolution_check` for the two small groups) | How the setting was chosen | That 0.8 is optimal; cluster counts and ARI depend on it |
| 14 RNA-only and 10 ATAC-only clusters (resolution 0.8, 20 neighbors, cosine) | `results/cell_metadata.tsv.gz` (`leiden_rna`, `leiden_atac`); `wf:build_embedding_graph` | Cluster counts at these settings | Cluster counts at other settings |

## 5. MultiVI

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| MultiVI treats zero counts in a modality as unmeasured | scvi-tools 1.4.2 `src/scvi/module/_multivae.py` lines 602–603 (`mask_expr = x_rna.sum(dim=1) > 0`) | Model behaviour | — |
| No nucleus with zero counts; minimum 60 RNA and 479 ATAC counts | `results/observed_library_check.tsv` | Every nucleus is treated as paired | — |
| Main `batch_key` should be the modality when mixing paired and single-modality data | Official scvi-tools MultiVI tutorial notebook (scverse/scvi-tutorials, `multimodal/MultiVI_tutorial.ipynb`, "important" note) | Guidance quoted | — |
| 128 hidden units, 11 latent dimensions; early stopping on `reconstruction_loss_validation`, patience 50; stopped at epoch 377; best validation reconstruction loss at 327 | `results/run_manifest.json`; `results/multivi_training_history.tsv` (argmin of `reconstruction_loss_validation`; total `validation_loss` and `elbo_validation` are lowest at 285); scvi-tools 1.4.2 `src/scvi/model/_multivi.py` lines 342–343 | Model size and convergence for this run | That the representation is meaningful |

## 6. Latent audit

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| Overall: dimension 3 with TSS fraction r = -0.70; dimensions 3 and 6 with log10 UMIs \|r\| ≈ 0.5 | `results/multivi_latent_qc_correlations.tsv` | Linear association across all nuclei | That it is technical (cell types differ in depth and TSS fraction) |
| Within populations: dimension 3 with TSS -0.42, with UMIs 0.10; dimension 9 with log10 fragments -0.54, TSS 0.48, FRiP 0.45 | `results/multivi_latent_qc_correlations_within_population.tsv` (`wf:latent_qc_correlations_within`, both variables centered on population means) | Which dimension tracks ATAC depth among nuclei of one type | Nonlinear structure |
| 17 joint clusters | `results/cell_metadata.tsv.gz` (`leiden_multivi`) | Count at resolution 0.8 | Stability across seeds (see seed check) |

## 7. Labels and accessibility evidence

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| Rule: score ≥ 0.3 and margin ≥ 0.2; all 17 clusters named; 1,074 / 902 / 342 / 269 / 102 / 29 nuclei | `wf:assign_populations`; `results/cluster_marker_scores.tsv`; population counts in `results/cell_metadata.tsv.gz` | Labels under this rule | Cell identity; subtypes inside a label |
| MS4A1 and CD79A promoters 44% and 42% in B cells, ≤ 5% elsewhere | `results/marker_evidence.tsv` (`promoter_fraction_accessible`) | Observed promoter detection | Regulation of those genes by those peaks |
| CLEC4C promoter 31% in pDCs, ≤ 6% elsewhere | Same | Same | Same |
| CD3D/E/G 12–14% in T cells, ≤ 9% elsewhere; NKG7 35% in NK, 24% in pDCs | Same | Same | Same |
| CD14-monocyte marker promoters ≤ 14% everywhere; LYZ, VCAN 7–14% in monocytes, ≤ 1% in other populations; CD14, FCN1 barely differ | Same | Same | Same |
| Most specific loci shown: MS4A1 0.40, IL7R 0.22, CLEC4C 0.19; monocyte markers ≤ 0.05 | `results/marker_locus_specificity.tsv`; rule `wf:marker_locus_specificity`, `wf:choose_marker_loci` | Single-peak specificity within ±20 kb of the annotated TSS span | Specificity of distal elements beyond 20 kb |
| For the monocyte loci shown (VCAN, MS4A7), at most 16% of the population's nuclei have a fragment in the most specific peak | `results/marker_locus_specificity.tsv` (`fraction_in_population` 0.114 and 0.157) | Specificity is capped by sparse detection | — |
| "Suggests that their chromatin differences lie mostly outside these marker loci" | Inference from ATAC-only clustering separating monocytes (section 8) while no marker locus is specific | Plausible interpretation | Where the differences lie; no region-level analysis was run |
| 10x `interval` is the span of annotated TSSs | Session check: FAM138A `chr1:36080-36081`, ACTB `chr7:5530600-5530602`; median width 1 bp | How the locus plot is anchored | — |

## 8. Agreement between representations

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| ARI 0.65 (RNA vs ATAC), 0.58 (joint vs RNA), 0.45 (joint vs ATAC); 14 / 10 / 17 clusters | `results/clustering_agreement.tsv` | Pairwise agreement at these settings | Agreement is confounded by cluster number, as stated |
| B cells 142/126 in RNA, 267/269 together in ATAC; FCGR3A 101/102 apart in RNA, 26 with CD14 in ATAC; NK 206/101 and 221/104; pDC 29/29 and 28/29 | Cross-tabulation of `population` with `leiden_rna` / `leiden_atac` in `results/cell_metadata.tsv.gz`; figure 13 | How single-measurement clusters partition the joint-defined populations | Independent accuracy (populations come from joint clusters) |
| Three T-cell groups in each measurement, with the marker fractions in the table | `results/t_cell_state_markers_by_rna_cluster.tsv`, `results/t_cell_state_markers_by_atac_cluster.tsv` | Marker detection by cluster among nuclei labelled T cells | Validated T-cell states; CD4 detection is low in single nuclei |
| RNA-only cluster 13: 27 nuclei, FOXP3 0.41, IL2RA 0.59, RTKN2 0.74 | Same RNA table | Regulatory-T-like marker profile | Treg identity (no protein, no suppression assay) |

## 9. Joint-specific divisions and merges

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| Each T-cell group split in two with a depth difference: naive CD8 11,845 vs 7,044 (1.7-fold), naive CD4 11,295 vs 8,029 (1.4-fold), memory CD4 10,510 vs 13,936 (1.3-fold) | Joint clusters 5/11, 3/13, 4/7: medians of `atac_fragments` by `leiden_multivi` in `results/cell_metadata.tsv.gz` | Depth difference between the paired clusters | That depth causes the split |
| Naive CD8 pair: similar detection of every panel marker; naive CD4 cluster 13: CCR7 0.48 vs 0.70, IKZF2 0.16 vs 0.03, IL2RA 0.09 vs 0.01, RTKN2 0.09 vs 0.01; holds 5 of 27 Treg-like nuclei | `results/t_cell_state_markers.tsv`; cross-tab RNA 13 × joint (13: 5) | Marker differences within the 16-gene panel | A gene-wide comparison (not run for T-cell pairs) |
| Single-measurement clusterings place most nuclei of each pair in the same cluster | Cross-tab: joint 5→RNA 2 (197/207), 11→RNA 2 (104/109), 5→ATAC 1 (198), 11→ATAC 1 (96); joint 3→RNA 1 (210/219), 13→RNA 1 (92/104), 3→ATAC 3 (206), 13→ATAC 3 (87) | Agreement | — |
| Memory CD4 cluster 7 vs 4: CCR7 0.20 vs 0.39, GZMK 0.10 vs 0.01, CCL5 0.16 vs 0.03; neither single measurement separates them; the two CD4 divisions are the candidates for new states | `results/t_cell_state_markers.tsv`; cross-tab joint 4→RNA 4 (159), 7→RNA 4 (142), 4→ATAC 2 (182), 7→ATAC 2 (151) | Marker difference between joint clusters | Effector-memory identity ("consistent with" only) |
| Cluster 15: 60 nuclei; CD4 detected in 17% and CD8A in 37% of nuclei; highest T-cluster median fragments (19,640) and doublet score (0.081) | `results/t_cell_state_markers.tsv`; medians by joint cluster | Doublet-like profile at cluster level | Co-detection of CD4 and CD8A in the same nuclei (not tested); that these are doublets |
| CD14 monocyte joint clusters table | `results/cd14_monocyte_cluster_qc.tsv` | Medians as printed | — |
| "About three-fold" | 13,739/4,224 = 3.25; 14,281/4,183 = 3.41 | Ratio range | — |
| RNA-only monocyte clusters mix depth groups in similar proportions | Cross-tab depth group × `leiden_rna`: RNA 0 = 164 higher / 348 lower; RNA 3 = 110 / 235 (both 32% higher) | Independence of the RNA split from depth | — |
| ATAC-only places 842 of 864 compared monocytes in one cluster | Cross-tab depth group × `leiden_atac`, DC-like nuclei removed: 586 lower + 256 higher in ATAC 0 | Same | — |
| Only MTRNR2L12 and MTRNR2L8 at adjusted p < 0.05, both higher in lower-depth group | `results/cd14_monocyte_depth_group_rna_comparison.tsv` (Wilcoxon, Benjamini–Hochberg, MT-* excluded, 38 DC-like nuclei removed; 588 vs 276 nuclei) | No other gene passes this test | Absence of any expression difference (power is limited; one donor) |
| MTRNR2L share correlates with mitochondrial fraction: ρ = 0.63 and 0.37 | `results/depth_group_gene_mito_link.tsv` (`wf:mitochondrial_link`) | Association across 2,718 nuclei | Mapping mechanism |
| DC-like group: RNA-only cluster 10, 38 nuclei, FCER1A 0.76, CD1C 0.76, CLEC10A 0.74; ATAC-only cluster 8 holds 32, 86% of that cluster; joint cluster 6 holds all 38, 23% of it | `results/dc_markers_by_rna_cluster.tsv`; `results/small_group_resolution_check.tsv` (ATAC-only row and resolution 0.8 row) | Placement in each clustering | DC identity beyond marker detection |
| Treg-like nuclei spread over joint CD4 clusters; ATAC-only does not separate them | Cross-tab RNA 13 × joint (4: 18, 7: 4, 13: 5); resolution check ATAC-only share 0.07 | Placement | — |
| Resolution 0.8→3.0 (17→41 clusters): from 1.6 upward all 38 DC-like nuclei share one joint cluster and make up 52–54% of it; Treg-like nuclei at most 34% | `results/small_group_resolution_check.tsv` | DC-like nuclei form only a slim majority of one cluster at higher resolution; Treg-like nuclei never a majority | Behaviour with a larger latent space or other graph settings |

## Seed check

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| Seeds 23, 42, 101 retrained on the same input; 215–353 epochs | `results/seed_check.tsv` (`seed:main`); `results/seed_check/seed_*.tsv` | Training variation only (graph and Leiden seed fixed at 17) | Variation from feature selection, QC, or data |
| Joint clusters change: 17–20 clusters, ARI 0.44–0.66 vs the main run | `results/seed_check.tsv` (`n_joint_clusters`, `ari_vs_main_joint`) | Instability of individual joint clusters | — |
| CD14 monocytes in 4–5 joint clusters of ≥50 nuclei, depth span 3.4–4.2-fold, in every seed | `results/seed_check.tsv` (`cd14_clusters_ge50`, `cd14_depth_ratio_max_min`) | The depth division recurs | That the same nuclei split the same way |
| DC-like at most 27%, Treg-like at most 19% of any joint cluster, in every seed | `results/seed_check.tsv` (`dc_share_of_best_cluster` 0.22–0.27, `treg_share_of_best_cluster` 0.085–0.19) at resolution 0.8 | The merges recur at resolution 0.8 | Behaviour at higher resolution in other seeds |
| Joint vs RNA-only ARI 0.36–0.58, vs ATAC-only 0.31–0.45, all below RNA vs ATAC 0.65 | `results/seed_check.tsv`; `results/clustering_agreement.tsv` | Ordering holds in every seed | — |
| The T-cell divisions were not checked across seeds | `seed:recurrence_metrics` has no T-cell metric | Stated on the page as a reason not to report individual joint clusters from one run | — |

## Going further

| Claim | Source | Establishes | Does not establish |
|---|---|---|---|
| The scATAC tutorial calls cluster-associated peaks with PEAKVI `differential_accessibility` | `source/epigenomics/scatac.rst` section 7 | Cross-reference is accurate | — |
| MultiVI provides `differential_accessibility` and `differential_expression` | scvi-tools 1.4.2 `src/scvi/model/_multivi.py` lines 724 and 853 | Methods exist in the pinned version | — |
| A comparison between the CD14 depth groups "would be confounded with depth" | The groups differ about three-fold in ATAC fragments (`results/cd14_monocyte_cluster_qc.tsv`) | Logical consequence of the design | No DA was run |
| MultiVI paper, "Integration of independent studies": three PBMC studies; without study/technology terms, cells separated by sample (accessibility) and technology (expression); correcting for both removed this while preserving cell types | Ashuach et al. 2023, Nature Methods, PMC10406609 full text (Europe PMC XML), section "Integration of independent studies" and Supplementary Fig. 4 reference | What the paper reports | Behaviour on other datasets |

## Interpretations stated on the page

- "Suggests nuclei whose own RNA was lost or diluted during isolation": from the
  ATAC-only medians; not tested.
- "Treat this division provisionally as most likely technical": rests on the
  depth difference, the absence of an RNA difference other than MTRNR2L genes,
  and the MTRNR2L correlation. The groups also differ in FRiP and TSS fraction.
  A refit with depth as a covariate was proposed, not run.
- "Consistent with a more effector-like memory state": marker detection only.
- "In the joint representation, these small groups sit close to their larger
  neighbors": from the resolution check and the seed check.
- "MultiVI is designed for a different use": describes the method's purpose
  (Ashuach et al. 2023; scvi-tools tutorial), not a result of this analysis.

## Introduction (reordered on the user's request)

The conclusion paragraph now comes first and the evidence is condensed to four
bullets. Each bullet restates a claim already listed in sections 8, 9, and the
seed check; "no detectable expression difference other than two
mitochondrial-derived transcripts" keeps the qualifier from section 9.

## Adversarial review

A separate reviewer read the page against the result files and the code,
looking only for wording stronger than the tests. It recomputed the numbers.
Corrections made from that review:

- FRiP comparison corrected to 0.67, with the comparison group named ("pass both").
- "Best validation loss" corrected to "best validation reconstruction loss".
- T-cell divisions reclassified: all three differ in depth; the two CD4 divisions also differ in markers.
- Section 6 rewritten around within-population correlations, because overall correlations mix in cell type.
- Resolution check reworded: DC-like nuclei form a slim majority (52–54%) of one cluster from resolution 1.6.
- Monocyte locus specificity explained by sparse detection, not shared markers.
- Region-level inference softened to "suggests ... does not test where".
- Intro states that populations were named from joint clusters.
- Monocyte division called "provisionally most likely technical".
- Mitochondrial-derived gene rationale stated as association plus "likely".
- Minor fixes: "most of the 3,009", upper-depth rule attributed to the scATAC tutorial, "largely separates (32 of 38)", co-detection wording, import-order wording, ambient-RNA wording.
- The 64,305 peak count and the early-stopping settings are now recorded in result files and cited.

The earlier first-run results (model with mitochondrial genes, then with MTRNR2L12) are not cited on the page. They are why the gene exclusion rule exists. A regulatory T-like joint cluster seen in the second run did not appear in the final run, which motivated the seed check.
