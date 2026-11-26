# vitiligo-panel-audit

Code and results for the paper:

**Do Machine-Learning Gene Signatures for Vitiligo Generalise? An Audit Against Random Baselines**
Aryan Misra and Yasha Hasija, Department of Biotechnology, Delhi Technological University (2026).

The study takes 14 published vitiligo gene panels and tests them, without retraining, on an independent bulk RNA-seq cohort (GEO: GSE311770; 105 skin biopsies from 28 patients and 27 healthy controls). Each panel is scored with a direction-aware within-sample rank score and compared with 1,000 expression-matched random gene sets, with Benjamini-Hochberg FDR control across panels.

## Main result

| Comparison | Replicates (FDR < 0.05 vs matched random genes) | Reversed (AUC CI below 0.5) |
|---|---|---|
| Lesional vs healthy (24 vs 27) | P05, P06, S04 (P02 borderline) | P04, P09, S02 |
| Perilesional vs healthy (27 vs 27) | P05, P06, S03 | P04, P09 |
| Non-lesional vs healthy (27 vs 27) | none | P09, S02 |

A curated pigment gene set with no data-driven selection reaches AUC 0.94 in lesional skin. Full tables are in `results/`.

Post hoc checks: in a site-matched comparison of lesional with perilesional skin from the same 23 patients, P02 and P06 stay strongly separated and P04, P09 and S02 stay reversed, while P05 and S04 weaken. Of the 12 inferred gene directions that could be checked in GSE75819, 11 agree. In a second, single-cell cohort (GSE298871, pseudobulk, 6 lesional vs 6 healthy), P05 and P06 replicate again, P04 and P09 again point in the opposite direction (AUC 0.17 and 0.11; exact p = 0.065 and 0.026), and S04 does not replicate.

Machine-learning comparison (post hoc): logistic regression trained on each panel's genes (5-fold CV, 5 repeats) reaches AUC 0.37 to 0.97 in lesional skin, but models trained on the panels' matched random gene sets reach median AUCs of 0.64 to 0.83, and after FDR correction no panel beats its random sets (lowest q = 0.098; SVM and random forest agree). P04 and P09, whose published directions are reversed, reach AUC 0.94 and 0.97 once retrained, because training flips the direction. A genome-wide LASSO (nested CV) reaches AUC 0.985; of the 9 genes chosen in at least half of the models, 7 are melanocyte-enriched in the GSE298871 single-cell data (vs 1.75% of all genes; binomial p = 1.8e-11). Trained on GSE311770 and tested once on GSE298871, only P06 holds up (AUC 0.97), and the training-free score is equal to or better than the retrained model for 9 of 14 panels. A small neural network (MLP) gives the same answer (lowest q = 0.12). In a third cohort (OMIX691, 10 vitiligo vs 5 healthy), no panel passes FDR (lowest q = 0.20), P04, P08 and P09 point the opposite way, the genome-wide LASSO trained on GSE311770 drops to AUC 0.44, and the training-free score is equal to or better than retraining for 10 of 14 panels.

## Pre-specified analysis

All panels, gene directions, null gene sets, reference sets, covariates and the random seed were written to `gene_sets_locked.json` before any panel was scored. Its SHA-256 fingerprint is in `gene_sets_locked.sha256`, and `code/stage45_evaluate.py` stops if the file has changed:

```
09addbb44fb273693a507cd042f9601cf50215134421bbe854b26ad30e4bbf04  gene_sets_locked.json
```

Decisions taken after locking (replication criterion, sebaceous threshold, pigment-removal rule, exploratory and post hoc analyses) are listed in `docs/Verification_log.md`.

## Repository layout

```
gene_sets_locked.json     locked gene sets, null sets and analysis rules
gene_sets_locked.sha256   fingerprint of the locked file
code/                     analysis scripts (run from the repository root)
data/                     sample metadata; raw inputs are downloaded (see data/README.md)
results/                  all result tables and figures, as used in the paper
docs/                     QC report, panel extraction notes, results summary, verification log
verification/             GEO records and GSE75819 probe values backing factual claims (see verification/README.md)
```

## Reproducing the analysis

Python 3.11 or later. Install dependencies:

```
pip install -r requirements.txt
```

Download the public inputs into `data/` as described in `data/README.md`, then run from the repository root:

| Step | Script | Output |
|---|---|---|
| 1. QC, filtering, log-CPM, PCA | `code/qc.py` | `data/meta.tsv`, `data/logcpm.tsv.gz`, `data/pca.tsv` |
| 2. Panel table | `code/build_panels.py` | `panels_summary.tsv`, `panels_gene_level.tsv` |
| 3. Infer missing directions from GSE65127 | `code/infer_directions.py` | `direction_check_GSE65127.tsv` |
| 4. Lock gene sets and nulls | `code/stage3_gene_sets.py` | `gene_sets_locked.json`, `.sha256` |
| 5. Score and evaluate | `code/stage45_evaluate.py` | `results_all.tsv`, `scores_per_sample.tsv` |
| 6. Melanocyte dependence (exploratory) | `code/stage6_exploratory.py` | `results_exploratory_melanocyte.tsv` |
| 7. Independent recomputation and audit | `code/audit.py` | `audit_*.tsv` |
| 8. Post hoc null and VP-008 checks | `code/sens_null_vp008.py`, `code/sens_null_10k.py` | `sens_*.tsv`, `results_exploratory_melanocyte_v2.tsv` |
| 9. Figures | `code/make_figure.py`, `code/make_figure_paper.py` | `fig_forest_AUC*.png/pdf` |
| 10. Post hoc: site-matched lesional vs perilesional | `code/posthoc_site_matched.py` | `posthoc_site_matched.tsv` |
| 11. Post hoc: check inferred directions in GSE75819 | `code/posthoc_direction_GSE75819.py`, then `code/posthoc_direction_rescore.py` | `posthoc_direction_*.tsv` |
| 12. Post hoc: replication-rule sensitivity (needs `results_all.tsv` and `sens_null_10k.tsv` in the working directory) | `code/posthoc_rule_sensitivity.py` | `posthoc_rule_sensitivity.tsv` |
| 13. Post hoc: second cohort GSE298871 (pseudobulk; needs `GSE298871_RAW.tar` from GEO extracted into a folder) | `code/sc_load_GSE298871.py`, `code/posthoc_GSE298871_pseudobulk.py`, `code/posthoc_GSE298871_celltypes.py` | `posthoc_GSE298871_*.tsv` |
| 14. Post hoc: machine-learning comparison (E1: panels vs random sets, LR; E1b: SVM and random forest; E2: genome-wide LASSO and melanocyte check; E3: GSE311770 to GSE298871). E2, E2b and E3 take the `sc.pkl` from step 13; E2b needs E2's output and E3 needs `posthoc_GSE298871_pseudobulk.tsv` in the working directory | `code/ml_e1_panels_vs_random.py`, `code/ml_e1b_other_models.py`, `code/ml_e2_lasso_genomewide.py`, `code/ml_e2b_melanocyte_specificity.py`, `code/ml_e3_cross_cohort.py` | `ml_e*.tsv`, `ml_e*_summary.json` |
| 14b. Post hoc: E1 with a small neural network (MLP, one hidden layer of 8 units, 200 random sets) | `code/ml_e1c_mlp.py` | `ml_e1c_mlp.tsv` |
| 14c. Post hoc third cohort: OMIX691 (NGDC; 10x scRNA-seq, 10 vitiligo, 5 healthy), pseudobulk training-free scoring and models trained on GSE311770 tested once. Download the 15 `OMIX691-20-*.tar.gz` files from https://download.cncb.ac.cn/OMIX/OMIX691/ and extract them into one folder | `code/posthoc_OMIX691.py <folder>` | `posthoc_OMIX691_pseudobulk.tsv`, `ml_e3_OMIX691.tsv`, `ml_e3_OMIX691_summary.json` |
| 15. Machine-learning figure and summary table | `code/make_figure_ml.py` | `fig_ml.pdf/png`, `ml_table_summary.tsv` |

Scripts write their outputs to the working directory; the copies used in the paper are in `results/`. To reproduce the published numbers exactly, keep the committed `gene_sets_locked.json` and run steps 1 and 5 onwards. Re-running step 4 rebuilds the lock from the current NCBI gene annotation, which may differ slightly from the version used here.

## Score

For sample *j*, genes are ranked within the sample and scaled to (0, 1]. With *U* and *D* the panel's up- and down-regulated genes:

s_j = mean rank of U genes − mean rank of D genes

A term is omitted when its set is empty. The score needs no training and no normalisation across samples.

## Data sources

- GSE311770: Berrevoet D. et al., Front. Immunol. 2026, doi:10.3389/fimmu.2026.1829996
- GSE298871: Brunner P. M. et al., J. Allergy Clin. Immunol. 2025 (second cohort, post hoc)
- GSE65127: Regazzetti C. et al., J. Invest. Dermatol. 2015, doi:10.1038/jid.2015.335
- MSigDB Hallmark Interferon-gamma response, v2026.2
- NCBI Gene (Homo sapiens gene_info) for symbol to Ensembl mapping

The 14 source papers and their DOIs are listed in `results/panels_summary.tsv` and `results/panels_gene_level.tsv`.

## License

Code: MIT (see `LICENSE`). Result tables and figures: CC BY 4.0.
