# vitiligo-panel-audit

Code and results for the paper:

**Do Published Transcriptomic Biomarkers of Vitiligo Replicate? An Independent RNA-seq Validation and Audit of Fourteen Gene Panels**
Aryan Misra, Delhi Technological University (2026).

The study takes 14 published vitiligo gene panels and tests them, without retraining, on an independent bulk RNA-seq cohort (GEO: GSE311770; 105 skin biopsies from 28 patients and 27 healthy controls). Each panel is scored with a direction-aware within-sample rank score and compared with 1,000 expression-matched random gene sets, with Benjamini-Hochberg FDR control across panels.

## Main result

| Comparison | Replicates (FDR < 0.05 vs matched random genes) | Reversed (AUC CI below 0.5) |
|---|---|---|
| Lesional vs healthy (24 vs 27) | P05, P06, S04 (P02 borderline) | P04, P09, S02 |
| Perilesional vs healthy (27 vs 27) | P05, P06, S03 | P04, P09 |
| Non-lesional vs healthy (27 vs 27) | none | P09, S02 |

A curated pigment gene set with no data-driven selection reaches AUC 0.94 in lesional skin. Full tables are in `results/`.

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

Scripts write their outputs to the working directory; the copies used in the paper are in `results/`. To reproduce the published numbers exactly, keep the committed `gene_sets_locked.json` and run steps 1 and 5 onwards. Re-running step 4 rebuilds the lock from the current NCBI gene annotation, which may differ slightly from the version used here.

## Score

For sample *j*, genes are ranked within the sample and scaled to (0, 1]. With *U* and *D* the panel's up- and down-regulated genes:

s_j = mean rank of U genes − mean rank of D genes

A term is omitted when its set is empty. The score needs no training and no normalisation across samples.

## Data sources

- GSE311770: Berrevoet D. et al., Front. Immunol. 2026, doi:10.3389/fimmu.2026.1829996
- GSE65127: Regazzetti C. et al., J. Invest. Dermatol. 2015, doi:10.1038/jid.2015.335
- MSigDB Hallmark Interferon-gamma response, v2026.2
- NCBI Gene (Homo sapiens gene_info) for symbol to Ensembl mapping

The 14 source papers and their DOIs are listed in `results/panels_summary.tsv` and `results/panels_gene_level.tsv`.

## License

Code: MIT (see `LICENSE`). Result tables and figures: CC BY 4.0.
