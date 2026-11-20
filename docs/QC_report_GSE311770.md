# GSE311770 data audit and go/no-go gate (20 Nov 2026)

## Data
- File: GSE311770_all_samples_rsem_gene_counts.tsv.gz (RSEM gene counts, Ensembl IDs), 63,140 genes x 105 samples. Raw reads not deposited (privacy).
- Platform: Element AVITI. Library size: median 4.2 M reads (range 1.6-6.7 M); no sample below 25% of median.
- 18,848 genes pass CPM > 1 in >= 20 samples.

## Sample structure (verified from library names)
| region | batch A | batch B | total |
|---|---|---|---|
| control | 9 | 18 | 27 |
| lesional | 13 | 11 | 24 |
| non-lesional | 15 | 12 | 27 |
| perilesional | 15 | 12 | 27 |

- 27 patient IDs on GEO vs 28 in the paper. VP-008 has two sample sets (an "-X" set with lesional/perilesional/non-lesional, and a perilesional + non-lesional pair without X). Both are female; correlation between them equals typical between-patient correlation, so identity cannot be resolved from data. Decision: treat as two individuals (consistent with the paper's 28). This does not affect the primary analysis: only the X set has a lesional sample.
- VP-001, VP-013 and VP-019 have no lesional sample (27 IDs - 3 = 24 lesional).
- Label check passed: melanocyte genes (TYR, DCT, PMEL, MLANA) are near zero in lesional and present in all other regions.

## Confounders found
1. Batch: both batches contain patients and controls (not fully confounded), but controls are 2:1 in batch B. Batch differs within controls on PC1/PC2 (p ~ 0.02-0.03). Action: batch as covariate; report per-batch AUC as robustness.
2. Sebaceous-gland content: PC1 (21% of variance) is driven by sebaceous genes (ELOVL3, AWAT2, FADS2, THRSP, MGST1, CIDEA, SOAT1, KRT79). 9 of 27 controls (all batch B) vs 1 of 24 lesional samples sit in the high-sebaceous group, consistent with controls being sampled from the upper arm and patients mostly from arm or leg. Action: compute a sebaceous score as a covariate; sensitivity analysis excluding the high-sebaceous group; flag any published panel containing lipid/sebaceous genes.
3. Sex: controls 20 F / 7 M; patients 13 F / 14 M. Action: exclude sex-chromosome genes from random gene sets; sex as covariate.
4. Technical noise: negative PC1 loadings dominated by RN7SL/snoRNA genes. Action: restrict analysis to protein-coding genes.

## Biological note
CXCL9/CXCL10 are low in lesional skin and highest in perilesional skin. Lesional-vs-control discrimination is therefore expected to be driven mainly by melanocyte loss, not immune activity; the perilesional comparison carries the immune signal.

## Novelty gate
The 2026 rank-based meta-analysis (Ge, medRxiv, posted 9 Feb 2026) predates the public release of GSE311770 (3 Jun 2026), so it cannot have used this cohort, and it analysed differential expression, not published biomarker panels.

## Verdict: GO, with the four confounder controls above built into the analysis.
