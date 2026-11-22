# Verification log (20 Nov 2026, second pass)

Each earlier claim re-checked against its primary source. Status: CONFIRMED, CORRECTED, or OPEN.

## Dataset (GSE311770)
| Claim | Check | Status |
|---|---|---|
| 105 samples: 24 L, 27 PL, 27 NL, 27 controls | Sample titles + library names; matches paper Methods (24/27/27 + 27 controls) | CONFIRMED |
| 28 patients in paper vs 27 IDs on GEO | VP-008 has two sample sets, same sex, inseparable by expression; treated as two people | CONFIRMED (identity OPEN, no effect on primary analysis) |
| Two batches | Paper states two library-prep batches; library names show two plates with independent numbering (1-52 cols 12->6; 1-53 cols 1->7); all 52 shared indices appear once per plate | CONFIRMED (batch = plate is an inference, not stated by authors) |
| Labels correct | Melanocyte genes ~0 in lesional, present elsewhere | CONFIRMED |
| Sebaceous confounder | PC1 top loadings ELOVL3, AWAT2, FADS2, THRSP, MGST1, CIDEA, SOAT1, KRT79; 9/27 controls vs 1/24 lesional in high-PC1 group | CONFIRMED |
| Controls not site-matched | Paper: controls upper arm; patients preferentially arm or leg | CONFIRMED |
| Novelty vs 2026 meta-analysis | Preprint posted 9 Feb 2026; GSE311770 public 3 Jun 2026 | CONFIRMED |
| No included panel used GSE311770 | All discovery datasets recorded (GSE65127, GSE53146, GSE75819, GSE30999, GSE98793, blood sets) | CONFIRMED |

## Panels
| Item | Status |
|---|---|
| Gene count: earlier stated "68 genes, 67 detected" | CORRECTED: 71 gene entries, 64 unique genes, 63 detected (IFNG not detected) |
| P07 autophagy directions | CONFIRMED: 44-gene list, first 20 = up (top-5 up and top-5 down named in text match the list ends) |
| P04 (Zhao 2023) DCT/KIF1A direction | CORRECTED: paper reports both "significantly higher in vitiligo" in GSE65127, which contradicts melanocyte loss. Now scored as reported (up). |
| P09 (AbdElneam 2024) direction | CORRECTED: genes were selected as "upregulated" from GEO2R healthy-vs-lesional; scored as reported (up). Likely a reversed GEO2R contrast. |
| P06 (Zhang 2021) TYR/TYRP1/DCT direction | CORRECTED: not reported for the authors' own data (only literature cited) -> inferred from discovery data |
| P09 discovery dataset | CORRECTED: GSE65127 only (GEO2R), not "integrated" |
| S02 discovery dataset | CORRECTED: GSE65127 + GSE53146 (vitiligo), GSE98793 (MDD) |
| S05 HMGA1 | CORRECTED: down (reported in blood by qPCR) |
| GSE53148 duplicates GSE53146 | CONFIRMED: GEO lists GSE53148 as a SuperSeries of GSE53146 (human) and GSE53147 (mouse) |
| GSE90880 is blood | CONFIRMED (series title) |
| GSE75819 has no healthy controls | CONFIRMED (15 patients x lesional/non-lesional; stated in P07's own Methods) |

## Direction rule (fixed now, before any scoring)
1. Primary analysis: direction as reported by the authors.
2. Not reported (13 gene entries): inferred from the paper's own discovery dataset with healthy controls (GSE65127 covers P05, P06, S02), never from GSE311770.
3. Sensitivity analysis: direction-agnostic AUC (max of AUC, 1-AUC), labelled optimistic. A large gap between 1 and 3 identifies direction errors in the source papers.

## Still OPEN
- INDIACom-2027 track scope (fetched page showed the 2026 edition); confirm the 2027 call for papers.
- Xia 2026 full text (excluded for now).
- Sample-level sebaceous/sex covariates to be added at scoring stage.

## Direction check against GSE65127 (lesional vs healthy, 10 vs 10)
- P04 (Zhao 2023) states DCT and KIF1A are higher in vitiligo in GSE65127. GSE65127 shows the opposite: DCT log2 diff -4.85, KIF1A -0.51 (both p = 1.8e-4). The reported direction is wrong for the paper's own validation dataset.
- P09 (AbdElneam 2024) selected 10 melanocyte genes as "upregulated" from GEO2R on GSE65127. All 10 are strongly lower in lesional skin (log2 diff -0.9 to -4.7, p <= 5.8e-4). Consistent with a reversed GEO2R contrast.
- Other mismatches (P01, P07, P08, S03, S05) are small and non-significant, and those papers used other discovery datasets; not treated as errors.
- 13 unreported directions filled from GSE65127 (direction_check_GSE65127.tsv).

## Stage 3 lock
- gene_sets_locked.json, SHA-256 09addbb44fb273693a507cd042f9601cf50215134421bbe854b26ad30e4bbf04, created before any panel was scored on GSE311770.
- Universe: 13,813 detected protein-coding autosomal genes. IFN-gamma: MSigDB Hallmark v2026.2 (191/200 detected). Melanogenesis: GO:0042438 (22 genes). Melanocyte identity: GO:0030318 + HPA melanocyte-enriched (20 genes). 1,000 expression-matched null sets per panel (seed 20261120).
- Note: the curated pigment lists do not contain GPR143, TRPM1, KIF1A or SNCA, which the provisional flag in stage 1 included. The curated definition is the one used.

## Audit of stages 4-6 (third pass)
| Item | Status |
|---|---|
| AUCs reproduced independently from raw counts with separate code | CONFIRMED (max difference 0.000000) |
| Region labels re-parsed from library names | CONFIRMED |
| Sex call | CONFIRMED (gap 0.51 vs 6.09 log2CPM) |
| "5 of 14 panels replicate" | CORRECTED: 3 replicate after FDR on the matched-null test (P05, P06, S04); P02 and P07 only before FDR |
| Logistic covariate adjustment | REPLACED: perfect separation made odds ratios meaningless; linear model (score ~ group + batch + sex + sebaceous) used instead; all replicating/reversed panels keep sign and significance |
| Scoring-method robustness | ADDED: z-score scoring leaves all replicate/reversed verdicts unchanged; P08 and P01 are method-sensitive |
| "Lesional skin is immunologically quieter than healthy" | CORRECTED: not supported (CI includes 0.5); only lesional < paired non-lesional is supported (p = 0.005) |
| "IFN-gamma hallmark works in perilesional skin" | CORRECTED: not supported (CI 0.47-0.79) |
| S02 reversal | QUALIFIED: 2 of its 3 gene directions were inferred from GSE65127 |
| Lu 2025 / GSE53148 statement | CORRECTED: GSE53148 mixes human and mouse SubSeries; not a duplicate of Lu's own training data |
| Replication criterion | DISCLOSED: locked file fixed BH correction but not which test; criterion (CI > 0.5 and q_null < 0.05) stated before writing; both versions reported |
| Sebaceous exclusion threshold (75th percentile) and minimum-2-genes rule | DISCLOSED: set at analysis time, not in the locked file |
| Melanocyte-dependence analysis | LABELLED EXPLORATORY |

## Second-cohort check (20 Nov 2026, evening)
| Item | Status |
|---|---|
| GSE298871 described earlier as "incomplete bulk RNA-seq (17 of 44 samples)" | CORRECTED: all 17 GEO samples are 10x Genomics single-cell libraries (Chromium 5', enzymatic whole-skin dissociation, Cell Ranger 7.2). The paper's bulk RNA-seq (15 patients, 14 controls) was not deposited. |
| Use GSE298871 as a second validation cohort? | NOT DONE (decision): pseudobulk from enzymatically dissociated cells does not preserve tissue cell composition, and the replicating panels are driven by melanocyte content; with 6 vs 6 samples, agreement or disagreement would be uninterpretable. |
| Use GSE65127 as a second cohort? | NOT DONE (decision): it was used to build or validate 8 of the 14 panels (and to infer directions); independent only for P01, P02, P07; adds little and risks circularity. |
| Consequence for the paper | GSE311770 is the only independent bulk cohort available; stated as a limitation and as evidence for the gap (all older vitiligo cohorts were reused to build the panels). |

## Pass 4 (20 Nov 2026): response to external critique (post-hoc, not in the locked plan)
- Null baselines: locked bins used mean/SD over all 105 samples; null medians drift (P02 0.65, P06 0.62, S04 0.39, P09 0.33). Re-ran with SD from controls only, mean+SD from controls, expression-only, and unrestricted random genes (sens_null_vp008.py -> sens_null_variants.tsv). L vs HC (re-run with the exact locked 13,813-gene universe): P05, P06, S04 pass FDR under all variants (q <= 0.042); P02 passes under all four alternatives (q = 0.007) and under a re-draw of the locked binning (q 0.045 vs locked 0.073). With 10,000 draws of the locked binning (sens_null_10k.py) P02 q = 0.052: borderline; headline stays 3 of 14. P07 fails under all. P04/P09/S02 reversed under all.
- VP-008: dropping the non-X set changes PL/NL AUCs by at most 0.025 (sens_vp008.tsv).
- Correction: lesional count is 27 IDs minus VP-001, VP-013, VP-019 = 24 (earlier text omitted VP-013).
- Melanocyte residual: P02 (DCT only) keeps AUC 0.66-0.81 after linear, cubic or rank adjustment, so non-linearity does not explain it; the melanocyte proxy score also contains the panel genes (P02, P06). Treat residual AUCs as imperfect adjustment, not melanocyte-independent signal.
- S04: robust to covariates (adj p 0.005) but batch-heterogeneous (batch A 0.95, batch B 0.73); explanation unknown.

- Expression filter: 61 of 64 unique panel genes pass CPM>1 in >=20 samples (SLC45A2, IFNG, IL1B fail); P09 scored on 9/10 genes, S03 on 2/4. Earlier '63 detected' used a looser any-count criterion.
- Melanocyte analysis re-done without circularity (proxy excludes each panel's own genes; results_exploratory_melanocyte_v2.tsv): rho P02 -0.90, P06 -0.85, P05 -0.53, S04 -0.30; residual AUCs 0.68-0.86 treated as proxy limitation.

## Pass 5 (20 Nov 2026): external review of draft v1
- Fixed: Eq. 1 empty-set clause; precise per-compartment missing-sample statement (counts were correct, wording "two complete sets" was not); Phipson-Smyth p-value formula and reference; scored/total genes marked in Table I; inferred-direction and P02 notes in abstract and conclusion; Table III hyphenation; references renumbered in order of first citation.
- Checked and not found: missing tilde and gamma symbols (render correctly; reviewer likely read text-extraction output); Table II overflow; "track loss melanocytes" wording.
- Rejected: acceptance-odds percentages (no basis).

## Pass 6 (22 Nov 2026): post hoc checks requested in external reviews

1. Site-matched comparison, lesional vs perilesional skin of the same patient (23 pairs, VP-008 X and non-X kept as separate people as in the lock). P02 and P06 higher in lesional skin in 20/23 pairs (q < 0.001); P04, P09, S02 remain reversed (q <= 0.011); P05 16/23 (q = 0.13); S04 17/23 (q = 0.051); S03 and the IFN-gamma reference lower in lesional skin. Script `posthoc_site_matched.py`.
2. The 13 directions inferred from GSE65127 were checked in GSE75819 (GPL6884, 15 paired L/N biopsies; values pulled from the GEO series matrix for the listed probes). CACTIN (C19orf29) is not on GPL6884. 11 of 12 agree (9 with Wilcoxon p < 0.05). HADH disagrees (up in all 15 pairs in GSE75819; down in GSE65127, log2 -0.28; in GSE311770 L vs HC median logCPM 5.47 vs 5.72, p = 0.078). Pairs share an intensity batch, checked. Script `posthoc_direction_GSE75819.py`.
3. Rescoring with confirmed directions only: P05 L vs HC AUC 0.903 -> 0.927 (HADH and CACTIN removed), 0.907 with HADH flipped; P06 and S02 unchanged (all their checked directions agree). Script `posthoc_direction_rescore.py`.
4. Replication-rule sensitivity: CI-only 6 (L) / 7 (PL); locked rule 3 / 3; 10,000 nulls q < 0.05 3 / 3; 10,000 nulls q < 0.01 1 (S04) / 0. The statement in one external guide that results are "unchanged at q < 0.01" is false and was not used. Script `posthoc_rule_sensitivity.py`.
5. GEO records re-read for GSE75819, GSE65127, GSE53148 and GSE90880 (`verification/GEO_records_checked_2026-11-22.txt`). GSE75819 non-lesional samples are labelled "Normal Skin" in GEO; the paper now says so.

## Pass 7 (22 Nov 2026): second cohort, GSE298871 (post hoc)

1. GEO record re-read: 17 samples (6 healthy, 6 lesional, 5 non-lesional), platform GPL24676, 10x matrices only; the bulk RNA-seq of the same study was not deposited. Three healthy samples are named 465, 632, 658 and three CON161/171/194, which suggests two sources of healthy skin; noted, not adjusted.
2. Pseudobulk: UMIs summed over all cells per biopsy (200 to 3,485 cells); CPM > 1 in >= 4 samples kept 18,360 genes. Locked panels, locked null sets and the same rank score applied unchanged; AUC p exact.
3. Lesional vs healthy (6 vs 6): P05 0.889 (q 0.035), P06 0.972 (q 0.035), pigment reference 0.917, P04 0.167, P09 0.111, P02 0.833 (q 0.16), S04 0.639 (q 0.61), S02 0.500, S03 0.861 (q 0.16).
4. Composition varies widely (keratinocytes 1.5% to 61%). Excluding 658 (3.5% keratinocytes, mostly endothelium) and 486_NL (200 cells, 3 keratinocytes): P02 1.00 (q 0.014), P06 1.00, P05 0.90, P04 0.00, P09 0.00, P03 0.83 (q 0.049), S04 0.57. P03 passing only in the sensitivity analysis is not claimed in the paper.
5. S04: most S04 gene UMIs come from keratinocytes (59% in controls, 83% in lesional skin). S04 score within keratinocytes AUC 0.56 (p 0.82); fraction of cycling keratinocytes AUC 0.61 (p 0.59). The keratinocyte-proliferation explanation is not supported in this cohort; the paper now treats S04 as unconfirmed.
6. Abstract shortened to 249 words to make room for one sentence on the second cohort. Table II gained a column with the second-cohort AUC.
