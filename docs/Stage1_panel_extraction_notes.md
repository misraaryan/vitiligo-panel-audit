# Stage 1: panel extraction notes (20 Nov 2026)

Sources: PubMed / PubMed Central full texts (efetch), gene symbols mapped to Ensembl via the Ensembl REST API.

## Inclusion criteria (fixed before scoring)
- Human skin transcriptomic data from GEO used for discovery.
- Paper names a specific gene panel (or single gene) proposed as vitiligo biomarker / hub / diagnostic genes.
- Gene list extractable from the main text or tables.
- Discovery data must not include GSE311770.
- Tier "primary": vitiligo-only skin studies. Tier "secondary": comorbidity (shared-gene) studies and a blood-derived panel tested in skin.

## Result
- 9 primary panels, 5 secondary panels, 2 excluded (see panels_excluded.tsv).
- 71 gene entries (64 unique genes), all present in GSE311770; 63/64 unique genes detected (IFNG below detection, expected for bulk skin).

## Gene direction
Direction not reported for 13 gene entries (P05 all 7; P06 TYR, TYRP1, DCT, LARP7; S02 EXOSC7, KLRG1). P04 and P09 report directions that contradict their own discovery data (see Verification_log.md).
Pre-specified rule: infer direction from each paper's own discovery dataset (lesional vs control), never from GSE311770. P05 and P04 need GSE65127/GSE53146/GSE75819 for this.

## Methodological problems observed in the source literature (usable in Introduction/Discussion)
1. GSE75819 contains lesional and non-lesional skin from the same patients and no healthy controls, yet P02 and P07 describe the comparison as vitiligo vs "normal/healthy control".
2. P06 (Zhang 2021) validated a melanocyte panel (TYR, TYRP1, DCT) on GSE90880, a peripheral blood dataset, reporting AUC 0.942.
3. Excluded Lu 2025 used GSE53148 as its test/validation set. GSE53148 is a GEO SuperSeries combining human GSE53146 with mouse GSE53147; it is not a separate human cohort, and the paper does not state that mouse samples were excluded. (Not a duplication of Lu's own training data.)
4. Most claimed AUCs were computed on the training/discovery data; only P04, P05 and P06 report any external dataset, and all of those reuse the same three overlapping microarray cohorts.
5. 5 of 9 primary panels contain melanocyte/pigmentation genes; P09 consists entirely of them (the paper describes them as "upregulating" genes, although these genes are lost in lesional skin).

## Provisional pigment flag
TYR, TYRP1, DCT, MLANA, PMEL, OCA2, SLC45A2, GPR143, TRPM1, EDNRB, KIF1A, SNCA. Final melanogenesis vs melanocyte-identity lists are defined in stage 3 from curated sources, before outcome analysis.
