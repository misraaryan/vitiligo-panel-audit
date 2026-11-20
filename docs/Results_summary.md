# Results summary (audited version, 20 Nov 2026)

Scoring: direction-aware within-sample rank score, as-reported gene directions, locked gene sets (SHA-256 09addbb4...). AUC > 0.5 = panel is higher in vitiligo skin, as the authors claimed. 95% CI: 2,000 bootstrap resamples of individuals. p_null: empirical p against 1,000 expression-matched random gene sets; q_null: Benjamini-Hochberg across the 14 panels within each comparison.

Replication criterion (stated before writing; the locked file fixed BH correction but not which test): AUC CI above 0.5 AND q_null < 0.05. Panels that pass before but not after FDR are reported separately.

## Lesional vs healthy (24 vs 27): primary
| Verdict | Panels (AUC, 95% CI) |
|---|---|
| Replicates (FDR) | P06 Zhang 2021 0.96 (0.90-1.00); P05 Lei 2020 0.90 (0.80-0.98); S04 Wang 2024 0.82 (0.69-0.93) |
| Beats random genes before FDR only | P02 Yang 2025 0.95 (0.89-1.00, q_null 0.073); P07 Zhao 2024 0.76 (0.62-0.88, q_null 0.084) |
| Above chance only (not better than random genes; not robust to scoring method) | P08 Yao 2026 0.73 |
| Fails | P01, P03, S01, S03, S05 |
| Reversed | P09 AbdElneam 2024 0.03; P04 Zhao 2023 0.10; S02 Xiang 2026 0.20 (2 of 3 directions inferred) |

## Perilesional vs healthy (27 vs 27): co-primary
Replicates (FDR): P05 0.83, P06 0.86, S03 (STAT1 + CXCL10) 0.76. Before FDR only: P02 0.84, P07 0.69, S01 0.68, S04 0.67. Reversed: P09 0.10, P04 0.23.

## Secondary (exploratory)
- Non-lesional vs healthy: no panel passes after FDR; P02, P05, P06, S03 pass before FDR (AUC 0.71-0.77).
- Lesional vs non-lesional within the same 24 patients: P02 and P06 higher in lesional in 92% of patients; P05, P07, S04 in 83%. IFN-gamma hallmark score is lower in lesional than in the patient's own non-lesional skin in 79% of patients (Wilcoxon p = 0.005).

## Robustness checks (audit)
- Independent recomputation from raw counts (separate ranking and AUC code): identical to 6 decimals.
- Sex call: clean gap (RPS4Y1 females <= 0.51, males >= 6.09 log2CPM).
- Covariate adjustment (linear model: score ~ group + batch + sex + sebaceous score): every replicating and reversed panel keeps its sign and stays significant after FDR. This replaces the earlier logistic models, which were unusable for pigment panels (perfect separation).
- Second scoring method (mean z-score): all replicating and reversed verdicts unchanged. Method-sensitive: P08 lesional AUC 0.73 -> 0.52 and P01 0.42 -> 0.65; neither replicates under either method.
- No second independent bulk cohort exists: GSE298871 on GEO is single-cell only, and older cohorts (GSE65127, GSE53146, GSE75819) were used to build the panels.
- Within each sequencing batch and after excluding high-sebaceous samples, the replicating and reversed panels keep their direction (results_all.tsv).

## Melanocyte dependence (EXPLORATORY, added after results)
- P02 and P06 scores track melanocyte content almost perfectly (Spearman rho -0.92, -0.93). P05's strongest gene is GPR143 (melanosome gene, absent from the curated pigment lists); removing it lowers P05's lesional AUC by 0.20.
- After regressing out melanocyte content, lesional AUCs fall: P02 0.95 -> 0.74, P06 0.96 -> 0.76, P05 0.90 -> 0.68.
- Pre-specified curated-pigment removal could only be applied to P09 (other pigment panels keep fewer than 2 genes).

## What the data support, and what they do not
Supported:
- 3 of 14 panels replicate in lesional skin after FDR; 2 more replicate before FDR. The strongest are dominated by melanocyte-content genes.
- 3 panels discriminate only in the direction opposite to their published claim; for P04 and P09 this matches their own discovery data (GSE65127).
- An interferon-related panel (S03) replicates in perilesional skin but not lesional skin; the IFN-gamma hallmark is lower in lesional than paired non-lesional skin.

Not supported (corrected from the first draft of this summary):
- "Lesional skin is immunologically quieter than healthy skin": IFN-gamma hallmark AUC 0.40, CI 0.24-0.55 includes 0.5. Only the paired lesional vs non-lesional difference is significant.
- "IFN-gamma hallmark separates perilesional skin": AUC 0.63, CI 0.47-0.79 includes 0.5.
- "5 of 14 panels replicate": true only before FDR.
