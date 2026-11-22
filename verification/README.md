# Verification evidence

Files backing the factual claims in the paper and the post hoc checks added on 22 Nov 2026.

| Claim in the paper | Evidence |
|---|---|
| GSE75819 has 30 samples from 15 patients (lesional and non-lesional), no healthy donors; GEO labels the non-lesional samples "Normal Skin" | `GEO_records_checked_2026-11-22.txt` |
| GSE90880 (used by P06 for validation) is blood, not skin | `GEO_records_checked_2026-11-22.txt` |
| GSE53148 is a SuperSeries (GSE53146 + GSE53147) | `GEO_records_checked_2026-11-22.txt` |
| P04 (DCT, KIF1A up) and P09 (ten pigment genes up) contradict GSE65127 | `direction_check_GSE65127.tsv`, made by `../infer_directions.py` |
| 13 inferred directions: 11 of 12 testable genes agree in GSE75819; HADH disagrees; CACTIN not on GPL6884 | `GSE75819_selected_probes.tsv` (GEO series matrix values for the listed probes), `posthoc_direction_GSE75819.tsv`, made by `../posthoc_direction_GSE75819.py` |
| P05 AUC 0.90 to 0.93 without HADH; P06 and S02 unchanged | `posthoc_direction_rescore.tsv`, made by `../posthoc_direction_rescore.py` |
| Site-matched lesional vs perilesional comparison (23 pairs) | `posthoc_site_matched.tsv`, made by `../posthoc_site_matched.py` |
| Replication count under alternative rules (6 / 3 / 3 / 1 lesional) | `posthoc_rule_sensitivity.tsv`, made by `../posthoc_rule_sensitivity.py` |
| Second cohort GSE298871: 17 single-cell samples (6 healthy, 6 lesional, 5 non-lesional), bulk data not deposited | GEO record GSE298871 (platform GPL24676); results in `posthoc_GSE298871_*.tsv`, made by `posthoc_GSE298871_pseudobulk.py` and `posthoc_GSE298871_celltypes.py` |

GSE75819 notes: values are BeadStudio average-normalised after background subtraction and can be negative, and the series has two intensity batches. Every lesional and non-lesional pair sits in the same batch, so only within-pair differences of log2(max(x, 0) + 16) are used. For genes with several probes, the probe with the highest mean signal is used.

All analyses in this folder were run after the locked gene-set file was scored and are reported as post hoc.
