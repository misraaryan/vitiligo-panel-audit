import pandas as pd, numpy as np

ENS = dict(x.split("=") for x in """GABARAPL2=ENSG00000034713,USP8=ENSG00000138592,RELA=ENSG00000173039,SP1=ENSG00000185591,TBC1D17=ENSG00000104946,DCT=ENSG00000080166,ALOX5=ENSG00000012779,SNCA=ENSG00000145335,SLC1A4=ENSG00000115902,IL33=ENSG00000137033,KIF1A=ENSG00000130294,CACTIN=ENSG00000105298,DCTN1=ENSG00000204843,GPR143=ENSG00000101850,HADH=ENSG00000138796,MRPL47=ENSG00000136522,NKTR=ENSG00000114857,NUF2=ENSG00000143228,TYR=ENSG00000077498,TYRP1=ENSG00000107165,LARP7=ENSG00000174720,PTEN=ENSG00000171862,TP53=ENSG00000141510,BCL2L1=ENSG00000171552,MAPK3=ENSG00000102882,ATG3=ENSG00000144848,ATG9A=ENSG00000198925,NFE2L2=ENSG00000116044,CCL2=ENSG00000108691,RB1CC1=ENSG00000023287,SLC1A5=ENSG00000105281,PARK7=ENSG00000116288,MLANA=ENSG00000120215,PMEL=ENSG00000185664,OCA2=ENSG00000104044,SLC45A2=ENSG00000164175,TRPM1=ENSG00000134160,EDNRB=ENSG00000136160,PTGS2=ENSG00000073756,CCL5=ENSG00000271503,PRDX4=ENSG00000123131,EXOSC7=ENSG00000075914,KLRG1=ENSG00000139187,MAPK14=ENSG00000112062,IFNG=ENSG00000111537,STAT1=ENSG00000115415,IL1B=ENSG00000125538,CXCL10=ENSG00000169245,BUB1=ENSG00000169679,CEP55=ENSG00000138180,CDK1=ENSG00000170312,TOP2A=ENSG00000131747,CENPF=ENSG00000117724,PBK=ENSG00000168078,MELK=ENSG00000165304,CCNB2=ENSG00000157456,MAD2L1=ENSG00000164109,NUSAP1=ENSG00000137804,TTK=ENSG00000112742,NEK2=ENSG00000117650,CDKN3=ENSG00000100526,PTTG1=ENSG00000164611,HMGA1=ENSG00000137309,PSMD13=ENSG00000185627""".split(","))

# panel_id, tier, first author/year, PMID, DOI, theme, discovery data, validation claimed, claimed AUC, genes "SYM:dir"
P = [
 ("P01","primary","Luo 2023","37287971","10.3389/fimmu.2023.1164124","Mitophagy","GSE53146+GSE75819 merged","qPCR only","0.89-0.99 per gene (training data)",
  "GABARAPL2:up USP8:up RELA:down SP1:down TBC1D17:down"),
 ("P02","primary","Yang 2025","40046940","10.2147/CCID.S496781","Oxidative stress","GSE75819 (lesional vs non-lesional; no healthy controls)","none","0.996 (training data)",
  "DCT:down"),
 ("P03","primary","Liu 2025","39814853","10.1038/s41598-025-86061-4","Ferroptosis","GSE65127+GSE53146 merged","none","model 1.000; genes 0.74-0.85 (training data)",
  "ALOX5:up SNCA:down SLC1A4:down IL33:down"),
 ("P04","primary","Zhao 2023","38084289","10.2147/CCID.S443165","ML (LASSO+SVM-RFE)","GSE53146+GSE75819 merged","GSE65127","0.998",
  "DCT:up KIF1A:up"),
 ("P05","primary","Lei 2020","32756109","10.1097/MD.0000000000021297","WGCNA hub genes","GSE65127","GSE75819","reported in figures only",
  "CACTIN:NR DCTN1:NR GPR143:NR HADH:NR MRPL47:NR NKTR:NR NUF2:NR"),
 ("P06","primary","Zhang 2021","34107850","10.1080/21655979.2021.1933743","ML (WGCNA+LASSO+SVM-RFE+RF)","GSE53146+GSE65127+GSE75819 (RRA)","GSE90880 (blood)","0.942 (blood validation)",
  "TYR:NR TYRP1:NR DCT:NR LARP7:NR"),
 ("P07","primary","Zhao 2024","38841253","10.4103/ijd.ijd_655_23","Autophagy hub genes","GSE75819 (lesional vs non-lesional)","qPCR only","not reported",
  "PTEN:up TP53:down BCL2L1:down MAPK3:down ATG3:up ATG9A:down RELA:down NFE2L2:up CCL2:up RB1CC1:up"),
 ("P08","primary","Yao 2026","42488926","10.2147/CCID.S617044","Programmed cell death + SMR","GSE65127+GSE53146+GSE75819 merged","none","SLC1A5 0.733",
  "SLC1A5:down PARK7:up"),
 ("P09","primary","AbdElneam 2024","39345987","10.2147/CCID.S480990","PPI hub genes","GSE65127 (GEO2R, healthy vs lesional)","none","not reported",
  "TYR:up MLANA:up TYRP1:up PMEL:up OCA2:up SLC45A2:up GPR143:up DCT:up TRPM1:up EDNRB:up"),
 ("S01","secondary","Zeng 2025","41158793","10.2147/CCID.S545677","Shared with periodontitis","GEO vitiligo + periodontitis","none","PTGS2 ~0.8",
  "PTGS2:up CCL5:up PRDX4:up"),
 ("S02","secondary","Xiang 2026","42418414","10.1371/journal.pone.0352672","Shared with depression","GSE65127+GSE53146 (vitiligo) + GSE98793 (MDD)","GSE52790, GSE80009","0.76-0.80",
  "EXOSC7:NR KLRG1:NR MAPK14:up"),
 ("S03","secondary","Lu 2024","38327551","10.2147/CCID.S451080","Shared with Hashimoto's","GEO vitiligo + HT","qPCR (blood)","IFNG 0.86, IL1B 0.80, STAT1 0.89, CXCL10 1.00",
  "IFNG:up STAT1:up IL1B:up CXCL10:up"),
 ("S04","secondary","Wang 2024","39258216","10.2147/CCID.S470149","Shared with psoriasis","GSE75819 + GSE30999","GSE65127","not reported",
  "BUB1:up CEP55:up CDK1:up TOP2A:up CENPF:up PBK:up MELK:up CCNB2:up MAD2L1:up NUSAP1:up TTK:up NEK2:up CDKN3:up PTTG1:up"),
 ("S05","secondary","Wang 2025","40443672","10.3389/fimmu.2025.1543355","Blood-derived, tested in skin","blood microarrays","GSE53146, GSE75819","not stated for skin",
  "HMGA1:down PSMD13:down"),
]
EXCLUDED = [
 ("Lu 2025","40948579","10.1155/drp/6672081","30-gene ANN signature listed only in a figure; gene list not extractable. Also tested on GSE53148, which duplicates GSE53146 human samples."),
 ("Xia 2026","-","Electronic Journal of Biotechnology","Full text not accessible here; discovery dataset unknown, so overlap with GSE311770 cannot be excluded. Revisit if full text obtained."),
]

logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
cnt = pd.read_csv("data/GSE311770_all_samples_rsem_gene_counts.tsv.gz", sep="\t", index_col=0)
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True); cnt = cnt.groupby(level=0).sum()
cpm = cnt / cnt.sum() * 1e6
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)

PIGMENT = {"TYR","TYRP1","DCT","MLANA","PMEL","OCA2","SLC45A2","GPR143","TRPM1","EDNRB","KIF1A","SNCA"}  # provisional flag only; final lists defined in stage 3

rows = []
for pid, tier, au, pmid, doi, theme, disc, val, auc, genes in P:
    for g in genes.split():
        sym, d = g.split(":")
        e = ENS.get(sym)
        present = e in cpm.index
        detected = present and ((cpm.loc[e] > 1).sum() >= 10)
        src = "inferred_pending" if d == "NR" else "reported"
        rows.append(dict(panel=pid, direction_source=src, tier=tier, study=au, pmid=pmid, doi=doi, theme=theme,
                         discovery=disc, validation_claimed=val, auc_claimed=auc, gene=sym,
                         ensembl=e, direction=d, in_matrix=present, detected=detected,
                         median_cpm=round(float(cpm.loc[e].median()), 2) if present else None,
                         pigment_flag=sym in PIGMENT))
gt = pd.DataFrame(rows)
gt.to_csv("panels_gene_level.tsv", sep="\t", index=False)

summ = gt.groupby(["panel","tier","study","theme","discovery","validation_claimed","auc_claimed"], sort=False).agg(
    n_genes=("gene","size"), n_detected=("detected","sum"),
    n_direction_reported=("direction", lambda s: (s != "NR").sum()),
    n_pigment=("pigment_flag","sum"), genes=("gene", lambda s: ", ".join(s))).reset_index()
summ.to_csv("panels_summary.tsv", sep="\t", index=False)
pd.DataFrame(EXCLUDED, columns=["study","pmid","doi","reason"]).to_csv("panels_excluded.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(summ[["panel","study","theme","n_genes","n_detected","n_direction_reported","n_pigment","discovery"]].to_string(index=False))
print("\nundetected genes:", list(gt[~gt.detected].gene))
