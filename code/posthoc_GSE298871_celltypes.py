"""Post hoc, GSE298871: coarse marker-based cell typing and a test of whether panel S04
(cell-cycle genes shared with psoriasis) reflects keratinocyte proliferation.
python code/posthoc_GSE298871_celltypes.py <sc.pkl> <gene_sets_locked.json> <outdir>"""
import json, pickle, sys, numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
PKL, LOCK, OUT = sys.argv[1:4]
D = pickle.load(open(PKL, "rb")); L = json.load(open(LOCK)); S, feat = D["samples"], D["feat"]
sym = pd.Series(feat[1].values); ens = feat[0].values
MK = {"keratinocyte": ["KRT5","KRT14","KRT1","KRT10","KRT15","KRT6A"],
      "melanocyte": ["PMEL","MLANA","TYRP1","DCT","TYR"],
      "T/NK cell": ["CD3D","CD3E","CD2","TRAC","NKG7","GNLY"],
      "myeloid": ["LYZ","CD68","CD14","AIF1","CD1C","HLA-DQA1"],
      "fibroblast": ["COL1A1","COL1A2","DCN","LUM","PDGFRA"],
      "endothelial": ["PECAM1","VWF","CDH5","CLDN5"],
      "pericyte/SMC": ["ACTA2","RGS5","TAGLN","MYH11"],
      "sweat/secretory": ["SCGB2A2","DCD","KRT7","KRT19"],
      "mast": ["TPSAB1","CPA3","KIT"]}
idx = {t: [i for i in np.where(sym.isin(g))[0]] for t, g in MK.items()}
s04 = set(L["panels"]["S04"]["up"]); s04_i = np.where(pd.Series(ens).isin(s04))[0]
grp = lambda k: "lesional" if k.endswith("_LS") else "nonlesional" if k.endswith("_NL") else "control"
cycle_i = np.where(sym.isin(["MKI67", "TOP2A"]))[0]
comp, s04src, kc = [], [], []
for k, s in S.items():
    m = s["m"].tocsr(); tot = np.asarray(m.sum(0)).ravel()
    norm = m.multiply(1e4 / tot).tocsr(); norm.data = np.log1p(norm.data)
    sc = pd.DataFrame({t: np.asarray(norm[i].mean(0)).ravel() for t, i in idx.items()})
    lab = sc.idxmax(1).where(sc.max(1) > 0.5, "unassigned").values
    comp.append(pd.Series(lab).value_counts(normalize=True).rename(k))
    s04u = np.asarray(m[s04_i].sum(0)).ravel()
    s04src.append(pd.Series(s04u).groupby(lab).sum().rename(k) / max(s04u.sum(), 1))
    ker = lab == "keratinocyte"; cyc = np.asarray((m[cycle_i] > 0).sum(0)).ravel() > 0
    kc.append(dict(sample=k, group=grp(k), n_cells=len(lab), n_keratinocytes=int(ker.sum()),
                   frac_keratinocytes_cycling=round(float(cyc[ker].mean()), 4) if ker.any() else np.nan,
                   frac_Tcells_cycling=round(float(cyc[lab == "T/NK cell"].mean()), 4) if (lab == "T/NK cell").any() else np.nan,
                   kc_pseudobulk=np.asarray(m[:, ker].sum(1)).ravel()))
comp = pd.DataFrame(comp).fillna(0).round(3); comp.insert(0, "group", [grp(k) for k in comp.index])
src = pd.DataFrame(s04src).fillna(0).round(3); src.insert(0, "group", [grp(k) for k in src.index])
# S04 scored within keratinocytes only (keratinocyte pseudobulk, same rank score)
KP = pd.DataFrame({r["sample"]: r["kc_pseudobulk"] for r in kc}, index=ens).groupby(level=0).sum()
cpm = KP / KP.sum() * 1e6; keep = (cpm > 1).sum(1) >= 4; R = np.log2(cpm[keep] + 1).rank() / keep.sum()
s04_kc = R.loc[[g for g in s04 if g in R.index]].mean()
kc = pd.DataFrame([{k: v for k, v in r.items() if k != "kc_pseudobulk"} for r in kc]).set_index("sample")
kc["S04_score_keratinocyte_pseudobulk"] = s04_kc.round(4)
def mw(col, g):
    a = kc.loc[kc.group == g, col]; b = kc.loc[kc.group == "control", col]
    a, b = a.dropna(), b.dropna()
    if len(a) == 0 or len(b) == 0: return np.nan, np.nan
    u = stats.mannwhitneyu(a, b, method="exact"); return round(float(u.statistic) / (len(a) * len(b)), 3), round(float(u.pvalue), 4)
summ = pd.DataFrame([dict(measure=c, comparison=f"{g} vs control", AUC=mw(c, g)[0], p_exact=mw(c, g)[1]) for c in
        ["S04_score_keratinocyte_pseudobulk", "frac_keratinocytes_cycling", "frac_Tcells_cycling"] for g in ["lesional", "nonlesional"]])
comp.to_csv(f"{OUT}/posthoc_GSE298871_celltype_composition.tsv", sep="\t")
src.to_csv(f"{OUT}/posthoc_GSE298871_S04_UMI_source.tsv", sep="\t")
kc.to_csv(f"{OUT}/posthoc_GSE298871_S04_keratinocytes.tsv", sep="\t")
summ.to_csv(f"{OUT}/posthoc_GSE298871_S04_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(comp.sort_values("group").to_string()); print(src.sort_values("group").to_string()); print(kc.sort_values("group").to_string()); print(summ.to_string(index=False))
