"""E2b (post hoc): are the genes chosen by genome-wide LASSO (E2) melanocyte genes? Melanocyte
specificity is defined independently in GSE298871 single-cell data (healthy and non-lesional biopsies,
which contain melanocytes): cells are typed with the same marker rule as posthoc_GSE298871_celltypes.py;
a gene is 'melanocyte-enriched' if its mean normalised expression in melanocytes is >= 4x that in all
other cells and it is detected in >= 20% of melanocytes. The marker genes used for typing are excluded
from the enrichment test to avoid circularity.
python ml_e2b_melanocyte_specificity.py <sc.pkl>"""
import pickle, sys, json, numpy as np, pandas as pd, scipy.sparse as sp
D = pickle.load(open(sys.argv[1], "rb")); S, feat = D["samples"], D["feat"]
sym = pd.Series(feat[1].values); ens = feat[0].values
MK = {"keratinocyte": ["KRT5","KRT14","KRT1","KRT10","KRT15","KRT6A"], "melanocyte": ["PMEL","MLANA","TYRP1","DCT","TYR"],
      "T/NK cell": ["CD3D","CD3E","CD2","TRAC","NKG7","GNLY"], "myeloid": ["LYZ","CD68","CD14","AIF1","CD1C","HLA-DQA1"],
      "fibroblast": ["COL1A1","COL1A2","DCN","LUM","PDGFRA"], "endothelial": ["PECAM1","VWF","CDH5","CLDN5"],
      "pericyte/SMC": ["ACTA2","RGS5","TAGLN","MYH11"], "sweat/secretory": ["SCGB2A2","DCD","KRT7","KRT19"], "mast": ["TPSAB1","CPA3","KIT"]}
idx = {t: np.where(sym.isin(g))[0] for t, g in MK.items()}
mel_sum = np.zeros(len(ens)); oth_sum = np.zeros(len(ens)); mel_det = np.zeros(len(ens)); n_mel = 0; n_oth = 0
for k, s in S.items():
    if k.endswith("_LS"): continue
    m = s["m"].tocsc().astype(float); tot = np.asarray(m.sum(0)).ravel()
    norm = m.multiply(1e4 / tot).tocsc(); lg = norm.copy(); lg.data = np.log1p(lg.data)
    sc = pd.DataFrame({t: np.asarray(lg[i].mean(0)).ravel() for t, i in idx.items()})
    lab = sc.idxmax(1).where(sc.max(1) > 0.5, "unassigned").values
    mel = lab == "melanocyte"; oth = (lab != "melanocyte") & (lab != "unassigned")
    mel_sum += np.asarray(norm[:, mel].sum(1)).ravel(); oth_sum += np.asarray(norm[:, oth].sum(1)).ravel()
    mel_det += np.asarray((m[:, mel] > 0).sum(1)).ravel(); n_mel += mel.sum(); n_oth += oth.sum()
mm = mel_sum / max(n_mel, 1); om = oth_sum / max(n_oth, 1)
spec = pd.DataFrame({"ensembl": ens, "symbol": sym.values, "mean_melanocyte": mm, "mean_other": om,
                     "fold": (mm + 1e-3) / (om + 1e-3), "det_frac_melanocyte": mel_det / max(n_mel, 1)}).groupby("ensembl").first()
markers = set(g for v in MK.values() for g in v)
spec["melanocyte_enriched"] = (spec.fold >= 4) & (spec.det_frac_melanocyte >= 0.2) & ~spec.symbol.isin(markers)
spec["is_typing_marker"] = spec.symbol.isin(markers)
lasso = pd.read_csv("ml_e2_lasso_selected_genes.tsv", sep="\t").set_index("ensembl")
lasso = lasso.join(spec[["fold", "det_frac_melanocyte", "melanocyte_enriched", "is_typing_marker"]], how="left")
lasso.to_csv("ml_e2_lasso_selected_genes_annotated.tsv", sep="\t")
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0, usecols=[0])
bg = spec.reindex(X.index)
bg_frac = float(bg.melanocyte_enriched.mean())
def frac(th):
    sub = lasso[lasso.selection_freq >= th]; sub2 = sub[~sub.is_typing_marker.fillna(False).astype(bool)]
    return len(sub), int(sub.is_typing_marker.fillna(False).sum()), int(sub2.melanocyte_enriched.fillna(False).sum()), len(sub2)
out = {"cells_melanocyte": int(n_mel), "cells_other": int(n_oth), "background_fraction_melanocyte_enriched_among_detected_genes": round(bg_frac, 4)}
for th in [0.5, 0.25, 0.04]:
    n, nm, ne, n2 = frac(th)
    out[f"freq>={th}"] = dict(n_genes=n, typing_markers=nm, melanocyte_enriched_non_marker=ne, of_non_marker=n2)
json.dump(out, open("ml_e2b_summary.json", "w"), indent=1); print(json.dumps(out, indent=1))
print(lasso.head(25).to_string())
