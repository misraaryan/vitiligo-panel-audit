"""Post hoc third cohort: OMIX691 (NGDC; Xu et al. 2022; 10x single-cell RNA-seq of skin, 5 healthy (H01-H05)
and 10 vitiligo (P01-P10) samples). Counts are summed per sample (pseudobulk) and analysed exactly as GSE298871:
(1) training-free scoring of the locked panels (within-sample ranks over detected genes; detected = CPM > 1 in
>= 20% of samples, i.e. 3 of 15) against the locked matched null sets, BH-FDR over the 14 panels;
(2) models trained once on GSE311770 (lesional vs healthy) and tested once here, as in ml_e3_cross_cohort.py.
Run from the analysis directory (needs data/, gene_sets_locked.json):
python code/posthoc_OMIX691.py <dir with H01..P10 Cell Ranger folders>"""
import json, sys, os, numpy as np, pandas as pd, scipy.io
from scipy import stats
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests
RAW = sys.argv[1]; L = json.load(open("gene_sets_locked.json"))
names = [f"H0{i}" for i in range(1, 6)] + [f"P{i:02d}" for i in range(1, 11)]
pb, ncell = {}, {}
for k in names:
    d = os.path.join(RAW, k, "outs", "filtered_gene_bc_matrices", "GRCh38")
    m = scipy.io.mmread(os.path.join(d, "matrix.mtx")).tocsr()
    g = pd.read_csv(os.path.join(d, "genes.tsv"), sep="\t", header=None)[0].values
    pb[k] = pd.Series(np.asarray(m.sum(1)).ravel(), index=g).groupby(level=0).sum(); ncell[k] = m.shape[1]
pb = pd.DataFrame(pb); y = np.array([1 if k.startswith("P") else 0 for k in pb.columns])
cpm = pb / pb.sum() * 1e6; det = (cpm > 1).sum(1) >= 3
logc = np.log2(cpm[det] + 1); R = logc.rank(axis=0) / len(logc)
print("cells per sample:", ncell, "| genes detected:", int(det.sum()))
def auc(score, yy):
    score = np.asarray(score); u = stats.mannwhitneyu(score[yy == 1], score[yy == 0], method="exact")
    return u.statistic / ((yy == 1).sum() * (yy == 0).sum()), u.pvalue
def score(up, dn):
    up = [x for x in up if x in R.index]; dn = [x for x in dn if x in R.index]; s = pd.Series(0.0, index=R.columns)
    if up: s += R.loc[up].mean()
    if dn: s -= R.loc[dn].mean()
    return s
rows = []
sets = {pid: (p["up"], p["down"], L["null_sets"][pid]) for pid, p in L["panels"].items()}
sets["REF_pigment"] = ([], L["pigment"]["union"], None); sets["REF_IFNG"] = (L["ifng_hallmark"]["genes"], [], None)
for name, (up, dn, ns) in sets.items():
    a, p = auc(score(up, dn).values, y); ndet = sum(x in R.index for x in up + dn)
    r = dict(set=name, n_case=int(y.sum()), n_ctrl=int((y == 0).sum()), genes_detected=f"{ndet}/{len(up)+len(dn)}", AUC=round(a, 3), p_MW_exact=p)
    if ns is not None:
        na = np.array([auc(score(x[:ns["n_up"]], x[ns["n_up"]:]).values, y)[0] for x in ns["sets"]])
        r["null_AUC_median"] = round(float(np.median(na)), 3); r["p_null"] = (1 + (na >= a).sum()) / (1 + len(na))
    rows.append(r)
o = pd.DataFrame(rows); m_ = o.p_null.notna(); o.loc[m_, "q_null"] = multipletests(o.loc[m_, "p_null"], method="fdr_bh")[1]
o.to_csv("posthoc_OMIX691_pseudobulk.tsv", sep="\t", index=False); print(o.to_string(index=False))
# trained on GSE311770, tested once here
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0); meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]; ytr = (meta.loc[s, "region"] == "lesional").astype(int).values
Rtr = (X[s].rank() / len(X)).T; Rte = R.T; common = Rtr.columns.intersection(Rte.columns)
def fit_test(genes):
    g = [x for x in genes if x in common]
    if not g: return np.nan, 0
    mdl = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)).fit(Rtr[g].values, ytr)
    return auc(mdl.predict_proba(Rte[g].values)[:, 1], y)[0], len(g)
tf = o.set_index("set"); rows = []
for pid, p in L["panels"].items():
    a, n = fit_test(p["up"] + p["down"])
    null = np.array([fit_test(x)[0] for x in L["null_sets"][pid]["sets"]]); null = null[~np.isnan(null)]
    rows.append(dict(panel=pid, genes_used=n, trained_LR_test_AUC=round(a, 3), null_median=round(float(np.median(null)), 3),
                     p_emp=(1 + (null >= a).sum()) / (1 + len(null)), training_free_test_AUC=tf.loc[pid, "AUC"]))
e = pd.DataFrame(rows); e["q_BH"] = multipletests(e.p_emp, method="fdr_bh")[1].round(4)
mdl = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-2, 1, 10), penalty="l1", solver="liblinear",
      cv=StratifiedKFold(5, shuffle=True, random_state=7), scoring="roc_auc", max_iter=5000)).fit(Rtr[common].values, ytr)
ga, gp = auc(mdl.predict_proba(Rte[common].values)[:, 1], y)
e.to_csv("ml_e3_OMIX691.tsv", sep="\t", index=False)
summ = dict(test_n=int(len(y)), test_vitiligo=int(y.sum()), genes_detected=int(det.sum()), common_genes=int(len(common)),
            genomewide_L1_test_AUC=round(float(ga), 3), genomewide_L1_test_p=float(f"{gp:.3g}"), cells=ncell)
json.dump(summ, open("ml_e3_OMIX691_summary.json", "w"), indent=1); print(summ); print(e.to_string(index=False))
