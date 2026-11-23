"""ML experiment E3 (post hoc): train on GSE311770 (lesional vs healthy, 24 vs 27), test once on the
independent GSE298871 pseudobulk (6 vs 6). Features: within-sample scaled ranks (rank / n detected genes,
computed in each dataset separately), restricted to genes detected in both datasets, so no scaling
information passes between cohorts. Models: (a) L2 logistic regression (C=1) on each panel's genes;
(b) the same on each of the panel's 1,000 locked matched random sets (one fit each) -> empirical p;
(c) genome-wide L1 logistic regression with C tuned by 5-fold CV on the training cohort only.
Compared with the training-free score of the same panel in GSE298871 (posthoc_GSE298871_pseudobulk.tsv).
python ml_e3_cross_cohort.py <sc.pkl>"""
import json, pickle, sys, numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests
L = json.load(open("gene_sets_locked.json"))
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]
ytr = (meta.loc[s, "region"] == "lesional").astype(int).values
D = pickle.load(open(sys.argv[1], "rb")); S, feat = D["samples"], D["feat"]
keep_s = [k for k in S if not k.endswith("_NL")]
pb = pd.DataFrame({k: np.asarray(S[k]["m"].sum(1)).ravel() for k in keep_s}, index=feat[0].values).groupby(level=0).sum()
cpm = pb / pb.sum() * 1e6; det = (cpm > 1).sum(1) >= 4          # same filter as the pseudobulk analysis (4 of the samples)
Xte = np.log2(cpm[det] + 1); yte = np.array([1 if k.endswith("_LS") else 0 for k in Xte.columns])
Rtr = (X[s].rank() / len(X)).T; Rte = (Xte.rank() / len(Xte)).T     # samples x genes, within-sample scaled ranks
common = Rtr.columns.intersection(Rte.columns)
def auc(score, y):
    u = stats.mannwhitneyu(score[y == 1], score[y == 0], method="exact"); return u.statistic / ((y == 1).sum() * (y == 0).sum()), u.pvalue
def fit_test(genes):
    g = [x for x in genes if x in common]
    if not g: return np.nan, 0
    m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)).fit(Rtr[g].values, ytr)
    return auc(m.predict_proba(Rte[g].values)[:, 1], yte)[0], len(g)
tf = pd.read_csv("posthoc_GSE298871_pseudobulk.tsv", sep="\t"); tf = tf[tf.comparison == "L_vs_HC"].set_index("set")
rows = []
for pid, p in L["panels"].items():
    a, n = fit_test(p["up"] + p["down"])
    null = np.array([fit_test(x)[0] for x in L["null_sets"][pid]["sets"]]); null = null[~np.isnan(null)]
    rows.append(dict(panel=pid, genes_used=n, trained_LR_test_AUC=round(a, 3), null_median=round(float(np.median(null)), 3),
                     p_emp=(1 + (null >= a).sum()) / (1 + len(null)), training_free_test_AUC=tf.loc[pid, "AUC"]))
o = pd.DataFrame(rows); o["q_BH"] = multipletests(o.p_emp, method="fdr_bh")[1].round(4)
# genome-wide L1, tuned on training cohort only
A = Rtr[common].values
m = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-2, 1, 10), penalty="l1", solver="liblinear",
                  cv=StratifiedKFold(5, shuffle=True, random_state=7), scoring="roc_auc", max_iter=5000)).fit(A, ytr)
ga, gp = auc(m.predict_proba(Rte[common].values)[:, 1], yte); nsel = int((m[-1].coef_ != 0).sum())
o.to_csv("ml_e3_cross_cohort.tsv", sep="\t", index=False)
summ = dict(train_n=len(ytr), test_n=int(len(yte)), test_lesional=int(yte.sum()), common_genes=int(len(common)),
            genomewide_L1_test_AUC=round(float(ga), 3), genomewide_L1_test_p=float(f"{gp:.3g}"), genomewide_L1_genes_selected=nsel)
json.dump(summ, open("ml_e3_summary.json", "w"), indent=1); print(summ); print(o.to_string(index=False))
