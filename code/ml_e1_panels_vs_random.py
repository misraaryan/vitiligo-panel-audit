"""ML experiment E1 (post hoc): do classifiers trained on published panel genes beat
classifiers trained on matched random gene sets? Lesional vs healthy skin, GSE311770.
Features: log2 CPM of the panel's detected genes, standardised with training-fold mean and SD only.
Model: L2 logistic regression (C=1, fixed in advance, no tuning). Evaluation: 5-fold stratified CV,
5 repeats (fixed seeds); AUC from out-of-fold predictions, averaged over repeats. Each sample is a
different individual in this comparison, so sample-level folds cannot leak a person.
Random baseline: the 1,000 locked expression-matched null sets of each panel, same procedure."""
import json, numpy as np, pandas as pd, sys, time
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
L = json.load(open("gene_sets_locked.json"))
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]
y = (meta.loc[s, "region"] == "lesional").astype(int).values
Xs = X[s].T  # samples x genes
SEEDS = [11, 22, 33, 44, 55]
def cv_auc(genes):
    g = [x for x in genes if x in Xs.columns]
    if not g: return np.nan
    A = Xs[g].values; aucs = []
    for sd in SEEDS:
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=sd).split(A, y):
            mu = A[tr].mean(0); sd_ = A[tr].std(0); sd_[sd_ == 0] = 1
            m = LogisticRegression(C=1.0, max_iter=2000).fit((A[tr] - mu) / sd_, y[tr])
            oof[te] = m.decision_function((A[te] - mu) / sd_)
        aucs.append(roc_auc_score(y, oof))
    return float(np.mean(aucs))
rows = []
t0 = time.time()
for pid, p in L["panels"].items():
    genes = p["up"] + p["down"]
    obs = cv_auc(genes)
    ns = L["null_sets"][pid]["sets"]
    null = np.array(Parallel(n_jobs=2, batch_size=25)(delayed(cv_auc)(x) for x in ns))
    rows.append(dict(panel=pid, n_genes=sum(x in Xs.columns for x in genes), cv_AUC=round(obs, 3),
                     null_median=round(float(np.median(null)), 3), null_q95=round(float(np.quantile(null, 0.95)), 3),
                     frac_null_ge_obs=round(float((null >= obs).mean()), 3), p_emp=(1 + (null >= obs).sum()) / (1 + len(null))))
    print(pid, rows[-1], f"{time.time()-t0:.0f}s", flush=True)
o = pd.DataFrame(rows); o["q_BH"] = multipletests(o.p_emp, method="fdr_bh")[1].round(4)
o.to_csv("ml_e1_panels_vs_random.tsv", sep="\t", index=False); print(o.to_string(index=False))
