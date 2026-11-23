"""E1b (post hoc): is the E1 result model-specific? Panel genes vs matched random sets with
(a) random forest (200 trees, 5-fold CV x 2 repeats; first 100 locked null sets) and
(b) RBF support vector machine (C=1, gamma='scale'; 5-fold CV x 3 repeats; first 200 null sets).
Standardisation uses training-fold statistics only. Same samples as E1."""
import json, numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
L = json.load(open("gene_sets_locked.json"))
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]
y = (meta.loc[s, "region"] == "lesional").astype(int).values; Xs = X[s].T
def cv_auc(genes, kind):
    g = [x for x in genes if x in Xs.columns]; A = Xs[g].values; out = []
    seeds = [11, 22] if kind == "rf" else [11, 22, 33]
    for sd in seeds:
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=sd).split(A, y):
            mu = A[tr].mean(0); s_ = A[tr].std(0); s_[s_ == 0] = 1
            if kind == "rf":
                m = RandomForestClassifier(200, random_state=sd, n_jobs=1).fit(A[tr], y[tr]); oof[te] = m.predict_proba(A[te])[:, 1]
            else:
                m = SVC(C=1.0, gamma="scale").fit((A[tr] - mu) / s_, y[tr]); oof[te] = m.decision_function((A[te] - mu) / s_)
        out.append(roc_auc_score(y, oof))
    return float(np.mean(out))
rows = []
for kind, nn in [("svm_rbf", 200), ("rf", 100)]:
    for pid, p in L["panels"].items():
        obs = cv_auc(p["up"] + p["down"], kind)
        null = np.array(Parallel(n_jobs=2)(delayed(cv_auc)(x, kind) for x in L["null_sets"][pid]["sets"][:nn]))
        rows.append(dict(model=kind, panel=pid, cv_AUC=round(obs, 3), null_median=round(float(np.median(null)), 3),
                         n_null=nn, p_emp=(1 + (null >= obs).sum()) / (1 + nn)))
        print(rows[-1], flush=True)
o = pd.DataFrame(rows)
o["q_BH"] = o.groupby("model").p_emp.transform(lambda p: multipletests(p, method="fdr_bh")[1]).round(4)
o.to_csv("ml_e1b_other_models.tsv", sep="\t", index=False); print(o.to_string(index=False))
