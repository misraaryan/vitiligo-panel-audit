"""E1c (post hoc): the E1 comparison with a small neural network. Multi-layer perceptron, one hidden
layer of 8 ReLU units, L2 penalty alpha=1.0 (fixed in advance, no tuning; strong because n = 51),
L-BFGS, up to 1,000 iterations. 5-fold stratified CV x 3 repeats; standardisation with training-fold
statistics only; first 200 locked null sets per panel. Same samples as E1."""
import json, warnings, numpy as np, pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
warnings.filterwarnings("ignore")
L = json.load(open("gene_sets_locked.json"))
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]
y = (meta.loc[s, "region"] == "lesional").astype(int).values; Xs = X[s].T
def cv_auc(genes):
    g = [x for x in genes if x in Xs.columns]; A = Xs[g].values; out = []
    for sd in [11, 22, 33]:
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=sd).split(A, y):
            mu = A[tr].mean(0); s_ = A[tr].std(0); s_[s_ == 0] = 1
            m = MLPClassifier((8,), alpha=1.0, solver="lbfgs", max_iter=1000, random_state=sd).fit((A[tr] - mu) / s_, y[tr])
            oof[te] = m.predict_proba((A[te] - mu) / s_)[:, 1]
        out.append(roc_auc_score(y, oof))
    return float(np.mean(out))
import os
PART = "ml_e1c_mlp.partial.tsv"
rows = pd.read_csv(PART, sep="\t").to_dict("records") if os.path.exists(PART) else []
done = {r["panel"] for r in rows}
for pid, p in L["panels"].items():
    if pid in done: continue
    obs = cv_auc(p["up"] + p["down"])
    null = np.array(Parallel(n_jobs=2)(delayed(cv_auc)(x) for x in L["null_sets"][pid]["sets"][:200]))
    rows.append(dict(model="mlp", panel=pid, cv_AUC=round(obs, 3), null_median=round(float(np.median(null)), 3),
                     n_null=200, p_emp=(1 + (null >= obs).sum()) / 201))
    print(rows[-1], flush=True); pd.DataFrame(rows).to_csv(PART, sep="\t", index=False)
o = pd.DataFrame(rows); o["q_BH"] = multipletests(o.p_emp, method="fdr_bh")[1].round(4)
o.to_csv("ml_e1c_mlp.tsv", sep="\t", index=False); print(o.to_string(index=False))
