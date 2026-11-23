"""ML experiment E2 (post hoc): what does an L1 (LASSO) logistic regression pick when it may choose
from all 18,848 detected genes? Lesional vs healthy, GSE311770. Nested CV: outer 5-fold x 5 repeats;
inside each outer training fold, genes are standardised and the L1 strength is tuned by an inner
5-fold CV (10 values). No gene filtering or selection outside the training fold.
Output: outer-fold AUC, selection frequency of each gene, and overlap with the curated pigment set.
Gene symbols from the 10x features file of GSE298871 (Ensembl -> symbol)."""
import json, pickle, sys, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
SC = sys.argv[1]
L = json.load(open("gene_sets_locked.json"))
X = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
s = meta.index[meta.region.isin(["lesional", "control"])]
y = (meta.loc[s, "region"] == "lesional").astype(int).values
A = X[s].T.values; genes = X.index.values
feat = pickle.load(open(SC, "rb"))["feat"]; sym = dict(zip(feat[0], feat[1]))
pig = set(L["pigment"]["union"])
SEEDS = [11, 22, 33, 44, 55]
sel_count = pd.Series(0, index=genes); n_models = 0; aucs = []; nsel = []
for sd in SEEDS:
    oof = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=sd).split(A, y):
        m = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-2, 1, 10), penalty="l1", solver="liblinear",
                          cv=StratifiedKFold(5, shuffle=True, random_state=sd + 1), scoring="roc_auc", max_iter=5000))
        m.fit(A[tr], y[tr]); oof[te] = m.predict_proba(A[te])[:, 1]
        coef = m[-1].coef_.ravel(); nz = genes[coef != 0]
        sel_count[nz] += 1; n_models += 1; nsel.append(len(nz))
    aucs.append(roc_auc_score(y, oof))
freq = (sel_count / n_models).sort_values(ascending=False)
top = freq[freq > 0]
out = pd.DataFrame({"ensembl": top.index, "symbol": [sym.get(g, "") for g in top.index], "selection_freq": top.values.round(3),
                    "in_curated_pigment_set": [g in pig for g in top.index]})
out.to_csv("ml_e2_lasso_selected_genes.tsv", sep="\t", index=False)
summ = dict(outer_AUC_mean=round(float(np.mean(aucs)), 3), outer_AUC_by_repeat=[round(a, 3) for a in aucs],
            genes_selected_per_model_median=int(np.median(nsel)), n_models=n_models,
            genes_selected_ever=int((freq > 0).sum()), genes_selected_ge50pct=int((freq >= 0.5).sum()),
            pigment_among_ge50pct=int(sum(g in pig for g in freq[freq >= 0.5].index)))
json.dump(summ, open("ml_e2_lasso_summary.json", "w"), indent=1)
print(summ); print(out.head(30).to_string(index=False))
