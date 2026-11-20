"""Stage 6 (EXPLORATORY, not pre-specified): how much of each panel's discrimination is explained by melanocyte content?"""
import json, numpy as np, pandas as pd
from scipy import stats

L = json.load(open("gene_sets_locked.json"))
sc = pd.read_csv("scores_per_sample.tsv", sep="\t", index_col=0)
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
R = logc.rank(axis=0) / len(logc)
pig = -sc["REF_pigment_union_down"]          # higher = more melanocyte content

def auc(a, b):
    return stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))

rows = []
for pid, p in L["panels"].items():
    s = sc[pid]
    r_pig = stats.spearmanr(s, pig).correlation
    # residual after regressing panel score on melanocyte content (all 105 samples)
    b = np.polyfit(pig, s, 1); resid = s - np.polyval(b, pig)
    out = dict(panel=pid, study=p["study"], rho_with_melanocyte=r_pig)
    for comp, g1 in [("L_vs_HC", "lesional"), ("PL_vs_HC", "perilesional")]:
        case, ctrl = sc.index[sc.region == g1], sc.index[sc.region == "control"]
        out[f"{comp}_AUC"] = auc(s[case], s[ctrl])
        out[f"{comp}_AUC_resid"] = auc(resid[case], resid[ctrl])
    # leave-one-gene-out: which gene drives the L_vs_HC AUC?
    p = dict(p); p["up"] = [x for x in p["up"] if x in R.index]; p["down"] = [x for x in p["down"] if x in R.index]
    genes = p["up"] + p["down"]
    if len(genes) > 1:
        case, ctrl = sc.index[sc.region == "lesional"], sc.index[sc.region == "control"]
        drops = {}
        for gname in genes:
            up = [x for x in p["up"] if x != gname]; dn = [x for x in p["down"] if x != gname]
            s2 = (R.loc[up].mean() if up else 0) - (R.loc[dn].mean() if dn else 0)
            drops[p["symbols"][gname]] = out["L_vs_HC_AUC"] - auc(s2[case], s2[ctrl])
        top = max(drops, key=lambda k: abs(drops[k]))
        out["driver_gene"] = top; out["AUC_change_without_driver"] = -drops[top]
    rows.append(out)
d = pd.DataFrame(rows)
d.to_csv("results_exploratory_melanocyte.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(d.round(3).to_string(index=False))
