"""Post hoc (added 22 Nov 2026, not in the locked plan): site-matched comparison.
Lesional vs perilesional biopsies from the same patient (perilesional skin is taken next to the lesion,
so the pair shares body site, sex, age and patient). Same locked panels, same rank score.
Output: posthoc_site_matched.tsv"""
import json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

L = json.load(open("gene_sets_locked.json"))
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
meta["pid"] = meta.patient + np.where(meta.x.astype(bool), "_X", "")   # VP-008 sets treated as two people (as locked)
R = logc.rank(axis=0) / len(logc)

def score(up, dn):
    up = [g for g in up if g in R.index]; dn = [g for g in dn if g in R.index]
    s = pd.Series(0.0, index=R.columns)
    if up: s += R.loc[up].mean()
    if dn: s -= R.loc[dn].mean()
    return s

sets = {pid: (p["up"], p["down"]) for pid, p in L["panels"].items()}
sets["REF_pigment_union_down"] = ([], L["pigment"]["union"])
sets["REF_IFNG_hallmark"] = (L["ifng_hallmark"]["genes"], [])

les = meta[meta.region == "lesional"]; per = meta[meta.region == "perilesional"]
pairs = sorted(set(les.pid) & set(per.pid))
li = [les.index[les.pid == p][0] for p in pairs]; pi = [per.index[per.pid == p][0] for p in pairs]
rows = []
for name, (up, dn) in sets.items():
    s = score(up, dn)
    d = s[li].values - s[pi].values
    rows.append(dict(set=name, n_pairs=len(pairs), frac_lesional_higher=(d > 0).mean(),
                     median_diff=np.median(d), p_wilcoxon=stats.wilcoxon(d).pvalue))
df = pd.DataFrame(rows)
panels = ~df.set.str.startswith("REF")
df.loc[panels, "q_BH"] = multipletests(df.loc[panels, "p_wilcoxon"], method="fdr_bh")[1]
df.to_csv("posthoc_site_matched.tsv", sep="\t", index=False)
print("pairs:", len(pairs), pairs)
print(df.round(4).to_string(index=False))
