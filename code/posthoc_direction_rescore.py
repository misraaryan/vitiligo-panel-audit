"""Post hoc: rescore the panels with inferred directions (P05, P06, S02) on GSE311770
using only genes whose GSE65127-inferred direction is confirmed in GSE75819
(see posthoc_direction_GSE75819.tsv). Descriptive only (AUC + individual bootstrap CI)."""
import json, numpy as np, pandas as pd
from scipy import stats
L = json.load(open("gene_sets_locked.json")); rng = np.random.default_rng(L["seed"])
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
R = logc.rank(axis=0) / len(logc)
def score(up, dn):
    up=[g for g in up if g in R.index]; dn=[g for g in dn if g in R.index]
    s = pd.Series(0.0, index=R.columns)
    if up: s += R.loc[up].mean()
    if dn: s -= R.loc[dn].mean()
    return s
def auc(a, b): u = stats.mannwhitneyu(a, b, alternative="two-sided"); return u.statistic/(len(a)*len(b)), u.pvalue
def ci(s, c1, c0):
    return np.percentile([auc(s[rng.choice(c1,len(c1))].values, s[rng.choice(c0,len(c0))].values)[0] for _ in range(2000)], [2.5,97.5])
P = L["panels"]; sym = {e:s for p in P.values() for e,s in p["symbols"].items()}; ens = {s:e for e,s in sym.items()}
chk = pd.read_csv("posthoc_direction_GSE75819.tsv", sep="\t")
drop = set(chk[(chk.panel!="ctrl") & (chk.agree==False)].gene)      # disagree in GSE75819 (HADH)
untested = set(chk[chk.probe=="not on GPL6884"].gene)                 # CACTIN
var = {}
for pid in ["P05","P06","S02"]:
    up, dn = P[pid]["up"], P[pid]["down"]
    var[f"{pid}_locked"] = (up, dn)
    var[f"{pid}_confirmed_only"] = ([g for g in up if sym[g] not in drop|untested], [g for g in dn if sym[g] not in drop|untested])
    if any(sym[g] in drop for g in up+dn):
        var[f"{pid}_disputed_flipped"] = ([g for g in up if sym[g] not in drop]+[g for g in dn if sym[g] in drop],
                                         [g for g in dn if sym[g] not in drop]+[g for g in up if sym[g] in drop])
rows=[]
for name,(up,dn) in var.items():
    s = score(up,dn)
    for comp,(g1,g0) in {"L_vs_HC":("lesional","control"),"PL_vs_HC":("perilesional","control")}.items():
        c1 = meta.index[meta.region==g1]; c0 = meta.index[meta.region==g0]
        a,p = auc(s[c1],s[c0]); lo,hi = ci(s,c1,c0)
        rows.append(dict(variant=name, comparison=comp, genes_up="|".join(sym[g] for g in up), genes_down="|".join(sym[g] for g in dn),
                         AUC=round(a,3), CI_lo=round(lo,3), CI_hi=round(hi,3), p_MW=float(f"{p:.2g}")))
o = pd.DataFrame(rows); o.to_csv("posthoc_direction_rescore.tsv", sep="\t", index=False); print(o.to_string(index=False))
# HADH in GSE311770 itself (descriptive)
h = ens["HADH"]; c1 = meta.index[meta.region=="lesional"]; c0 = meta.index[meta.region=="control"]
print("\nHADH GSE311770 L vs HC: median logCPM", round(logc.loc[h,c1].median(),2), "vs", round(logc.loc[h,c0].median(),2),
      "MW p", f"{stats.mannwhitneyu(logc.loc[h,c1],logc.loc[h,c0]).pvalue:.2g}")
