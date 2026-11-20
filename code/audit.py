"""Independent audit: recompute from raw counts with separate code paths and check every reported claim."""
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

L = json.load(open("gene_sets_locked.json"))
res = pd.read_csv("results_all.tsv", sep="\t")
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)

# --- A. rebuild expression independently from raw counts
cnt = pd.read_csv("data/GSE311770_all_samples_rsem_gene_counts.tsv.gz", sep="\t", index_col=0)
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True); cnt = cnt.groupby(level=0).sum()
cpm = cnt.div(cnt.sum(axis=0), axis=1) * 1e6
keep = (cpm > 1).sum(axis=1) >= 20
lg = np.log2(cpm[keep] + 1)
old = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
print("A. rebuilt matrix identical to analysis matrix:", lg.shape == old.shape and np.allclose(lg.loc[old.index, old.columns].values, old.values))

# --- B. label sanity: every sample label parsed from library name
reg_from_name = meta.index.str.lower().map(lambda s: "control" if "contro" in s else "perilesional" if "perilesioneel" in s else "lesional" if "lesioneel" in s else "nonlesional" if "gezond" in s else "?")
print("B. region labels re-parsed consistently:", (reg_from_name == meta.region).all(), "| counts", meta.region.value_counts().to_dict())

# --- C. independent AUC with sklearn on an independently computed score (scipy rankdata)
R = pd.DataFrame(np.apply_along_axis(stats.rankdata, 0, lg.values) / len(lg), index=lg.index, columns=lg.columns)
def sc(up, dn):
    up = [g for g in up if g in R.index]; dn = [g for g in dn if g in R.index]
    return (R.loc[up].mean() if up else 0) - (R.loc[dn].mean() if dn else 0)
mx = 0
for pid, p in L["panels"].items():
    s = sc(p["up"], p["down"])
    for comp, g1 in [("L_vs_HC", "lesional"), ("PL_vs_HC", "perilesional"), ("NL_vs_HC", "nonlesional")]:
        m = meta.region.isin([g1, "control"])
        a = roc_auc_score((meta.region[m] == g1).astype(int), s[m.index[m]])
        b = res[(res.set == pid) & (res.comparison == comp)].AUC.iloc[0]
        mx = max(mx, abs(a - b))
print("C. max |AUC difference| between independent recomputation and reported:", round(mx, 6))

# --- D. sex call bimodality
y = lg.loc["ENSG00000129824"] if "ENSG00000129824" in lg.index else None
print("D. RPS4Y1 log2CPM: females max", round(y[y < 3].max(), 2), "| males min", round(y[y >= 3].min(), 2), "(clear gap = reliable sex call)")

# --- E. replication criteria: two-sided MW vs chance AND matched-null, both BH
r = res[res.set.isin(L["panels"].keys()) & res.comparison.isin(["L_vs_HC", "PL_vs_HC", "NL_vs_HC"])].copy()
for comp in r.comparison.unique():
    mm = r.comparison == comp
    r.loc[mm, "q_null_chk"] = multipletests(r.loc[mm, "p_null"], method="fdr_bh")[1]
def verdict(x):
    if x.AUC < 0.5 and x.CI_hi < 0.5: return "REVERSED"
    if x.AUC > 0.5 and x.CI_lo > 0.5 and x.q_null_chk < 0.05: return "REPLICATES (beats chance and matched random genes, FDR<0.05)"
    if x.AUC > 0.5 and x.CI_lo > 0.5 and x.p_null < 0.05: return "beats random genes only before FDR"
    if x.AUC > 0.5 and x.CI_lo > 0.5: return "beats chance only"
    return "FAILS"
r["verdict"] = r.apply(verdict, axis=1)
print("\nE. verdicts with FDR-corrected matched-null criterion")
print(r[["set", "comparison", "AUC", "CI_lo", "CI_hi", "p_null", "q_null_chk", "verdict"]].round(3).to_string(index=False))

# --- F. covariate adjustment WITHOUT separation problems: OLS score ~ group + batch + sex + sebaceous
meta["male"] = (y >= 3).astype(int); meta["seb"] = sc(L["sebaceous_covariate"], [])
F = []
for pid, p in L["panels"].items():
    s = sc(p["up"], p["down"])
    for comp, g1 in [("L_vs_HC", "lesional"), ("PL_vs_HC", "perilesional")]:
        m = meta.region.isin([g1, "control"])
        d = pd.DataFrame({"score": (s[m.index[m]] - s[m.index[m]].mean()) / s[m.index[m]].std(),
                          "case": (meta.region[m] == g1).astype(int), "batch": meta.batch[m], "male": meta.male[m], "seb": meta.seb[m]})
        f0 = smf.ols("score ~ case", d).fit(); f1 = smf.ols("score ~ case + C(batch) + male + seb", d).fit()
        F.append(dict(set=pid, comparison=comp, effect_unadj=f0.params["case"], effect_adj=f1.params["case"], p_adj=f1.pvalues["case"]))
F = pd.DataFrame(F)
F["q_adj"] = F.groupby("comparison").p_adj.transform(lambda p: multipletests(p, method="fdr_bh")[1])
F["sign_kept"] = np.sign(F.effect_unadj) == np.sign(F.effect_adj)
print("\nF. covariate-adjusted group effect (SD units); sign preserved for all:", F.sign_kept.all())
print(F.round(3).to_string(index=False))
F.to_csv("audit_adjusted_ols.tsv", sep="\t", index=False)

# --- G. second scoring method: mean z-score (direction-aware)
Z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
G = []
for pid, p in L["panels"].items():
    up = [g for g in p["up"] if g in Z.index]; dn = [g for g in p["down"] if g in Z.index]
    s = (Z.loc[up].mean() if up else 0) - (Z.loc[dn].mean() if dn else 0)
    for comp, g1 in [("L_vs_HC", "lesional"), ("PL_vs_HC", "perilesional")]:
        m = meta.region.isin([g1, "control"])
        G.append(dict(set=pid, comparison=comp, AUC_rank=res[(res.set == pid) & (res.comparison == comp)].AUC.iloc[0],
                      AUC_zscore=roc_auc_score((meta.region[m] == g1).astype(int), s[m.index[m]])))
G = pd.DataFrame(G); G["diff"] = G.AUC_zscore - G.AUC_rank
print("\nG. rank vs z-score scoring: max |diff| =", round(G["diff"].abs().max(), 3), "| same side of 0.5 for all:",
      ((G.AUC_rank - 0.5) * (G.AUC_zscore - 0.5) > 0).all())
print(G[G["diff"].abs() > 0.05].round(3).to_string(index=False))
G.to_csv("audit_scoring_method.tsv", sep="\t", index=False)
r.to_csv("audit_verdicts.tsv", sep="\t", index=False)

# --- H. claims about the IFN-gamma reference
x = res[res.set == "REF_IFNG_hallmark"][["comparison", "AUC", "CI_lo", "CI_hi", "p_MW"]]
print("\nH. IFN-gamma reference:\n", x.round(3).to_string(index=False))
