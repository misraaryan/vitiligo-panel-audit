"""Post hoc second cohort: GSE298871 (Brunner et al. 2025; 10x scRNA-seq; 6 healthy, 6 lesional,
5 non-lesional skin biopsies). Counts are summed per biopsy (pseudobulk) and the locked panels are
scored exactly as in GSE311770 (within-sample ranks over detected genes; no retraining), with the
locked expression-matched null sets. Run from the analysis directory:
python code/posthoc_GSE298871_pseudobulk.py <sc.pkl> <gene_sets_locked.json> <outdir> [comma-separated samples to exclude]"""
import json, pickle, sys, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
PKL, LOCK, OUT = sys.argv[1:4]; EXCL = sys.argv[4].split(",") if len(sys.argv) > 4 else []
D = pickle.load(open(PKL, "rb")); L = json.load(open(LOCK))
S, feat = {k: v for k, v in D["samples"].items() if k not in EXCL}, D["feat"]
TAG = "_excl_lowKC" if EXCL else ""
ens = feat[0].values
grp = {k: ("lesional" if k.endswith("_LS") else "nonlesional" if k.endswith("_NL") else "control") for k in S}
# pseudobulk: sum UMIs over all cells of each biopsy (duplicate Ensembl IDs summed)
pb = pd.DataFrame({k: np.asarray(s["m"].sum(1)).ravel() for k, s in S.items()}, index=ens).groupby(level=0).sum()
cpm = pb / pb.sum() * 1e6
keep = (cpm > 1).sum(1) >= 4                      # ~20% of samples, as in GSE311770 (20 of 105)
logc = np.log2(cpm[keep] + 1)
R = logc.rank(axis=0) / len(logc)
print("genes detected:", keep.sum(), "| samples:", pd.Series(grp).value_counts().to_dict())
def score(up, dn):
    up = [g for g in up if g in R.index]; dn = [g for g in dn if g in R.index]
    s = pd.Series(0.0, index=R.columns)
    if up: s += R.loc[up].mean()
    if dn: s -= R.loc[dn].mean()
    return s
def auc(a, b):
    u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact"); return u.statistic/(len(a)*len(b)), u.pvalue
cases = {"L_vs_HC": "lesional", "NL_vs_HC": "nonlesional"}
ctrl = [k for k in S if grp[k] == "control"]
pairs = sorted(k[:-3] for k in S if k.endswith("_LS") and k[:-3]+"_NL" in S)
sets = {pid: (p["up"], p["down"], L["null_sets"][pid]) for pid, p in L["panels"].items()}
sets["REF_pigment"] = ([], L["pigment"]["union"], None); sets["REF_IFNG"] = (L["ifng_hallmark"]["genes"], [], None)
rows, scores = [], {}
for name, (up, dn, ns) in sets.items():
    s = score(up, dn); scores[name] = s
    ndet = sum(g in R.index for g in up + dn)
    for comp, g in cases.items():
        case = [k for k in S if grp[k] == g]
        a, p = auc(s[case], s[ctrl])
        r = dict(set=name, comparison=comp, n_case=len(case), n_ctrl=len(ctrl), genes_detected=f"{ndet}/{len(up)+len(dn)}", AUC=round(a, 3), p_MW_exact=p)
        if ns is not None:
            na = np.array([auc(sc[case], sc[ctrl])[0] for sc in (score(x[:ns["n_up"]], x[ns["n_up"]:]) for x in ns["sets"])])
            r["null_AUC_median"] = round(float(np.median(na)), 3); r["p_null"] = (1 + (na >= a).sum()) / (1 + len(na))
        rows.append(r)
    d = np.array([s[k+"_LS"] - s[k+"_NL"] for k in pairs])
    rows.append(dict(set=name, comparison="L_vs_NL_paired", n_case=len(pairs), n_ctrl=len(pairs), genes_detected=f"{ndet}/{len(up)+len(dn)}",
                     AUC=round(float((d > 0).mean()), 3), p_MW_exact=stats.wilcoxon(d).pvalue if len(set(np.sign(d))) > 0 else np.nan))
o = pd.DataFrame(rows)
for comp in cases:
    m = (o.comparison == comp) & o.p_null.notna()
    o.loc[m, "q_null"] = multipletests(o.loc[m, "p_null"], method="fdr_bh")[1]
o.to_csv(f"{OUT}/posthoc_GSE298871_pseudobulk{TAG}.tsv", sep="\t", index=False)
pd.DataFrame(scores).assign(group=pd.Series(grp)).to_csv(f"{OUT}/posthoc_GSE298871_scores{TAG}.tsv", sep="\t")
pd.set_option("display.width", 200); print(o.to_string(index=False))
