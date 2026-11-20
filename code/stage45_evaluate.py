"""Stages 4-5: score locked panels on GSE311770 and evaluate. Reads only the locked gene-set file."""
import json, hashlib, numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

raw = open("gene_sets_locked.json").read()
assert hashlib.sha256(raw.encode()).hexdigest() == open("gene_sets_locked.sha256").read().split()[0], "lock changed"
L = json.loads(raw)
rng = np.random.default_rng(L["seed"])
NBOOT = 2000

logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
meta["pid"] = np.where(meta.x.fillna(False).astype(bool), meta.patient + "X", meta.patient)  # VP-008 X / non-X = two people
meta.loc[meta.region == "control", "pid"] = meta.patient

# per-sample within-sample ranks (0..1) over all detected genes
R = logc.rank(axis=0) / len(logc)

def score(up, down):
    up = [g for g in up if g in R.index]; down = [g for g in down if g in R.index]
    s = pd.Series(0.0, index=R.columns)
    if up: s += R.loc[up].mean()
    if down: s -= R.loc[down].mean()
    return s

def auc(x_case, x_ctrl):
    u = stats.mannwhitneyu(x_case, x_ctrl, alternative="two-sided")
    return u.statistic / (len(x_case) * len(x_ctrl)), u.pvalue

# covariates
cnt_sex = logc.loc["ENSG00000129824"] if "ENSG00000129824" in logc.index else None  # RPS4Y1
meta["male"] = (cnt_sex > 3).astype(int)
meta["sebaceous"] = score(L["sebaceous_covariate"], [])
meta["batchB"] = (meta.batch == "B").astype(int)
seb_cut = meta.sebaceous.quantile(0.75)

COMPS = {"L_vs_HC": ("lesional", "control"), "PL_vs_HC": ("perilesional", "control"),
         "NL_vs_HC": ("nonlesional", "control")}

def boot_ci(s, case, ctrl):
    a = []
    for _ in range(NBOOT):  # resample individuals (each contributes one sample per comparison)
        c1 = rng.choice(case, len(case)); c0 = rng.choice(ctrl, len(ctrl))
        a.append(auc(s[c1].values, s[c0].values)[0])
    return np.percentile(a, [2.5, 97.5])

def evaluate(name, up, down, null_sets=None, n_up=None):
    s = score(up, down)
    rows = []
    for comp, (g1, g0) in COMPS.items():
        case = meta.index[meta.region == g1]; ctrl = meta.index[meta.region == g0]
        a, p = auc(s[case], s[ctrl])
        lo, hi = boot_ci(s, case, ctrl)
        r = dict(set=name, comparison=comp, n_case=len(case), n_ctrl=len(ctrl), AUC=a, CI_lo=lo, CI_hi=hi,
                 p_MW=p, AUC_agnostic=max(a, 1 - a))
        if null_sets is not None:
            na = np.array([auc(sc[case], sc[ctrl])[0] for sc in
                           (score(ns[:n_up], ns[n_up:]) for ns in null_sets)])
            r["null_AUC_median"] = np.median(na)
            r["p_null"] = (1 + (na >= a).sum()) / (1 + len(na))
        # covariate-adjusted: disease ~ z(score) + batch + sex + sebaceous
        idx = case.append(ctrl)
        X = pd.DataFrame({"z": (s[idx] - s[idx].mean()) / s[idx].std(), "batchB": meta.batchB[idx],
                          "male": meta.male[idx], "seb": (meta.sebaceous[idx] - meta.sebaceous[idx].mean()) / meta.sebaceous[idx].std()})
        y = idx.isin(case).astype(int)
        try:
            m = sm.Logit(y, sm.add_constant(X)).fit(disp=0, maxiter=200)
            r["adj_OR_per_SD"] = float(np.exp(m.params["z"])); r["adj_p"] = float(m.pvalues["z"])
        except Exception:
            r["adj_OR_per_SD"] = np.nan; r["adj_p"] = np.nan  # e.g. perfect separation
        # sensitivity: batch B only; excluding high-sebaceous samples
        cb, kb = case[meta.batch[case] == "B"], ctrl[meta.batch[ctrl] == "B"]
        r["AUC_batchB"] = auc(s[cb], s[kb])[0]
        ca, ka = case[meta.batch[case] == "A"], ctrl[meta.batch[ctrl] == "A"]
        r["AUC_batchA"] = auc(s[ca], s[ka])[0]
        cs, ks = case[meta.sebaceous[case] <= seb_cut], ctrl[meta.sebaceous[ctrl] <= seb_cut]
        r["AUC_lowSeb"] = auc(s[cs], s[ks])[0]; r["n_lowSeb"] = f"{len(cs)}v{len(ks)}"
        rows.append(r)
    # paired lesional vs non-lesional within the same patients
    l = meta[meta.region == "lesional"].set_index("pid"); n = meta[meta.region == "nonlesional"].set_index("pid")
    both = l.index.intersection(n.index)
    d = s[l.loc[both, "sample"] if "sample" in l else l.loc[both].index].values if False else None
    ls = pd.Series(s[meta.index[meta.region == "lesional"]].values, index=meta.pid[meta.region == "lesional"].values)
    ns_ = pd.Series(s[meta.index[meta.region == "nonlesional"]].values, index=meta.pid[meta.region == "nonlesional"].values)
    both = ls.index.intersection(ns_.index)
    diff = ls[both] - ns_[both]
    rows.append(dict(set=name, comparison="L_vs_NL_paired", n_case=len(both), n_ctrl=len(both),
                     AUC=float((diff > 0).mean() + 0.5 * (diff == 0).mean()), p_MW=stats.wilcoxon(diff).pvalue,
                     AUC_agnostic=max((diff > 0).mean(), (diff < 0).mean())))
    return rows, s

res, scores = [], {}
for pid, p in L["panels"].items():
    ns = L["null_sets"][pid]
    rows, s = evaluate(pid, p["up"], p["down"], ns["sets"], ns["n_up"])
    for r in rows: r.update(study=p["study"], tier=p["tier"], n_genes=len(p["up"]) + len(p["down"]))
    res += rows; scores[pid] = s
    # pigment-removal sensitivity
    for lab, key in [("noMelanogenesis", "melanogenesis_GO0042438"), ("noIdentity", "melanocyte_identity_GO0030318_HPA"), ("noPigment", "union")]:
        rm = set(L["pigment"][key])
        up2 = [g for g in p["up"] if g not in rm]; dn2 = [g for g in p["down"] if g not in rm]
        if len(up2) + len(dn2) == len(p["up"]) + len(p["down"]):
            continue  # no pigment genes in this panel
        if len(up2) + len(dn2) < 2:
            res.append(dict(set=f"{pid}_{lab}", comparison="L_vs_HC", study=p["study"], tier=p["tier"],
                            n_genes=len(up2) + len(dn2), note="fewer than 2 genes left; not scored"))
            continue
        rows, _ = evaluate(f"{pid}_{lab}", up2, dn2)
        for r in rows: r.update(study=p["study"], tier=p["tier"], n_genes=len(up2) + len(dn2))
        res += rows

# reference sets
for name, up, down in [("REF_IFNG_hallmark", L["ifng_hallmark"]["genes"], []),
                       ("REF_pigment_union_down", [], L["pigment"]["union"])]:
    rows, s = evaluate(name, up, down)
    for r in rows: r.update(study="reference", tier="reference", n_genes=len(up) + len(down))
    res += rows; scores[name] = s

df = pd.DataFrame(res)
# BH across the 14 published panels within each comparison (primary, as-reported scoring)
base = df.set.isin(L["panels"].keys())
for comp in df.comparison.dropna().unique():
    m = base & (df.comparison == comp)
    df.loc[m, "q_MW"] = multipletests(df.loc[m, "p_MW"], method="fdr_bh")[1]
    if "p_null" in df and df.loc[m, "p_null"].notna().all():
        df.loc[m, "q_null"] = multipletests(df.loc[m, "p_null"], method="fdr_bh")[1]
df.to_csv("results_all.tsv", sep="\t", index=False)
pd.DataFrame(scores).join(meta[["patient", "pid", "region", "batch", "male", "sebaceous"]]).to_csv("scores_per_sample.tsv", sep="\t")

pd.set_option("display.width", 260); pd.set_option("display.max_rows", 300)
cols = ["set", "study", "n_genes", "AUC", "CI_lo", "CI_hi", "AUC_agnostic", "null_AUC_median", "p_null", "q_MW", "adj_OR_per_SD", "adj_p", "AUC_batchA", "AUC_batchB", "AUC_lowSeb"]
for comp in ["L_vs_HC", "PL_vs_HC", "NL_vs_HC"]:
    print(f"\n=== {comp} ===")
    print(df[(df.comparison == comp) & ~df.set.str.contains("_no")][cols].round(3).to_string(index=False))
print("\n=== L_vs_NL_paired (fraction of patients with lesional score > non-lesional) ===")
print(df[(df.comparison == "L_vs_NL_paired") & ~df.set.str.contains("_no")][["set", "study", "n_case", "AUC", "AUC_agnostic", "p_MW", "q_MW"]].round(3).to_string(index=False))
print("\n=== pigment removal (L_vs_HC and PL_vs_HC) ===")
print(df[df.set.str.contains("_no") & df.comparison.isin(["L_vs_HC", "PL_vs_HC"])][["set", "comparison", "n_genes", "AUC", "CI_lo", "CI_hi", "note"]].round(3).to_string(index=False))
print("\nsex: controls male", int(meta.male[meta.region == 'control'].sum()), "/ 27; sebaceous 75th pct cut", round(seb_cut, 3))
