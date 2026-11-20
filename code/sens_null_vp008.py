"""Post-hoc sensitivity (added 20 Nov 2026 after external critique; NOT part of the locked plan).
1) Null-baseline check: re-draw matched random gene sets with (a) the locked binning reproduced,
   (b) SD bins computed in healthy controls only, (c) mean+SD bins from controls only,
   (d) unrestricted random genes (expression decile only), (e) fully unrestricted random genes.
2) VP-008 check: drop the VP-008 non-X sample set and recompute perilesional/non-lesional AUCs."""
import json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

L = json.load(open("gene_sets_locked.json"))
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
R = logc.rank(axis=0) / len(logc)
Rv = R.values; gidx = {g: i for i, g in enumerate(R.index)}
N = 1000

# universe = genes used for the locked nulls
uni = sorted({g for p in L["null_sets"].values() for s in p["sets"] for g in s} )
# rebuild the true universe the same way the lock did: use the size recorded to sanity check
print("universe recorded:", L["universe_size"])
gm = pd.read_csv("data/ncbi_gene_ensembl_map.tsv", sep="\t", dtype=str).drop_duplicates("ensembl").set_index("ensembl")
pc = gm.index[(gm.type == "protein-coding") & (~gm.chromosome.isin(["X", "Y", "MT", "X|Y"]))]
universe = list(logc.index.intersection(pc))  # identical rule to stage3_gene_sets.py
print("universe rebuilt:", len(universe)); assert len(universe)==L["universe_size"]

ctrl = meta.index[meta.region == "control"]
def make_bins(mu, sd, nsd=5):
    return pd.qcut(mu, 10, labels=False).astype(str) + "_" + (pd.qcut(sd, nsd, labels=False).astype(str) if nsd else "0")
U = logc.loc[universe]
variants = {
    "locked_all_samples": make_bins(U.mean(1), U.std(1)),
    "SD_in_controls": make_bins(U.mean(1), U[ctrl].std(1)),
    "mean_SD_in_controls": make_bins(U[ctrl].mean(1), U[ctrl].std(1)),
    "expression_only": make_bins(U.mean(1), None, 0),
    "unrestricted": pd.Series("0", index=universe),
}

def score_idx(up_i, dn_i):
    s = np.zeros(Rv.shape[1])
    if up_i: s += Rv[up_i].mean(0)
    if dn_i: s -= Rv[dn_i].mean(0)
    return s

def auc(s, case, ctl):
    x, y = s[case], s[ctl]
    return stats.mannwhitneyu(x, y).statistic / (len(x) * len(y))

col = {c: i for i, c in enumerate(R.columns)}
def idx(names): return [col[n] for n in names]
comps = {"L_vs_HC": ("lesional", "control"), "PL_vs_HC": ("perilesional", "control")}

out = []
for vname, bins in variants.items():
    by_bin = {b: np.array([gidx[g] for g in ix]) for b, ix in bins.groupby(bins).groups.items()}
    rng = np.random.default_rng(20261120)
    for comp, (a, b) in comps.items():
        case = idx(meta.index[meta.region == a]); ct = idx(meta.index[meta.region == b])
        for pid, p in L["panels"].items():
            up = [g for g in p["up"] if g in gidx]; dn = [g for g in p["down"] if g in gidx]
            genes = up + dn
            obs = auc(score_idx([gidx[g] for g in up], [gidx[g] for g in dn]), case, ct)
            def gbin(g):
                if g in bins.index: return bins[g]
                # sex-linked/non-universe gene: nearest bin by mean decile, middle SD
                return bins.iloc[(U.mean(1) - logc.loc[g].mean()).abs().argmin()]
            gb = [gbin(g) for g in genes]
            excl = {gidx[g] for g in genes}
            na = []
            for _ in range(N):
                pick = []
                for bb in gb:
                    pool = by_bin[bb]
                    while True:
                        c = int(rng.choice(pool))
                        if c not in excl: break
                    pick.append(c)
                na.append(auc(score_idx(pick[:len(up)], pick[len(up):]), case, ct))
            na = np.array(na)
            out.append(dict(variant=vname, comparison=comp, set=pid, AUC=obs, null_median=np.median(na),
                            p_null=(1 + (na >= obs).sum()) / (1 + N)))
df = pd.DataFrame(out)
for (v, c), d in df.groupby(["variant", "comparison"]):
    df.loc[d.index, "q_null"] = multipletests(d.p_null, method="fdr_bh")[1]
df["passes_FDR"] = df.q_null < 0.05
df.to_csv("sens_null_variants.tsv", sep="\t", index=False)
piv = df.pivot_table(index=["comparison", "set"], columns="variant", values=["null_median", "q_null"])
pd.set_option("display.width", 250)
print(piv.round(3).to_string())

# ---- VP-008 sensitivity
print("\nVP-008 non-X removal")
meta["x"] = meta.x.fillna(False).astype(bool)
drop = meta.index[(meta.patient == "VP008") & (~meta.x)]
print("dropping", list(drop))
rows = []
for comp, reg in [("PL_vs_HC", "perilesional"), ("NL_vs_HC", "nonlesional")]:
    for keep_all in [True, False]:
        m2 = meta if keep_all else meta.drop(drop)
        case = idx(m2.index[m2.region == reg]); ct = idx(m2.index[m2.region == "control"])
        for pid, p in L["panels"].items():
            up = [gidx[g] for g in p["up"] if g in gidx]; dn = [gidx[g] for g in p["down"] if g in gidx]
            rows.append(dict(comparison=comp, version="all" if keep_all else "noVP008nonX", n_case=len(case),
                             set=pid, AUC=auc(score_idx(up, dn), case, ct)))
v = pd.DataFrame(rows).pivot_table(index=["comparison", "set"], columns="version", values="AUC")
v["diff"] = v["noVP008nonX"] - v["all"]
print(v.round(3).to_string()); print("max |diff|:", v["diff"].abs().max().round(4))
v.to_csv("sens_vp008.tsv", sep="\t")
