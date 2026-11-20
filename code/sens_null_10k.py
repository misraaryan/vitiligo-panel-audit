"""Post-hoc: locked binning (mean/SD over all samples) re-drawn with 10,000 null sets per panel,
to remove Monte Carlo noise from p_null (locked analysis used 1,000). Also exploratory melanocyte analysis
without circularity (melanocyte score excludes each panel's own genes)."""
import json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

L = json.load(open("gene_sets_locked.json"))
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta.tsv", sep="\t", index_col=0)
R = logc.rank(axis=0) / len(logc); Rv = R.values; gidx = {g: i for i, g in enumerate(R.index)}
gm = pd.read_csv("data/ncbi_gene_ensembl_map.tsv", sep="\t", dtype=str).drop_duplicates("ensembl").set_index("ensembl")
pc = gm.index[(gm.type == "protein-coding") & (~gm.chromosome.isin(["X", "Y", "MT", "X|Y"]))]
universe = list(logc.index.intersection(pc))  # identical rule to stage3_gene_sets.py
assert len(universe) == L["universe_size"], (len(universe), L["universe_size"])
U = logc.loc[universe]
mu, sd = U.mean(axis=1), U.std(axis=1)
bins = pd.qcut(mu, 10, labels=False).astype(str) + "_" + pd.qcut(sd, 5, labels=False).astype(str)
qmu = np.quantile(mu, np.linspace(0, 1, 11)[1:-1]); qsd = np.quantile(sd, np.linspace(0, 1, 6)[1:-1])
def gbin(g):  # identical rule to stage3_gene_sets.py
    if g in bins.index: return bins[g]
    m, s_ = logc.loc[g].mean(), logc.loc[g].std()
    return f"{int(np.clip(np.searchsorted(qmu, m), 0, 9))}_{int(np.clip(np.searchsorted(qsd, s_), 0, 4))}"
by_bin = {b: np.array([gidx[g] for g in ix]) for b, ix in bins.groupby(bins).groups.items()}
col = {c: i for i, c in enumerate(R.columns)}
def idx(n): return np.array([col[x] for x in n])
def aucv(S, case, ct):  # S: (nsets, nsamples) -> AUC per row, via ranks
    x = S[:, np.r_[case, ct]]; r = stats.rankdata(x, axis=1)
    return (r[:, :len(case)].sum(1) - len(case) * (len(case) + 1) / 2) / (len(case) * len(ct))
NN = 10000; rng = np.random.default_rng(20261121)
rows = []
for comp, reg in [("L_vs_HC", "lesional"), ("PL_vs_HC", "perilesional"), ("NL_vs_HC", "nonlesional")]:
    case = idx(meta.index[meta.region == reg]); ct = idx(meta.index[meta.region == "control"])
    for pid, p in L["panels"].items():
        up = [gidx[g] for g in p["up"] if g in gidx]; dn = [gidx[g] for g in p["down"] if g in gidx]
        obs_s = (Rv[up].mean(0) if up else 0) - (Rv[dn].mean(0) if dn else 0)
        obs = aucv(np.atleast_2d(obs_s), case, ct)[0]
        genes = up + dn; gb = [gbin(R.index[i]) for i in genes]
        picks = np.stack([rng.choice(by_bin[b][~np.isin(by_bin[b], genes)], NN) for b in gb], 1)
        S = np.zeros((NN, Rv.shape[1]))
        if up: S += Rv[picks[:, :len(up)]].mean(1)
        if dn: S -= Rv[picks[:, len(up):]].mean(1)
        na = aucv(S, case, ct)
        rows.append(dict(comparison=comp, set=pid, AUC=round(obs, 4), null_median_10k=round(np.median(na), 3),
                         p_null_10k=(1 + (na >= obs).sum()) / (1 + NN)))
d = pd.DataFrame(rows)
for c, g in d.groupby("comparison"): d.loc[g.index, "q_null_10k"] = multipletests(g.p_null_10k, method="fdr_bh")[1]
d.to_csv("sens_null_10k.tsv", sep="\t", index=False)
print(d.round(4).to_string())

# ---- exploratory melanocyte analysis without circularity
pig = L["pigment"]["union"]
L_ = (meta.region == "lesional").values; C = (meta.region == "control").values; PL = (meta.region == "perilesional").values
def a(s, m1, m0): return stats.mannwhitneyu(s[m1], s[m0]).statistic / (m1.sum() * m0.sum())
ex = []
for pid, p in L["panels"].items():
    up = [g for g in p["up"] if g in gidx]; dn = [g for g in p["down"] if g in gidx]
    s = (R.loc[up].mean().values if up else 0) - (R.loc[dn].mean().values if dn else 0)
    own = set(up) | set(dn); pg = [g for g in pig if g in gidx and g not in own]
    mel = R.loc[pg].mean().values
    lin = s - np.polyval(np.polyfit(mel, s, 1), mel)
    rs, rm = stats.rankdata(s), stats.rankdata(mel); rr = rs - np.polyval(np.polyfit(rm, rs, 1), rm)
    ex.append(dict(set=pid, n_overlap_with_pigment_list=len(own & set(pig)),
                   rho_melanocyte=round(stats.spearmanr(s, mel)[0], 2), AUC_L=round(a(s, L_, C), 3),
                   AUC_L_resid_linear=round(a(lin, L_, C), 3), AUC_L_resid_rank=round(a(rr, L_, C), 3),
                   AUC_PL_resid_linear=round(a(lin, PL, C), 3)))
e = pd.DataFrame(ex); e.to_csv("results_exploratory_melanocyte_v2.tsv", sep="\t", index=False)
print(e.to_string())
