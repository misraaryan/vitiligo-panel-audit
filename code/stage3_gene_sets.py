"""Stage 3: define and lock all gene sets and controls BEFORE any panel is scored."""
import json, hashlib, numpy as np, pandas as pd

rng = np.random.default_rng(20261120)
N_NULL = 1000

ann = pd.read_csv("data/ncbi_gene_ensembl_map.tsv", sep="\t", dtype=str)
ann = ann.drop_duplicates("ensembl")
sym2ens = dict(zip(ann.symbol, ann.ensembl))
logc = pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0)

# universe: detected, protein-coding, autosomal
a = ann.set_index("ensembl")
pc = a.index[(a.type == "protein-coding") & (~a.chromosome.isin(["X", "Y", "MT", "X|Y"]))]
universe = logc.index.intersection(pc)
print("detected genes:", len(logc), "| protein-coding autosomal universe:", len(universe))

def ens(symbols):
    out, miss = [], []
    for s in symbols:
        e = sym2ens.get(s)
        (out if e else miss).append(e or s)
    return out, miss

# published panels with final directions
g = pd.read_csv("panels_gene_level.tsv", sep="\t")
panels = {}
for pid, d in g.groupby("panel", sort=False):
    panels[pid] = {"study": d.study.iloc[0], "tier": d.tier.iloc[0],
                   "up": [e for e, x in zip(d.ensembl, d.direction_final) if x == "up"],
                   "down": [e for e, x in zip(d.ensembl, d.direction_final) if x == "down"],
                   "symbols": dict(zip(d.ensembl, d.gene))}

# IFN-gamma comparator (MSigDB Hallmark, v2026.2)
j = json.load(open("data/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2026.2.Hs.json"))
ifng_syms = list(j.values())[0]["geneSymbols"]
ifng, ifng_miss = ens(ifng_syms)

# pigment lists (GO via QuickGO, human, incl. descendants; HPA single-cell 'melanocytes, cell type enriched')
GO_MELANIN = ["APPL1","ASIP","ATP7A","CDH3","CITED1","CTNS","DCT","GIPC1","MC1R","MFSD12","MGRN1","OCA2","OPN3",
              "PMEL","RAB38","RAPGEF2","SLC45A2","TRPC1","TYR","TYRP1","WNT5A","ZEB2"]          # GO:0042438
GO_MELANOCYTE_DIFF = ["CITED1","EDNRB","ENPP1","KIT","LRMDA","MEF2C","MITF","MREG","OCA2","OR51E2","TYRP1","USP13","ZEB2"]  # GO:0030318
HPA_MELANOCYTE = ["ACCSL","ALX1","CEP72","DCT","LCE2A","MLANA","TYR"]
melanogenesis, _ = ens(GO_MELANIN)
identity, _ = ens(sorted(set(GO_MELANOCYTE_DIFF + HPA_MELANOCYTE)))
pigment_union = sorted(set(melanogenesis) | set(identity))

# covariate signatures
SEBACEOUS = ["ELOVL3","AWAT1","AWAT2","DGAT2L6","FADS2","THRSP","MGST1","KRT79","SOAT1","PM20D1"]
sebaceous, _ = ens(SEBACEOUS)

# expression-matched null sets: each panel gene replaced by a random universe gene from the same
# (mean expression decile x SD quintile) bin
mu, sd = logc.loc[universe].mean(axis=1), logc.loc[universe].std(axis=1)
bins = pd.qcut(mu, 10, labels=False).astype(str) + "_" + pd.qcut(sd, 5, labels=False).astype(str)
by_bin = {b: list(ix) for b, ix in bins.groupby(bins).groups.items()}
def gene_bin(e):
    if e in bins.index: return bins[e]
    if e in logc.index:  # detected but not in universe (e.g. sex-linked): bin by its own stats
        m, s_ = logc.loc[e].mean(), logc.loc[e].std()
        return f"{int(np.clip(np.searchsorted(np.quantile(mu, np.linspace(0,1,11)[1:-1]), m),0,9))}_" \
               f"{int(np.clip(np.searchsorted(np.quantile(sd, np.linspace(0,1,6)[1:-1]), s_),0,4))}"
    return None
nulls = {}
for pid, p in panels.items():
    genes = [x for x in p["up"] + p["down"] if x in logc.index]
    gb = [gene_bin(x) for x in genes]
    sets = []
    for _ in range(N_NULL):
        sets.append([str(rng.choice([c for c in by_bin[b] if c not in genes])) for b in gb])
    nulls[pid] = {"n_up": len([x for x in p["up"] if x in logc.index]), "sets": sets}

locked = {"created": "2026-11-20", "seed": 20261120, "n_null": N_NULL,
          "universe_size": len(universe), "panels": panels,
          "ifng_hallmark": {"genes": ifng, "missing_symbols": ifng_miss,
                            "n_in_data": int(sum(e in logc.index for e in ifng))},
          "pigment": {"melanogenesis_GO0042438": melanogenesis, "melanocyte_identity_GO0030318_HPA": identity,
                      "union": pigment_union},
          "sebaceous_covariate": sebaceous, "null_sets": nulls,
          "rules": {"direction": "as reported; NR inferred from discovery data (GSE65127); sensitivity = direction-agnostic",
                    "primary_comparison": ["lesional vs healthy control", "perilesional vs healthy control"],
                    "secondary": ["non-lesional vs healthy control", "lesional vs non-lesional (paired)"],
                    "covariates": ["batch", "sex", "sebaceous score"],
                    "multiple_testing": "Benjamini-Hochberg across 14 panels within each comparison"}}
s = json.dumps(locked, indent=1)
open("gene_sets_locked.json", "w").write(s)
h = hashlib.sha256(s.encode()).hexdigest()
open("gene_sets_locked.sha256", "w").write(h + "  gene_sets_locked.json\n")
print("IFN-gamma hallmark:", len(ifng_syms), "symbols,", locked["ifng_hallmark"]["n_in_data"], "detected; missing:", ifng_miss)
print("melanogenesis:", len(melanogenesis), "| identity:", len(identity), "| union:", len(pigment_union))
print("sebaceous detected:", sum(e in logc.index for e in sebaceous), "/", len(sebaceous))
print("panel genes overlapping pigment union:",
      {pid: [p["symbols"][e] for e in p["up"] + p["down"] if e in pigment_union] for pid, p in panels.items() if any(e in pigment_union for e in p["up"] + p["down"])})
print("sha256:", h)
