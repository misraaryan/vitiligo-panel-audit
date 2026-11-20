import re, gzip, numpy as np, pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

cnt = pd.read_csv("data/GSE311770_all_samples_rsem_gene_counts.tsv.gz", sep="\t", index_col=0)
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True)
cnt = cnt.groupby(level=0).sum()

def parse(col):
    c = col.lower()
    m = re.search(r"contro(?:l|le)_(\d+)", c)
    if m:
        n = int(m.group(1))
        return dict(sample=col, patient=f"C{n:02d}", region="control", batch="A" if n <= 9 else "B")
    m = re.search(r"vp_(\d{3})", c)
    pid = m.group(1)
    region = ("perilesional" if "perilesioneel" in c else
              "lesional" if "lesioneel" in c else
              "nonlesional" if "gezond" in c else "?")
    batch = "A" if pid in {"001","002","003","004","005","006","007","009","010","011","012","014","016","018","019"} else "B"
    return dict(sample=col, patient=f"VP{pid}", region=region, batch=batch, x="_x_" in c)

meta = pd.DataFrame([parse(c) for c in cnt.columns]).set_index("sample")
meta.to_csv("data/meta.tsv", sep="\t")
print("genes x samples:", cnt.shape)
print(pd.crosstab(meta.region, meta.batch, margins=True))
print("patients:", meta[meta.region != "control"].patient.nunique(),
      "| VP008 samples:", list(meta[meta.patient == "VP008"].index))

lib = cnt.sum()
print("\nlibrary size (M reads): min %.2f median %.2f max %.2f" % (lib.min()/1e6, lib.median()/1e6, lib.max()/1e6))
print(lib.groupby(meta.batch).median().div(1e6).round(2).rename("median lib M by batch").to_string())
print(lib.groupby(meta.region).median().div(1e6).round(2).rename("median lib M by region").to_string())
low = lib[lib < 0.25 * lib.median()]
print("very low-depth samples (<25% median):", list(low.index))

cpm = cnt / lib * 1e6
keep = (cpm > 1).sum(axis=1) >= 20
logc = np.log2(cpm[keep] + 1)
print("\ngenes kept (CPM>1 in >=20 samples):", keep.sum())

X = (logc.sub(logc.mean(1), axis=0)).T
top = logc.var(1).sort_values(ascending=False).index[:2000]
pca = PCA(5).fit(X[top])
pcs = pd.DataFrame(pca.transform(X[top]), index=X.index, columns=[f"PC{i+1}" for i in range(5)])
print("variance explained:", (pca.explained_variance_ratio_ * 100).round(1))
for pc in pcs.columns[:4]:
    pb = stats.mannwhitneyu(pcs.loc[meta.batch == "A", pc], pcs.loc[meta.batch == "B", pc]).pvalue
    pr = stats.kruskal(*[pcs.loc[meta.region == r, pc] for r in meta.region.unique()]).pvalue
    print(f"{pc}: batch p={pb:.1e}  region p={pr:.1e}")

markers = {"TYR": "ENSG00000077498", "DCT": "ENSG00000080166", "PMEL": "ENSG00000185664",
           "MLANA": "ENSG00000120215", "CXCL9": "ENSG00000138755", "CXCL10": "ENSG00000169245"}
print("\nmedian log2CPM by region (sanity check of labels):")
tab = pd.DataFrame({g: logc.loc[e].groupby(meta.region).median() if e in logc.index else np.nan
                    for g, e in markers.items()}).round(2)
print(tab.to_string())

ctrl = meta.region == "control"
print("\nbatch effect within controls only (same biology):")
for pc in ["PC1", "PC2"]:
    print(pc, "controls A vs B p=%.1e" % stats.mannwhitneyu(pcs.loc[ctrl & (meta.batch == "A"), pc],
                                                           pcs.loc[ctrl & (meta.batch == "B"), pc]).pvalue)
pcs.join(meta).to_csv("data/pca.tsv", sep="\t")
logc.to_csv("data/logcpm.tsv.gz", sep="\t", compression="gzip")
