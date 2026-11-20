import gzip, io, numpy as np, pandas as pd
from scipy import stats

# GSE65127 series matrix
lines = gzip.open("data/GSE65127_series_matrix.txt.gz", "rt").read().splitlines()
titles = next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]
titles = [t.strip('"') for t in titles]
s = lines.index("!series_matrix_table_begin") + 1
e = lines.index("!series_matrix_table_end")
mat = pd.read_csv(io.StringIO("\n".join(lines[s:e])), sep="\t", index_col=0)
mat.columns = titles
grp = pd.Series(["control" if t.startswith("Healthy") else
                 "lesional" if t.startswith("Lesional") else
                 "perilesional" if t.startswith("Peri") else "nonlesional" for t in titles], index=titles)
print(grp.value_counts().to_dict(), "| value range", float(mat.min().min()), float(mat.max().max()))

ann = pd.read_csv("data/GPL570-55999.txt", sep="\t", comment="#", usecols=["ID", "Gene Symbol"], dtype=str)
ann = ann.dropna()
ann = ann[~ann["Gene Symbol"].str.contains("///")]
mat = mat.loc[mat.index.intersection(ann.ID)]
mat["sym"] = ann.set_index("ID").loc[mat.index, "Gene Symbol"]
# collapse probes: keep the probe with highest mean expression per gene
mat["mean"] = mat[titles].mean(axis=1)
gene = mat.sort_values("mean", ascending=False).drop_duplicates("sym").set_index("sym")[titles]
gene.to_csv("data/GSE65127_gene.tsv.gz", sep="\t", compression="gzip")

L = grp[grp == "lesional"].index; C = grp[grp == "control"].index
def direction(sym):
    if sym not in gene.index:
        return "absent", np.nan, np.nan
    d = gene.loc[sym, L].mean() - gene.loc[sym, C].mean()
    p = stats.mannwhitneyu(gene.loc[sym, L], gene.loc[sym, C]).pvalue
    return ("up" if d > 0 else "down"), round(d, 2), p

g = pd.read_csv("panels_gene_level.tsv", sep="\t")
out = []
for _, r in g.iterrows():
    dd, lfc, p = direction(r.gene)
    out.append(dict(panel=r.panel, gene=r.gene, reported=r.direction, gse65127_dir=dd, log2diff=lfc, p=p))
o = pd.DataFrame(out)
o["p"] = o.p.map(lambda x: f"{x:.1e}" if pd.notna(x) else "")
print("\nNR genes, inferred from GSE65127 (lesional vs healthy):")
print(o[o.reported == "NR"].to_string(index=False))
print("\nReported directions that CONTRADICT GSE65127:")
c = o[(o.reported != "NR") & (o.gse65127_dir.isin(["up", "down"])) & (o.reported != o.gse65127_dir)]
print(c.to_string(index=False))
o.to_csv("direction_check_GSE65127.tsv", sep="\t", index=False)

# apply rule: fill NR with GSE65127 direction
g["direction_final"] = g.direction
m = g.direction == "NR"
g.loc[m, "direction_final"] = [direction(x)[0] for x in g.loc[m, "gene"]]
g.loc[m, "direction_source"] = "inferred_GSE65127"
g.to_csv("panels_gene_level.tsv", sep="\t", index=False)
