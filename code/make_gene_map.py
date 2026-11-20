"""Build data/ncbi_gene_ensembl_map.tsv (symbol, ensembl, chromosome, type, synonyms)
from NCBI Homo_sapiens.gene_info.gz. Usage: python code/make_gene_map.py data/Homo_sapiens.gene_info.gz"""
import sys
import pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "data/Homo_sapiens.gene_info.gz"
g = pd.read_csv(src, sep="\t", dtype=str, compression="gzip" if src.endswith(".gz") else None)
g["ensembl"] = g["dbXrefs"].str.extract(r"Ensembl:(ENSG\d+)", expand=False)
out = (g.dropna(subset=["ensembl"])
        .rename(columns={"Symbol": "symbol", "chromosome": "chromosome", "type_of_gene": "type", "Synonyms": "synonyms"})
        [["symbol", "ensembl", "chromosome", "type", "synonyms"]])
out.to_csv("data/ncbi_gene_ensembl_map.tsv", sep="\t", index=False)
print(len(out), "genes written to data/ncbi_gene_ensembl_map.tsv")
