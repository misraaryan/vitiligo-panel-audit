"""Read the GSE298871 10x matrices (GSE298871_RAW.tar extracted into <raw_dir>) and save them as one pickle.
python code/sc_load_GSE298871.py <raw_dir> sc.pkl"""
import gzip, glob, re, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, pickle, sys
RAW = sys.argv[1]; OUT = sys.argv[2]
samples = {}
for f in sorted(glob.glob(f"{RAW}/*_matrix.mtx.gz")):
    pre = f[:-len("matrix.mtx.gz")]
    gsm, name = re.match(r".*/(GSM\d+)_(.+)_$", pre).groups()
    m = scipy.io.mmread(f).tocsc()                      # genes x cells
    feat = pd.read_csv(pre+"features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(pre+"barcodes.tsv.gz", sep="\t", header=None)[0].values
    samples[name] = dict(gsm=gsm, m=m, feat=feat, bc=bc)
f0 = samples[next(iter(samples))]["feat"]
assert all((s["feat"][0].values == f0[0].values).all() for s in samples.values())
print(f0.head(3).to_string()); print(f0[2].value_counts().to_dict() if 2 in f0 else "")
for k,s in samples.items():
    tot = np.asarray(s["m"].sum(0)).ravel(); ng = np.asarray((s["m"]>0).sum(0)).ravel()
    print(k, s["gsm"], "cells", s["m"].shape[1], "median UMI", int(np.median(tot)), "median genes", int(np.median(ng)), "min UMI", int(tot.min()))
pickle.dump(dict(samples=samples, feat=f0), open(OUT,"wb"))
