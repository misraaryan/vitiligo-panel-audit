"""Post hoc: check the 13 directions inferred from GSE65127 against GSE75819
(Illumina GPL6884; 15 patients, paired lesional (V) vs non-lesional (N) skin).
Probe per gene = highest mean signal. Values are BeadStudio average-normalised
(can be negative), so expression = log2(max(x,0)+16). Paired L-N difference per
patient, Wilcoxon signed-rank test. CACTIN (C19orf29) is not on GPL6884."""
import pandas as pd, numpy as np
from scipy.stats import wilcoxon
d = pd.read_csv("verification/GSE75819_selected_probes.tsv", sep="\t")
d["pair"] = d.title.str[:-1]; d["grp"] = d.title.str[-1]
probes = {"GPR143":["ILMN_1756261"],"DCTN1":["ILMN_1789089","ILMN_2412807","ILMN_1793563"],
 "KLRG1":["ILMN_2124920","ILMN_1658399"],"MAPK14":["ILMN_2388090","ILMN_1720656","ILMN_1788002","ILMN_1737627"],
 "MRPL47":["ILMN_2391522","ILMN_1687036"],"EXOSC7":["ILMN_2055477"],"MLANA":["ILMN_1710217"],
 "NUF2":["ILMN_1658695","ILMN_2406724"],"TYR":["ILMN_1788774"],"NKTR":["ILMN_2330994","ILMN_1721621"],
 "TYRP1":["ILMN_2054652"],"LARP7":["ILMN_2392717"],"DCT":["ILMN_1701783"],"PMEL":["ILMN_1665994"],"HADH":["ILMN_1719906"]}
inferred = [("P05","DCTN1","down"),("P05","GPR143","down"),("P05","HADH","down"),("P05","MRPL47","up"),
 ("P05","NKTR","up"),("P05","NUF2","up"),("P05","CACTIN","up"),("P06","TYR","down"),("P06","TYRP1","down"),
 ("P06","DCT","down"),("P06","LARP7","up"),("S02","EXOSC7","up"),("S02","KLRG1","up")]
controls = [("ctrl","MLANA","down"),("ctrl","PMEL","down"),("ctrl","MAPK14","up")]
L = d[d.grp=="V"].set_index("pair"); N = d[d.grp=="N"].set_index("pair")
pairs = sorted(set(L.index)&set(N.index)); print("pairs:", len(pairs), "groups:", d.grp.value_counts().to_dict())
rows=[]
for panel,g,inf in inferred+controls:
    if g not in probes:
        rows.append(dict(panel=panel,gene=g,gse65127_dir=inf,probe="not on GPL6884")); continue
    pr = max(probes[g], key=lambda p: d[p].mean())
    f = lambda x: np.log2(np.clip(x,0,None)+16)
    diff = f(L.loc[pairs,pr].values) - f(N.loc[pairs,pr].values)
    p = wilcoxon(diff).pvalue
    gdir = "up" if np.median(diff)>0 else "down"
    rows.append(dict(panel=panel,gene=g,gse65127_dir=inf,probe=pr,n_pairs=len(pairs),
        median_log2_LminusN=round(float(np.median(diff)),3),frac_L_higher=round(float((diff>0).mean()),3),
        wilcoxon_p=float(f"{p:.2g}"),gse75819_dir=gdir,agree=(gdir==inf),
        agree_and_p05=(gdir==inf and p<0.05)))
o = pd.DataFrame(rows); o.to_csv("posthoc_direction_GSE75819.tsv",sep="\t",index=False)
print(o.to_string(index=False))
t = o[(o.panel!="ctrl")&o.probe.ne("not on GPL6884")]
print(f"\nTestable inferred genes: {len(t)}; agree: {t.agree.sum()}; agree with p<0.05: {t.agree_and_p05.sum()}; disagree with p<0.05: {((~t.agree)&(t.wilcoxon_p<0.05)).sum()}")
