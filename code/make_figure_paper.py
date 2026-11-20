import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("results_all.tsv", sep="\t")
order = ["P01","P02","P03","P04","P05","P06","P07","P08","P09","S01","S02","S03","S04","S05","REF_IFNG_hallmark","REF_pigment_union_down"]
lab = {r.set: (f"{r.set} {r.study}" if not r.set.startswith("REF") else ("IFN-γ hallmark (ref)" if "IFNG" in r.set else "Pigment genes (ref)"))
       for r in d.itertuples() if r.set in order}
# mark panels scored on fewer genes than published (genes below the expression filter)
import json
L = json.load(open("gene_sets_locked.json")); detected = set(pd.read_csv("data/logcpm.tsv.gz", sep="\t", index_col=0, usecols=[0]).index)
for pid, p in L["panels"].items():
    g = p["up"] + p["down"]; k = sum(x in detected for x in g)
    if k < len(g): lab[pid] += f" ({k}/{len(g)} genes)"
fig, axs = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
for ax, comp, title in [(axs[0], "L_vs_HC", "Lesional vs healthy (24 vs 27)"), (axs[1], "PL_vs_HC", "Perilesional vs healthy (27 vs 27)")]:
    x = d[(d.comparison == comp) & d.set.isin(order)].set_index("set").loc[order]; y = np.arange(len(order))[::-1]
    for yi, (s, r) in zip(y, x.iterrows()):
        v = str(r.verdict)
        col = "#5F5E5A" if s.startswith("REF") else "#A32D2D" if v.startswith("REVERSED") else "#0F6E56" if v.startswith("REPLICATES") \
            else "#5DCAA5" if v.startswith("beats random") else "#B4B2A9"
        ax.plot([r.CI_lo, r.CI_hi], [yi, yi], color=col, lw=2); ax.plot(r.AUC, yi, "o", color=col, ms=6)
        if pd.notna(r.get("null_AUC_median")): ax.plot(r.null_AUC_median, yi, "|", color="#185FA5", ms=10, mew=1.5)
    ax.axvline(0.5, ls="--", color="#888780", lw=0.8); ax.set_xlim(0, 1); ax.set_title(title, fontsize=11)
    ax.set_xlabel("AUC (95% bootstrap CI)"); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_yticks(np.arange(len(order))[::-1]); axs[0].set_yticklabels([lab[o] for o in order], fontsize=9)
#fig.text(0.01, 0.035, "Dark green: CI > 0.5 and beats 1,000 expression-matched random gene sets after FDR (q < 0.05).  Light green: beats random sets before FDR only.", fontsize=8)
#fig.text(0.01, 0.01, "Red: CI < 0.5 (reported direction contradicted).  Grey: fails.  Blue tick: median AUC of the matched random gene sets.", fontsize=8)
plt.tight_layout(w_pad=2.5); plt.savefig("fig_forest_AUC_paper.pdf", bbox_inches="tight", pad_inches=0.08); plt.savefig("fig_forest_AUC_paper.png", dpi=200, bbox_inches="tight", pad_inches=0.08)
