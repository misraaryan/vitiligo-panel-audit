"""Machine-learning figure and summary table. (a) retrained LR CV AUC per panel vs matched random gene sets,
with the training-free AUC for contrast; (b) genes most often chosen by genome-wide LASSO, coloured by
melanocyte enrichment in the independent single-cell data.
Also writes ml_table_summary.tsv (training-free AUC and q, retrained LR/SVM/RF q vs random sets, and
the GSE311770 -> GSE298871 test). Run from the repository root:
python code/make_figure_ml.py [results_dir]   (default: results)"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans", "font.size": 7})
import matplotlib.pyplot as plt
import sys
V = (sys.argv[1] if len(sys.argv) > 1 else "results") + "/"
e1 = pd.read_csv(V + "ml_e1_panels_vs_random.tsv", sep="\t").set_index("panel")
ra = pd.read_csv(V + "results_all.tsv", sep="\t"); ra = ra[ra.comparison == "L_vs_HC"].set_index("set")
e1b = pd.read_csv(V + "ml_e1b_other_models.tsv", sep="\t")
e3 = pd.read_csv(V + "ml_e3_cross_cohort.tsv", sep="\t").set_index("panel")
t = pd.DataFrame({"tf_AUC": ra.AUC, "tf_q": ra.q_null}).loc[e1.index]
t["lr_cv"] = e1.cv_AUC; t["rand_med"] = e1.null_median; t["q_lr"] = e1.q_BH
for m, c in [("svm_rbf", "q_svm"), ("rf", "q_rf")]:
    t[c] = e1b[e1b.model == m].set_index("panel").q_BH
t["x_trained"] = e3.trained_LR_test_AUC; t["x_tf"] = e3.training_free_test_AUC; t["x_q"] = e3.q_BH
t.to_csv("ml_table_summary.tsv", sep="\t")
las = pd.read_csv(V + "ml_e2_lasso_selected_genes_annotated.tsv", sep="\t")
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.3), gridspec_kw={"width_ratios": [1.15, 1]})
order = list(e1.index)[::-1]; y = np.arange(len(order))
for i, p in enumerate(order):
    a.plot([e1.loc[p, "null_median"], e1.loc[p, "null_q95"]], [i, i], color="#B4B2A9", lw=4, solid_capstyle="butt")
    a.plot(e1.loc[p, "cv_AUC"], i, "o", color="#185FA5", ms=4.5)
    a.plot(t.loc[p, "tf_AUC"], i, "o", mfc="none", mec="#A32D2D", ms=4.5, mew=1)
a.axvline(0.5, ls="--", lw=0.7, color="#888780")
a.set_yticks(y); a.set_yticklabels(order); a.set_xlim(0, 1); a.set_xlabel("AUC, lesional vs healthy (GSE311770)")
a.set_title("(a) Retrained panels vs matched random gene sets", fontsize=7.5, loc="left")
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], marker="o", color="#185FA5", ls="", ms=4.5, label="Retrained LR (5-fold CV)"),
                  Line2D([], [], marker="o", mfc="none", mec="#A32D2D", ls="", ms=4.5, label="Training-free score"),
                  Line2D([], [], color="#B4B2A9", lw=4, label="Random sets, median to 95th pct")],
         loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=6, frameon=False, handlelength=1.4, columnspacing=1.0)
top = las[las.selection_freq >= 0.44].copy().sort_values("selection_freq")
top["symbol"] = top.symbol.fillna(top.ensembl)
cat = np.where(top.is_typing_marker == True, "Melanocyte marker", np.where(top.melanocyte_enriched == True, "Melanocyte-enriched (scRNA)", "Other"))
col = {"Melanocyte marker": "#0F6E56", "Melanocyte-enriched (scRNA)": "#5DCAA5", "Other": "#B4B2A9"}
b.barh(range(len(top)), top.selection_freq, color=[col[c] for c in cat], height=0.75)
b.set_yticks(range(len(top))); b.set_yticklabels(top.symbol, fontsize=6); b.set_xlim(0, 1.0)
b.set_xlabel("Fraction of 25 nested-CV models selecting the gene")
b.set_title("(b) Genes chosen by genome-wide LASSO", fontsize=7.5, loc="left")
from matplotlib.patches import Patch
b.legend(handles=[Patch(color=v, label=k) for k, v in col.items()], loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=6, frameon=False, handlelength=1.4, columnspacing=1.0)
for ax in (a, b):
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(w_pad=2.0)
fig.savefig("fig_ml.pdf", bbox_inches="tight"); fig.savefig("fig_ml.png", dpi=200, bbox_inches="tight")
print(len(top), list(zip(top.symbol, cat)))
