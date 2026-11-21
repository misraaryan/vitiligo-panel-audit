"""Post hoc: how the replication count depends on the decision rule.
Locked rule = bootstrap CI lower bound > 0.5 AND BH q (matched null, 1,000 sets) < 0.05."""
import pandas as pd
r = pd.read_csv("results_all.tsv", sep="\t"); k = pd.read_csv("sens_null_10k.tsv", sep="\t")
r = r[r.set.str.match(r"^[PS]\d\d$") & r.comparison.isin(["L_vs_HC","PL_vs_HC"])]
r = r.merge(k[["comparison","set","q_null_10k"]], on=["comparison","set"])
rules = {"CI lower bound > 0.5 only": lambda d: d.CI_lo > 0.5,
         "Locked: CI > 0.5 and q < 0.05 (1,000 nulls)": lambda d: (d.CI_lo > 0.5) & (d.q_null < 0.05),
         "CI > 0.5 and q < 0.05 (10,000 nulls)": lambda d: (d.CI_lo > 0.5) & (d.q_null_10k < 0.05),
         "CI > 0.5 and q < 0.01 (10,000 nulls)": lambda d: (d.CI_lo > 0.5) & (d.q_null_10k < 0.01)}
rows = []
for name, f in rules.items():
    row = {"rule": name}
    for comp, lab in [("L_vs_HC","lesional_vs_healthy"),("PL_vs_HC","perilesional_vs_healthy")]:
        d = r[r.comparison == comp]; s = sorted(d.set[f(d)])
        row[lab+"_n"] = len(s); row[lab+"_panels"] = ", ".join(s) if s else "none"
    rows.append(row)
o = pd.DataFrame(rows); o.to_csv("posthoc_rule_sensitivity.tsv", sep="\t", index=False); print(o.to_string(index=False))
