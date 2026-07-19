#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg
from statannotations.Annotator import Annotator

df = pd.read_csv("input_expression.csv"); df["genotype"] = df["genotype"].astype(str)
order = ['het', 'ko', 'rescue', 'wt']
row = pg.anova(data=df, dv="expression", between="genotype").iloc[0]; omni_p = float(row.p_unc)
ph = pg.pairwise_tukey(data=df, dv="expression", between="genotype"); pairs = [(str(r.A), str(r.B), float(r.p_tukey)) for r in ph.itertuples()]
stars = lambda p: "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*"
sig = [(a, b, p) for a, b, p in pairs if p < 0.05]

fig, ax = plt.subplots(figsize=(4.6, 4.2))
sns.boxplot(data=df, x="genotype", y="expression", hue="genotype", order=order, palette=['#3B6DB3', '#C1432B', '#2E8B57', '#7A5195'], ax=ax, legend=False, width=0.6, fliersize=0)
sns.stripplot(data=df, x="genotype", y="expression", hue="genotype", order=order, palette=['#3B6DB3', '#C1432B', '#2E8B57', '#7A5195'], ax=ax,
              legend=False, size=3.5, alpha=0.7, edgecolor="black", linewidth=0.3, jitter=0.12)
if sig:
    annot = Annotator(ax, [(a, b) for a, b, _ in sig], data=df, x="genotype", y="expression", order=order)
    annot.configure(line_width=0.9, fontsize=10)
    annot.set_custom_annotations([stars(p) for _, _, p in sig])
    annot.annotate()
ax.set_xlabel("genotype"); ax.set_ylabel("Expression (a.u.)")
ax.set_title("One-way ANOVA, P = %s" % format(omni_p, ".2g"), fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_expression_bar_py.pdf")
print("Reproduced figure.")
