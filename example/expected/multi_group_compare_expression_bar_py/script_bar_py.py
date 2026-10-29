#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*"


def bracket_key(row):
    i1, i2 = order.index(row[0]), order.index(row[1])
    return (abs(i2 - i1), min(i1, i2))


df = pd.read_csv("input_expression.csv"); df["genotype"] = df["genotype"].astype(str)
order = ['het', 'ko', 'rescue', 'wt']
row = pg.anova(data=df, dv="expression", between="genotype").iloc[0]; omni_p = float(row.p_unc)
ph = pg.pairwise_tukey(data=df, dv="expression", between="genotype"); pairs = [(str(r.A), str(r.B), float(r.p_tukey)) for r in ph.itertuples()]
sig = [(a, b, p) for a, b, p in pairs if p < 0.05]

fig, ax = plt.subplots(figsize=(5.4, 4.2))
sns.boxplot(data=df, x="genotype", y="expression", hue="genotype", order=order, hue_order=order, palette=['#3B6DB3', '#C1432B', '#2E8B57', '#7A5195'], ax=ax, legend=False, width=0.6, fliersize=0)
sns.stripplot(data=df, x="genotype", y="expression", hue="genotype", order=order, hue_order=order, palette=['#3B6DB3', '#C1432B', '#2E8B57', '#7A5195'], ax=ax,
              legend=False, size=3.5, alpha=0.7, edgecolor="black", linewidth=0.3, jitter=0.12)
span = df["expression"].max() - df["expression"].min(); top = df["expression"].max()
pos = {g: i for i, g in enumerate(order)}
sig.sort(key=bracket_key)
for k, (a, b, p) in enumerate(sig):
    yk = top + span * (0.08 + 0.13 * k); x1, x2 = sorted((pos[a], pos[b]))
    ax.plot([x1, x1, x2, x2], [yk - 0.01 * span, yk, yk, yk - 0.01 * span], color="black", lw=0.8)
    ax.text((x1 + x2) / 2, yk, stars(p), ha="center", va="bottom", fontsize=10)
ax.set_xlabel("genotype"); ax.set_ylabel("Expression (a.u.)")
ax.set_title("One-way ANOVA, P = %s" % format(omni_p, ".2g"), fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_expression_bar_py.pdf")
print("Reproduced figure.")
