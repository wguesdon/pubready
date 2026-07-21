# Assumption checks

pubplot measures the parametric assumptions before it chooses a test, and records
what it measured. This file documents what is checked, with which test, at which
threshold, and how the result steers the choice. It pairs with
[decision_tree.md](decision_tree.md).

## What `figkit inspect` reports
`figkit inspect` describes the shape of the data so you can pick x and y and reason
about the design. For the whole file it reports the row and column counts. For each
column it reports the type, the missing count, the number of unique values, the
category levels when there are 20 or fewer, and the min, max, and mean of numeric
columns. Group sizes come from the levels of the grouping column.

`inspect` does not run the normality or variance tests. Those are computed by the
recipe at analysis time, under `--test auto`, on the exact rows being analyzed, and
their p-values are written into the stats table and the manifest. Measuring inside
the recipe is deliberate: the check runs on the same data as the test, so the two
cannot disagree.

## Normality: Shapiro-Wilk
Normality is tested with Shapiro-Wilk. A group or variable is treated as normal when
its Shapiro-Wilk P is above 0.05. Where the test is applied depends on the recipe:

- `two_group_compare`, `multi_group_compare` — each group's values.
- `correlation` — each of the two variables.
- `paired_compare` — the within-subject differences.
- `factorial_anova` — the residuals of the fitted ANOVA model.

Sample-size guards. With fewer than 3 finite values the test cannot run, and the
data is treated as not meeting the normality assumption, so the recipe takes the
nonparametric branch. Above 5000 values Shapiro-Wilk is unreliable, so it is skipped
and the data is treated as normal. Both guards are recorded.

## Equal variance
For the continuous group comparisons the recipe also checks whether the groups share
a variance, and treats them as equal when the test's P is above 0.05:

- `two_group_compare` — F test (`var.test`) on the two groups. Equal variance keeps
  the Student's t-test; unequal variance switches to Welch.
- `multi_group_compare` — Levene's test across all groups. Equal variance keeps
  one-way ANOVA with Tukey; unequal variance switches to Welch's ANOVA with
  Games-Howell.

When a variance test cannot be computed, the recipe treats the variances as equal
rather than failing. `paired_compare` and `correlation` key off normality alone, and
`factorial_anova` keys off residual normality rather than group variance.

## Seeing the assumption: the QQ plot

A p-value is a poor sole judge of normality. Shapiro-Wilk is underpowered at small n,
so it passes visibly non-normal small samples, and it over-rejects at large n, where
the central limit theorem already makes the mean-based tests robust. So pubplot shows
the assumption rather than only asserting it. `figkit diagnose` writes a QQ plot for the
checked quantity, labelled with the Shapiro-Wilk p, skewness, and excess kurtosis, and
prints the recommended test without drawing the figure, for the choose-a-test
conversation. Every `figkit plot` bundle carries the same panel as `qc_normality_*`.
Read the QQ: points near the line support normality, and mild, judged deviation is
usually acceptable.

## How the checks steer the choice

| Recipe | Normal, equal variance | Normal, unequal variance | Not normal |
|---|---|---|---|
| `two_group_compare` | Student's t-test | Welch t-test | Mann-Whitney U |
| `multi_group_compare` | ANOVA + Tukey | Welch ANOVA + Games-Howell | Kruskal-Wallis + Dunn (BH) |
| `paired_compare` | paired t-test | paired t-test | Wilcoxon signed-rank |
| `correlation` | Pearson | Pearson | Spearman |
| `factorial_anova` | Type II ANOVA | Type II ANOVA | aligned rank transform |

`factorial_anova` reads the last column from the normality of the model residuals,
not from a group-variance test.

## Where the numbers land
Every run records the checks it made. The stats table carries the Shapiro-Wilk
p-values (per group, per variable, or on the differences) and the variance-test
p-value. The manifest repeats them under an `assumptions` block next to the test
that was chosen. A reviewer sees the check, the threshold, and the decision without
rerunning anything.

## Overriding a check
`--test` forces the method regardless of what the checks say. The checks still run
and are still recorded, so the bundle shows both what the data suggested and what was
actually run. Use this when domain knowledge outweighs a borderline assumption test;
the methods paragraph will state the forced choice.
