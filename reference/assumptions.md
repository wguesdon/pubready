# Assumption checks

pubready measures the parametric assumptions before it chooses a test, and records
what it measured. This file documents what is checked, with which test, at which
threshold, and how the result steers the choice. It pairs with
[decision_tree.md](decision_tree.md).

## What `figkit inspect` reports
`figkit inspect` describes the shape of the data so you can pick x and y and reason
about the design. For the whole file it reports the row and column counts. For each
column it reports the type, the missing count, the number of unique values, the
category levels when there are 20 or fewer, and the min, max, and mean of numeric
columns. Each level is in quotes with its count of rows, so the group sizes are on
the page. A numeric column with 2 to 6 distinct values is a candidate group, marked
as a numeric code, and a column whose values repeat across the levels of another
column is a possible subject id. Both engines print the same words.

`inspect` does not run the normality or variance tests. Those are computed by the
recipe at analysis time, under `--test auto`, on the exact rows being analyzed, and
their p-values are written into the stats table and the manifest. Measuring inside
the recipe is deliberate: the check runs on the same data as the test, so the two
cannot disagree.

## Normality: Shapiro-Wilk
Normality is tested with Shapiro-Wilk. A group or variable is treated as normal when
its Shapiro-Wilk P is above 0.05. Where the test is applied depends on the recipe:

- `two_group_compare`, `multi_group_compare`, `feature_panel` — each group's values.
- `two_group_compare --paired`, `paired_compare` — the within-pair differences.
- `correlation` — each of the two variables.
- `factorial_anova` — the residuals of the fitted ANOVA model.

Sample-size guards. With fewer than 3 finite values the test cannot run, and the
data is treated as not meeting the normality assumption. A constant sample is also
not normal, in both engines. Above 5000 values Shapiro-Wilk is unreliable, so it is
skipped and the data is treated as normal. The guards are recorded.

The attainable P. A failed normality check sends `--test auto` to a rank test, but
at a small n a rank test cannot reach P < 0.05: 3 against 3 gives at best P = 0.10.
Before it takes the rank test, the recipe computes the smallest P that the rank
test can reach at the group sizes. When that P is above 0.05, the parametric test
stays, `test_meta` records the attainable P, and the methods text says why. The rule
applies to `two_group_compare`, `paired_compare`, `multi_group_compare` and
`feature_panel`.

## Equal variance
For the continuous group comparisons the recipe also checks whether the groups share
a variance, and treats them as equal when the test's P is above 0.05:

- `two_group_compare` — Levene's test on the two groups, centered on the median.
  The P goes to the stats table, but it does not choose the test: two normal
  groups get Welch's t-test, which does not assume equal variances.
  `--test student_t` forces Student's t-test.
- `multi_group_compare` — Levene's test across all groups. Equal variance keeps
  one-way ANOVA with Tukey; unequal variance switches to Welch's ANOVA with
  Games-Howell.

When a variance test cannot be computed, the recipe treats the variances as equal
rather than failing. `paired_compare` and `correlation` key off normality alone, and
`factorial_anova` keys off residual normality rather than group variance.
`feature_panel` runs no variance test: it uses Welch's t-test for two groups and
the one-way ANOVA for three groups or more.

## Seeing the assumption: the QQ plot

A p-value is a poor sole judge of normality. Shapiro-Wilk is underpowered at small n,
so it passes visibly non-normal small samples, and it over-rejects at large n, where
the central limit theorem already makes the mean-based tests robust. So pubready shows
the assumption rather than only asserting it. `figkit diagnose` writes a QQ plot for the
checked quantity, labelled with the Shapiro-Wilk p, skewness, and excess kurtosis, and
prints the recommended test without drawing the figure, for the choose-a-test
conversation. Each `figkit plot` bundle of a recipe that chooses its test from a
normality check carries the same panel as `qc_normality_*`.
Read the QQ: points near the line support normality, and mild, judged deviation is
usually acceptable.

## How the checks steer the choice

| Recipe | Normal, equal variance | Normal, unequal variance | Not normal |
|---|---|---|---|
| `two_group_compare` | Welch t-test | Welch t-test | Mann-Whitney U |
| `multi_group_compare` | ANOVA + Tukey | Welch ANOVA + Games-Howell | Kruskal-Wallis + Dunn (BH) |
| `paired_compare` | paired t-test | paired t-test | Wilcoxon signed-rank |
| `correlation` | Pearson | Pearson | Spearman |
| `factorial_anova` | Type II ANOVA | Type II ANOVA | aligned rank transform (R); Python stops and routes to R |
| `feature_panel` | Welch t-test or one-way ANOVA | Welch t-test or one-way ANOVA | Mann-Whitney U or Kruskal-Wallis |

`factorial_anova` reads the last column from the normality of the model residuals,
not from a group-variance test. In every row, the "Not normal" column becomes the
parametric test when the rank test cannot reach P < 0.05 at the group sizes.

## Where the numbers land
Every run records the checks it made. The stats table of `two_group_compare`
carries the Shapiro-Wilk P of each group, the Levene P and, under `--paired`, the P
of the differences. The table of `paired_compare` carries the P of the differences,
the table of `correlation` the P of each variable, and the table of
`factorial_anova` the P of the residuals. The manifest records the checks under an
`assumptions` block next to the test that was chosen. For `multi_group_compare`,
the manifest gives whether all groups are normal and the Levene P, and the QC panel
shows each group. A reviewer sees the check, the threshold, and the decision
without rerunning anything.

## Overriding a check
`--test` forces the method regardless of what the checks say. The checks still run
and are still recorded, so the bundle shows both what the data suggested and what was
actually run. Use this when domain knowledge outweighs a borderline assumption test;
the methods paragraph will state the forced choice.

## Advice before the test

`figkit advise` reuses the inspect functions for data facts. It does not run
Shapiro-Wilk or another statistical test. Its next action is diagnose when the
recipe chooses its test from normality. The recipe applies the attainable-P
rule to the usable observations.
