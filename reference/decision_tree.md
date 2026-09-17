# Statistical decision tree

This is the reference the agent walks with the scientist to choose a test. Every
rule here matches the recipe code, so the doc and the analysis cannot drift. For
each design it names the default that `--test auto` reaches, the assumption checks
behind that default, and the flag to override it.

The order of questions:

1. What is the outcome? Continuous, a proportion or category, or time to an event.
2. How many groups, and are they paired or repeated?
3. Do the parametric assumptions hold? Measured, not assumed. See
   [assumptions.md](assumptions.md).
4. Is there more than one comparison? If so, correct for it.

`--test auto` measures the assumptions and picks the default. The scientist can
force any test with `--test`, and the recipe records both the choice and the
checks that led to it.

## Continuous outcome

### Two groups → `two_group_compare`
Check normality of each group (Shapiro-Wilk) and equal variance (Levene's test).

- Both groups normal, equal variance → Student's two-sample t-test.
- Both groups normal, unequal variance → Welch two-sample t-test.
- Either group not normal → Mann-Whitney U test.
- Paired design (`--paired`): paired t-test when the differences are normal, else
  the Wilcoxon signed-rank test. Paired data that carries a subject id is better
  served by `paired_compare`, which draws the connecting lines.

Effect size is Cohen's d for the t-tests and rank-biserial r for the rank tests.
Force with `--test student_t | welch_t | wilcoxon`.

### Three or more groups → `multi_group_compare`
Check normality of every group (Shapiro-Wilk) and equal variance across groups
(Levene's test).

- All groups normal, equal variance → one-way ANOVA, post hoc Tukey's HSD.
- All groups normal, unequal variance → Welch's one-way ANOVA, post hoc
  Games-Howell.
- Any group not normal → Kruskal-Wallis test, post hoc Dunn's test with
  Benjamini-Hochberg adjustment.

The omnibus test is the figure subtitle. Only pairwise comparisons that reach
significance (adjusted P < 0.05) get a bracket; the full pairwise table is saved.
Force with `--test anova | welch | kruskal`. This recipe also handles two groups.

### Paired, two conditions with a subject id → `paired_compare`
Check normality of the within-subject differences (Shapiro-Wilk).

- Differences normal → paired t-test (effect size Cohen's dz).
- Otherwise → Wilcoxon signed-rank test (rank-biserial r).

A line connects each subject's two points. Force with `--test t | wilcoxon`.

### Two or three factors → `factorial_anova`
Type II ANOVA by default. Under `--test auto`, when the fitted model's residuals
are not normal (Shapiro-Wilk on the residuals, P < 0.05) the recipe switches to
the aligned rank transform (ARTool), the nonparametric factorial path. Every
effect's p-value is drawn on the figure. Force with `--test anova | art`. The
aligned rank transform is R only; `--test art --engine python` errors and routes
you to R.

## Association between two continuous variables → `correlation`
Check normality of each variable (Shapiro-Wilk).

- Both variables normal → Pearson correlation (r).
- Otherwise → Spearman rank correlation (rho).
- Kendall rank correlation (tau) only when forced.

Scatter with a linear fit and a 95% confidence band, the coefficient and p
annotated. Force with `--test pearson | spearman | kendall`.

## Proportion or category → `proportions`
Raw rows with a group column and an outcome column. Build the contingency table,
then:

- Any expected cell below 5 → Fisher's exact test.
- Otherwise → Pearson's chi-squared test (no continuity correction).

A 2x2 table also reports an odds ratio. The figure is a 100% stacked bar of the
outcome proportions per group. Force with `--test chisq | fisher`.

## Time to event

### Group comparison → `survival_km`
Kaplan-Meier curves with censoring ticks and a number-at-risk table. Groups are
compared with the log-rank test, its p-value on the plot. Input `--time --event
--x`, with the event column coded 1 = event, 0 = censored.

### Adjusted hazard → `cox_forest`
Cox proportional-hazards model over one or more covariates, drawn as a forest of
hazard ratios with 95% confidence intervals. A global cox.zph test checks the
proportional-hazards assumption and its p-value is reported. Input `--time
--event --covariates "a,b,c"`.

## Descriptive figures (no test chosen from assumptions)
These recipes show structure rather than run an inferential test:

- `heatmap` — clustered heatmap; adds a per-feature Welch t-test only when the
  annotation has exactly two groups (stars on the R figure; the Python engine
  keeps the stars in the stats table).
- `correlation_heatmap` — pairwise correlation across the numeric columns,
  clustered, with Benjamini-Hochberg-adjusted significance stars.
- `volcano`, `enrichment_dot`, `upset` — descriptive; any differential or
  enrichment test is upstream of the figure.
- `pca` — PCA scatter with 95% ellipses; group separation is tested with a seeded
  PERMANOVA and its p is reported.
- `spider_response` — the change from baseline of each patient over time, one
  line per patient, with the RECIST 1.1 reference lines at +20% and -30%. The
  stats table carries the best and the worst change of each patient and the
  category those two put the patient in. That reads the change from baseline
  alone, and a full RECIST 1.1 assessment also reads the nadir, the non-target
  lesions and any new lesion.

## Multiple comparisons
Correct whenever there is more than one comparison. The multi-group post-hoc tests
carry their own correction: Tukey's HSD and Games-Howell are inherently adjusted,
and Dunn's test uses Benjamini-Hochberg. The correlation heatmap adjusts across the
matrix with Benjamini-Hochberg. A single two-group comparison needs no correction,
so `--p_adjust` defaults to none there.

## Overriding
The tree gives the default; the scientist decides. When the choice is not obvious, run
`figkit diagnose` first: it shows the QQ plot and the assumption numbers and prints the
recommended test, so the scientist can confirm it or force another with `--test`.
Forcing a test still runs and records the assumption checks, and the methods paragraph
states exactly what was run.
