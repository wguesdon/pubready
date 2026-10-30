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

### Repeated measures over time → `timecourse`
Each subject is measured at every time point, and a subject can belong to a group.
The test is a mixed ANOVA with the group between subjects and time within subjects,
or a one-way repeated-measures ANOVA when there is no group. The ANOVA needs
complete subjects, so the recipe drops a subject with a missing visit and records
the drop.

- Mauchly's test rejects sphericity (P < 0.05) → the time effect and the
  interaction carry the Greenhouse-Geisser correction.
- Otherwise → no correction.

Both P values are in the stats table. Force with `--test gg | uncorrected`. Two
time points with a subject id are better served by `paired_compare`. The
ANOVA treats time as a factor, so it does not test a rate of change. A question
about a growth rate needs a mixed model with time as a number, which no recipe
fits. The QC panel shows the residuals. Growth data such as tumour volumes
often give residuals that spread with the mean. In that condition, analyse the
log of the value.

### Response against dose → `dose_response`
A four-parameter logistic curve is fitted to each curve (`--group`). The midpoint
is an EC50 when the curve rises and an IC50 when it falls, with a 95% CI from the
fit. With two or more curves, an extra-sum-of-squares F test asks whether one
shared midpoint fits as well as one midpoint per curve. A dose of zero or below
cannot go on a log axis, so the recipe drops it and records the drop.

## Agreement and accuracy

### Two methods that measure the same quantity → `bland_altman`
Do not use a correlation to show agreement: two methods can correlate well and
still disagree. The Bland-Altman plot gives the bias and the 95% limits of
agreement, each with a confidence interval. A regression of the difference on
the mean tests for a proportional bias. The difference is `--x` minus `--y`,
so put the new method in `--x` and the reference method in `--y`. The limits
assume normal differences, so the QC panel shows their QQ plot.

### A marker that discriminates cases from controls → `roc_curve`
The outcome must have exactly two levels. Confirm the case level with the
scientist and pass `--positive` when the default is wrong. The default is the
second level in sorted order. The recipe reports the AUC with a DeLong 95% CI and
the threshold that maximises the Youden index. Two or more markers on the same
subjects are compared with the first marker by DeLong's test for correlated ROC
curves.

### Predicted risks against the observed outcome → `calibration`
Use this recipe when a model gives each subject a predicted probability and the
question is whether the predictions match the observed event rates. The ROC curve
measures discrimination, and calibration measures agreement. The outcome must
have exactly two levels, and `--positive` names the event level as for
`roc_curve`. Each `--x` column holds the predictions of one model, between 0 and
1. The recipe reports the Brier score, the calibration intercept and the
calibration slope with 95% Wald CIs, and the Hosmer-Lemeshow test with 10 groups.
An intercept near 0 and a slope near 1 show good calibration. A slope below 1
shows predictions that are too extreme, which is the usual sign of overfitting.

## Combining studies → `forest_meta`
Each row is one study with its estimate and 95% CI, or its standard error. A
ratio (odds, risk or hazard ratio) is pooled on the log scale. `--measure auto`
picks the ratio scale when every value is positive and the CIs are more symmetric
on the log scale. The recipe reports the fixed-effect and the DerSimonian-Laird
random-effects estimates. Read the heterogeneity (Cochran's Q, I2 and tau2) before
you report a pooled value. When I2 is high, report the random-effects estimate.

## Association between two continuous variables → `correlation`
Check normality of each variable (Shapiro-Wilk).

- Both variables normal → Pearson correlation (r).
- Otherwise → Spearman rank correlation (rho).
- Kendall rank correlation (tau) only when forced.

A scatter with the coefficient and p annotated. A Pearson scatter also has a
linear fit and a 95% confidence band. A Spearman or a Kendall scatter has no
line, because a rank correlation measures a monotonic association and not a
straight line. Force with `--test pearson | spearman | kendall`.

## Proportion or category → `proportions`
Raw rows with a group column and an outcome column. Build the contingency table,
then:

- Any expected cell below 5 → Fisher's exact test.
- Otherwise → Pearson's chi-squared test (no continuity correction).

A 2x2 table also reports an odds ratio. When a third variable splits the data
into strata, do not pool the table. See "Designs that no recipe fits". The
figure is a 100% stacked bar of the outcome proportions per group. Force with
`--test chisq | fisher`.

### The same patients in a category at each visit → `sankey`
Use this recipe when each patient has a category at two or more visits, for
example a response category or a toxicity grade. Give one row for each patient and
visit with `--id --x --y`. The figure shows the transitions between consecutive
visits. With exactly two visits that carry the same categories, Bowker's test of
symmetry asks whether the patients move in one direction more than in the other.
With two categories, this test is McNemar's test. Neither test applies a
continuity correction. With three visits or more, the figure is descriptive.

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

### Competing events → `cumulative_incidence`
Use this recipe when another event can prevent the event of interest, for example
a death before a relapse. Kaplan-Meier treats the competing event as censoring and
overestimates the probability of the event of interest. Code the status as 0 for
a censored subject and as one code for each event type, and name the cause of
interest with `--cause`. The recipe estimates the cumulative incidence of each
cause in each group with the Aalen-Johansen estimator, and Gray's test compares
the groups for each cause. Input `--time --event [--x --cause]`.

## Descriptive figures (no test chosen from assumptions)
These recipes show structure rather than run an inferential test:

- `heatmap` — clustered heatmap; adds a per-feature Welch t-test only when the
  annotation has exactly two groups (stars on the R figure; the Python engine
  keeps the stars in the stats table).
- `correlation_heatmap` — pairwise correlation across the numeric columns,
  clustered, with Benjamini-Hochberg-adjusted significance stars.
- `volcano`, `enrichment_dot`, `upset` — descriptive; any differential or
  enrichment test is upstream of the figure.
- `pca` — PCA scatter with 95% normal-theory ellipses for each group of four
  samples or more. Group separation is tested with a PERMANOVA. The test is exact
  when the design has 5000 labelings or fewer, and it uses 999 seeded permutations
  otherwise.
- `spider_response` — the change from baseline of each patient over time, one
  line per patient, with the RECIST 1.1 reference lines at +20% and -30%. The
  stats table carries the best and the worst change of each patient and the
  category those two put the patient in. That reads the change from baseline
  alone, and a full RECIST 1.1 assessment also reads the nadir, the non-target
  lesions and any new lesion.
- `waterfall` draws the best change from baseline of each patient as one sorted
  bar, with the same RECIST reference lines and the same categories as
  `spider_response`. It also accepts one precomputed change per patient.
- `swimmer` draws the time on treatment of each patient as one sorted bar, with an
  arrow for a course that is still ongoing and a symbol for each event time named
  in `--markers`. The stats table gives the median duration of each arm.

## Multiple comparisons
Correct whenever there is more than one comparison. The multi-group post-hoc tests
carry their own correction: Tukey's HSD and Games-Howell are inherently adjusted,
and Dunn's test uses Benjamini-Hochberg. The correlation heatmap adjusts across the
matrix with Benjamini-Hochberg. A single two-group comparison needs no correction,
so `--p_adjust` defaults to none there.

## Designs that no recipe fits
No recipe fits the designs below. Tell the scientist, and do not force the
nearest recipe, because a pooled test can give the wrong answer.

- A third variable splits the data into strata, such as a department or a
  center. Pooled, the UC Berkeley admissions give an odds ratio of 1.84 against
  women. Within the departments, the Mantel-Haenszel odds ratio is 0.90. The
  analysis needs a Mantel-Haenszel test or a regression model.
- The comparison must adjust for a covariate, such as age. The analysis needs
  ANCOVA or a regression model.
- The question is about a rate of change, or the repeated measures have
  missing visits. The analysis needs a mixed model.
- Each group is compared with one control only. `multi_group_compare` compares
  all pairs, and pubready has no Dunnett test.

## Overriding
The tree gives the default; the scientist decides. When the choice is not obvious, run
`figkit diagnose` first: it shows the QQ plot and the assumption numbers and prints the
recommended test, so the scientist can confirm it or force another with `--test`.
Forcing a test still runs and records the assumption checks, and the methods paragraph
states exactly what was run.
