# How pubready works

pubready turns tidy data into a publication-ready figure: drop the data in, decide
on the right statistical test, and get a figure with the significance annotation
drawn on top. It runs as an agent skill and ships a full reproducible record with
every figure.

This document explains the mental model, the workflow, the pieces, and how to
drive it.

## The idea

An agent (Claude Code, Codex, or opencode) sits between the scientist and a small
set of fixed, tested functions. The agent's job is to talk through the analysis
and then call the right function with the right arguments. It does not write the
plotting or the statistics itself. That code already exists, is version
controlled, and is the same every time.

Two consequences follow. The figure a scientist gets is drawn by reviewed code,
not by a model improvising ggplot in the moment. And every figure carries the
exact script, data, and package versions that produced it, so a reviewer can
rerun it.

## The workflow

1. The scientist drops a tidy CSV or Excel file into a project folder.
2. `figkit inspect` reports the columns, their types, group sizes, and candidate
   x and y columns. The data shape is measured, not guessed.
3. The agent chats through the design: how many groups, paired or unpaired,
   parametric assumptions, one factor or several, time-to-event or a matrix. The
   scientist confirms.
4. The agent calls `figkit plot` with the chosen recipe and columns.
5. A recipe checks the assumptions, picks the test (or honors an override), runs
   it, and draws the figure with the statistics on top.
6. The output is one self-contained bundle folder: the figure, the numbers, the
   script that made them, a copy of the input, and a methods paragraph.

If the scientist wants a different look, they edit `plot_config.yaml` in the
bundle and run `figkit render`, which redraws from the edited spec and recomputes
the statistics so the figure and the numbers never drift apart.

## Architecture

Four pieces, kept deliberately separate.

- **`figkit`** (`cli/figkit`) is the stable command surface. The agent only ever
  calls this. It runs everything inside the pinned container.
- **Recipes** (`core/r/R`, `core/python/pubready/recipes`) are the fixed
  functions. Each one owns its test choice, its figure, its stats table, its
  methods paragraph, and the standalone script it emits. Both engines are proper
  packages (`core/r` is a source R package, `core/python` has a `pyproject.toml`)
  with their own test suites, but the code is loaded from the runtime mount, not
  installed into the image.
- **The container** (`container/Containerfile`) is one Podman image with R and
  Python and pinned versions. The environment is pinned by the image digest; the
  code is pinned by the pubready git commit. Both are recorded in every bundle.
- **The bundle writer** takes what a recipe returns and writes the artifact
  folder. It is engine agnostic.

The agent stays thin because the recipes are thick. When a recipe covers the
question, the agent fills in arguments. When no recipe fits, the agent may write
original code, and that code is flagged so the scientist knows it left the tested
path.

## Using figkit

```
figkit inspect --data FILE [--format text|json]
figkit diagnose --recipe NAME --data FILE [options]   # QQ + assumption check, recommend a test
figkit plot --recipe NAME --data FILE [options]
figkit render --config plot_config.yaml [--data FILE]
figkit build          # build or rebuild the container image
figkit shell          # open a shell in the image
figkit test           # run both engines' test suites in the image
figkit version        # print version, image tag/digest, podman, git commit
```

Common options:

| Option | Meaning |
|---|---|
| `--engine r\|python` | which engine draws the figure (default r) |
| `--x`, `--y` | the axis columns for tidy data |
| `--fill`, `--facet` | second and third factors (factorial) |
| `--time`, `--event` | survival columns (event coded 1 = event, 0 = censored) |
| `--covariates "a,b,c"` | Cox model covariates |
| `--annotation FILE` | heatmap column annotation CSV |
| `--test` | force a test instead of `auto` |
| `--paired` | paired design |
| `--geom box\|violin\|bar\|raincloud` | figure geometry |
| `--scale`, `--cluster` | heatmap scaling and clustering |
| `--id`, `--group`, `--label` | the subject, the grouping and the label columns of the recipes that use them |
| `--pd_threshold`, `--pr_threshold` | the two RECIST reference lines of `spider_response` and `waterfall`, +20 and -30 by default |
| `--positive LEVEL` | the case level of the `roc_curve` outcome and the event level of the `calibration` outcome (default: the second level in sorted order) |
| `--ci_low`, `--ci_high`, `--se` | the 95% CI bounds, or the standard error, of each study in `forest_meta` |
| `--measure auto\|ratio\|difference` | the pooling scale of `forest_meta` (a ratio is pooled on the log scale) |
| `--markers "a,b"` | the columns of event times that `swimmer` draws on each bar |
| `--xlab`, `--ylab`, `--title`, `--palette` | labels and colors |
| `--out DIR` | output root (default `pubready_output`) |
| `--cause CODE` | the event code of the cause of interest of `cumulative_incidence` (default 1) |

## The recipes

Each recipe is a named function with a fixed argument surface, present in both
engines. Pick the engine per figure with `--engine`.

| Recipe | Tidy input | Test | Figure |
|---|---|---|---|
| `two_group_compare` | `group, value` | t-test or Mann-Whitney, paired or unpaired | box or violin with a bracket |
| `multi_group_compare` | `group, value` | one-way ANOVA / Welch / Kruskal-Wallis + post hoc | box with brackets for significant pairs, omnibus as a subtitle |
| `factorial_anova` | `f1, f2[, f3], value` | two-way and three-way ANOVA (R also does the aligned rank transform) | grouped box or bar, faceted, effects on top |
| `survival_km` | `time, event, group` | Kaplan-Meier with a log-rank test | step curves, log-rank p, number-at-risk table |
| `cox_forest` | `time, event, covariate…` | Cox proportional hazards | forest plot of hazard ratios with 95% CI |
| `heatmap` | matrix CSV + annotation CSV | per-feature Welch t-test when the annotation has two groups | clustered heatmap, z-scored, group annotation bar |
| `correlation` | `x, y` (two numeric columns) | Pearson / Spearman / Kendall, chosen from normality | scatter with a linear fit, CI band, r and p annotated |
| `correlation_heatmap` | table of numeric variables | pairwise correlation, BH-adjusted | clustered correlation matrix, coefficients + significance stars |
| `upset` | binary membership matrix | descriptive (intersection sizes, no test) | UpSet plot of set intersections |
| `volcano` | DE table (log2FC + p) | thresholds on fold change and p (test is upstream) | EnhancedVolcano-style four-color scatter, top hits labeled |
| `enrichment_dot` | GSEA / ORA result table | descriptive (enrichment is upstream) | dot plot, size = count, color = -log10 adjusted p |
| `paired_compare` | `condition, value, id` | paired t-test or Wilcoxon signed-rank | before/after points with connecting lines and a bracket |
| `proportions` | two categoricals | chi-square, or Fisher when a cell is sparse | 100% stacked bar of proportions with the p |
| `pca` | samples × features + group | PERMANOVA for group separation | PC1/PC2 scatter, 95% ellipses, variance on the axes |
| `spider_response` | `patient, time, value` | descriptive (change from baseline, no test) | one line per patient, RECIST reference lines at +20% and -30% |
| `waterfall` | `patient, time, value`, or one best change per patient | descriptive (best change from baseline, no test) | one sorted bar per patient, RECIST reference lines |
| `dose_response` | `dose, response[, curve]` | four-parameter logistic fit; extra-sum-of-squares F test of one shared EC50 | fitted curves on a log dose axis, mean +/- SEM at each dose, EC50 or IC50 in the legend |
| `roc_curve` | binary outcome, numeric markers | AUC with a DeLong 95% CI; DeLong test of each marker against the first | ROC curves with the Youden point, AUCs in the legend |
| `bland_altman` | `method_a, method_b`, one row per subject | bias and 95% limits of agreement with CIs; regression test of proportional bias | difference against mean, shaded CIs of the bias and the limits |
| `forest_meta` | one row per study: estimate and 95% CI or SE | fixed-effect and DerSimonian-Laird pooling; Cochran's Q, I2, tau2 | forest plot with weighted squares, pooled diamonds and an estimate column |
| `swimmer` | one row per patient: `duration[, status, arm, event times]` | descriptive (no test) | one bar per patient sorted by duration, arrows for ongoing courses, symbols for events |
| `timecourse` | `subject, time, value[, group]` | mixed ANOVA (group between, time within) or repeated-measures ANOVA; Greenhouse-Geisser when Mauchly's test rejects sphericity | mean +/- SEM per group over time, one thin line per subject |
| `cumulative_incidence` | `time, status[, group]`, status 0 = censored, other codes = causes | Aalen-Johansen cumulative incidence of each cause; Gray's test across groups for each cause | one panel per cause, step curves per group, Gray's P in the panel title, number-at-risk table |
| `calibration` | binary outcome, one column of predicted probabilities per model | Brier score, calibration intercept and slope with Wald CIs, Hosmer-Lemeshow test, observed/expected ratio | decile points with Wilson CIs, lowess smooth, ideal line, spike histogram of the predictions |
| `sankey` | `patient, visit, state`, one row per patient visit | transitions table; Bowker's test of symmetry, or McNemar's test, when there are two visits | stacked nodes per visit, ribbons coloured by the state of origin, labels with counts |

### Test selection

`--test auto` walks a short decision tree rather than picking a test from a single
prompt. For the two-group and multi-group recipes it checks normality per group
with the Shapiro-Wilk test and equal variance with Levene's test, then chooses:
parametric when both hold, non-parametric otherwise, with the Welch correction
when variances differ. The `correlation` recipe uses the same normality check to
pick Pearson (both variables normal) or Spearman. The scientist can override any
step by passing `--test` explicitly.

Both engines apply the p-value rule of R's `stats::wilcox.test` to the rank tests.
The exact null distribution is used when each group holds fewer than 50 values and
no value is tied, and the normal approximation with a continuity correction is used
in every other condition. SciPy takes the exact route only for very small samples,
so the Python engine states the rule instead of accepting the SciPy default. The
rank-biserial correlation is computed from the U statistic, and it carries a sign:
it is negative when the first group ranks below the second.

`figkit diagnose` makes that choice visible before any figure is drawn. It writes a
QQ plot for the quantity the test depends on (each group's values, the paired
differences, the model residuals, or each correlated variable) labelled with the
Shapiro-Wilk p, skewness, and kurtosis, and prints the recommended test. The QQ plot
is the tiebreaker a p-value cannot be alone: Shapiro-Wilk is underpowered at small n
and over-rejects at large n, so read the plot and tolerate mild deviation. Every
`figkit plot` bundle carries the same panel as `qc_normality_*`.

### Competing risks, calibration and state transitions

`cumulative_incidence` reads a status column in which 0 marks a censored subject
and each other code marks one event type, or cause. Kaplan-Meier overestimates the
probability of an event when another event can prevent it. The recipe therefore
estimates the cumulative incidence of each cause in each group with the
Aalen-Johansen estimator. Gray's test with rho = 0 compares the groups for each
cause. The R engine calls `cmprsk::cuminc`, and the Python engine restates the
`cinc` and `crstm` routines of cmprsk, so tied times give the same numbers.

- The cause of interest is the code that `--cause` names, and its panel comes
  first. Without `--cause`, the recipe uses code 1, or the first code in sorted
  order when no event has code 1.
- A subject is at risk at time t when its follow-up time is t or more. The
  number-at-risk table below each panel uses this rule.
- The x step is 1, 2, 2.5 or 5 times a power of ten, and it gives about six
  intervals up to the largest follow-up time.
- The stats table gives the cumulative incidence at 25%, 50% and 75% of the
  largest follow-up time. Each time point is rounded to a multiple of half the x
  step. A time point after the last follow-up of a group reads NA, as
  `cmprsk::timepoints` gives it.

`calibration` compares predicted probabilities with the observed outcome of each
subject. Each `--x` column holds the predictions of one model.

- The Hosmer-Lemeshow groups follow `ResourceSelection::hoslem.test`. The breaks
  are the unique type-7 quantiles of the predictions at 0, 0.1, ..., 1, and
  `cut(include.lowest = TRUE)` assigns each subject to a group.
- Tied quantiles merge, so a model with few distinct predictions has fewer than
  10 groups. The degrees of freedom are the count of groups minus 2.
- Two breaks can fall between the same two predictions and leave a group empty.
  `hoslem.test` then returns NaN, and the recipe leaves out the empty group.
- The calibration intercept is the intercept of a logistic regression of the
  outcome with logit(p) as an offset. The calibration slope is the coefficient of
  logit(p) in a logistic regression. Each estimate has a 95% Wald CI.
- The recipe clips each prediction to the range from 1e-8 to 1 - 1e-8 before the
  logit transform. The Brier score and the Hosmer-Lemeshow test use the
  predictions as given.
- The lowess smooth uses f = 2/3, a delta of 0.01 times the range of the
  predictions, and no robustness iterations. Robustness iterations treat the
  events of a binary outcome as outliers and pull the curve toward 0.
- Below the calibration panel, a spike histogram shows the predictions in bins of
  0.01. Events point up and non-events point down.

`sankey` draws the state of each patient at each visit, and a ribbon for each
transition between two consecutive visits.

- Visits sort as numbers when every visit is a number. Other visits sort in
  code-point order.
- States sort as numbers when every state is a number. Other states keep the
  order of their first appearance in the data, so sort the rows to set the order.
- A patient who is missing at a visit takes no part in the transitions into or
  out of that visit.
- With exactly two visits that carry the same set of states, Bowker's test of
  symmetry runs on the patients present at both visits. With two states, this
  test is McNemar's test. Neither engine applies a continuity correction.
- A pair of cells that holds no patient adds nothing to the statistic and no
  degree of freedom. `stats::mcnemar.test` returns NaN for such a table, so both
  engines then compute the statistic from its formula.

### Examples

```bash
# two groups, test chosen automatically
figkit plot --recipe two_group_compare --data data.csv --x group --y volume

# four groups, one-way ANOVA with Tukey post hoc, significant pairs only
figkit plot --recipe multi_group_compare --data data.csv --x genotype --y expression

# two-way ANOVA, grouped box coloured by the second factor
figkit plot --recipe factorial_anova --data data.csv --y response --x genotype --fill treatment

# Kaplan-Meier, drawn with the Python engine
figkit plot --engine python --recipe survival_km --data trial.csv --time months --event status --x arm

# Cox hazard ratios
figkit plot --recipe cox_forest --data trial.csv --time months --event status --covariates "arm,age,sex,stage"

# clustered heatmap with a sample annotation
figkit plot --recipe heatmap --data matrix.csv --annotation samples.csv

# correlation scatter, Pearson or Spearman chosen from the data
figkit plot --recipe correlation --data data.csv --x gene_a --y gene_b

# clustered correlation heatmap across a panel of numeric variables
figkit plot --recipe correlation_heatmap --data cytokines.csv

# UpSet plot from a binary membership matrix
figkit plot --recipe upset --data memberships.csv

# any comparison recipe in the GraphPad Prism style, or as a raincloud or bar
figkit plot --recipe two_group_compare --data data.csv --x group --y volume --theme prism
figkit plot --recipe two_group_compare --data data.csv --x group --y volume --geom raincloud

# volcano plot from a DESeq2 / limma results table
figkit plot --recipe volcano --data de_results.csv

# GSEA / ORA enrichment dot plot
figkit plot --recipe enrichment_dot --data enrichment.csv

# paired before/after with connecting lines
figkit plot --recipe paired_compare --data paired.csv --x condition --y value --id subject

# proportions with a chi-square or Fisher test
figkit plot --recipe proportions --data response.csv --x arm --y response

# PCA scatter with ellipses and a PERMANOVA p
figkit plot --recipe pca --data samples.csv --group group

# four-parameter logistic dose-response curves, one per compound
figkit plot --recipe dose_response --data viability.csv --x concentration --y viability --group compound

# ROC curves of two markers, with the case level stated
figkit plot --recipe roc_curve --data biomarkers.csv --y status --x marker_a,marker_b --positive disease

# Bland-Altman agreement between two assays
figkit plot --recipe bland_altman --data assays.csv --x new_assay --y reference_assay

# meta-analysis forest plot from hazard ratios and their CIs
figkit plot --recipe forest_meta --data trials.csv --y hazard_ratio --ci_low ci_low --ci_high ci_high --label trial

# repeated-measures time course with a mixed ANOVA
figkit plot --recipe timecourse --data growth.csv --x day --y volume --id mouse --group arm

# waterfall of the best change from baseline
figkit plot --recipe waterfall --data lesions.csv --x week --y lesion_mm --id patient --group arm

# swimmer plot of the time on treatment, with response and progression marks
figkit plot --recipe swimmer --data patients.csv --id patient --time months --event ended \
  --group arm --markers response_month,progression_month

# cumulative incidence of relapse (status 1) and of death without relapse (status 2)
figkit plot --recipe cumulative_incidence --data transplant.csv --time months --event status \
  --x arm --cause 1

# calibration of two risk models in a validation cohort
figkit plot --recipe calibration --data validation.csv --y event_1y \
  --x risk_clinical,risk_extended --positive yes

# transitions of the response category across assessments
figkit plot --recipe sankey --data response.csv --id patient --x week --y response
```

## The artifact bundle

Every `figkit plot` run writes one folder. That folder is the reproducible unit.

- `figure_*.pdf`, `.png`, `.svg` — the figure in each format.
- `qc_normality_*.png`, `.pdf` — a QQ panel per quantity the normality check ran on,
  for the recipes that choose a test from normality (absent for the rest).
- `script_*.R` or `script_*.py` — a standalone script that redraws the figure.
- `stats_*.csv` — the test result: statistic, p-value, adjusted p, effect size,
  group sizes.
- `input_*.csv` — a copy of the exact data used (heatmaps also copy the
  annotation), so the bundle does not depend on a file that may change.
- `data_log.md` — the original filename, md5 and sha256 checksums, column types,
  and any cleaning step applied.
- `plot_config.yaml` — the editable spec that drives `figkit render`.
- `session_info.txt` — the full environment: R `sessionInfo()` or a Python
  package list.
- `packages_*.txt` — the package versions for this figure: the ones the script loads
  directly first, then the complete installed set, headed with the R or Python version
  and the container OS.
- `manifest_*.json` — structured provenance: recipe and arguments, the container
  pinned by digest, the pubready git commit, and every package version.
- `methods_*.md` — a manuscript-ready methods paragraph, templated from the same
  manifest so the reported numbers always match the figure.
- `REPRODUCE.md` — the one Podman command to regenerate the figure.

The reproducibility guarantee is precise: the same data in the same pinned
environment gives the same figure visually and the same statistics exactly. It
does not promise a byte-identical PDF, because plotting libraries embed
non-deterministic bits like SVG element ids and PDF creation dates.

## The two engines

Both engines produce the same bundle layout and, where a matching method exists,
the same numbers. The manifest records which engine drew the figure.

| | R engine | Python engine |
|---|---|---|
| plotting | ggplot2 / ggpubr / ggprism | matplotlib / seaborn |
| brackets | rstatix + ggpubr | statannotations |
| ANOVA / post hoc | rstatix | pingouin + scikit-posthocs + statsmodels |
| survival | survival / survminer | lifelines |
| heatmap | ComplexHeatmap | PyComplexHeatmap (seaborn clustermap fallback) |
| correlation | ggplot2 + ggpubr `stat_cor` | seaborn `regplot` |
| correlation heatmap | ComplexHeatmap | seaborn `clustermap` |
| UpSet | ComplexHeatmap `UpSet` | `upsetplot` |
| Prism style (`--theme prism`) | ggprism `theme_prism` | matched matplotlib style |
| volcano / labels | ggplot2 + ggrepel | matplotlib + adjustText |
| raincloud geom | ggdist half-eye | matplotlib half-violin |
| PCA / PERMANOVA | prcomp + exact or seeded permutation | numpy SVD + exact or seeded permutation |
| dose-response (4PL) | `stats::nls` | `scipy.optimize.curve_fit` |
| ROC / DeLong | pROC | the rules of pROC, restated with numpy |
| meta-analysis | metafor `rma` | statsmodels `combine_effects` |
| mixed ANOVA | split-plot sums of squares, the GG epsilon and Mauchly's test of car | the same formulas with numpy |
| competing risks | cmprsk `cuminc` | the `cinc` and `crstm` routines of cmprsk, restated with numpy |
| calibration | `stats::glm`, `stats::lowess` | statsmodels `GLM` and `lowess`, the groups of `hoslem.test` restated |
| Sankey diagram | ggplot2 polygons, `stats::mcnemar.test` | matplotlib polygons, statsmodels `SquareTable` |

Where Python has no faithful equivalent, pubready routes you to R rather than
running something different. The one current case is the aligned rank transform
for non-parametric factorial ANOVA: `--test art --engine python` errors and points
to `--engine r`. The Python `heatmap` draws the clustered matrix and the group
bar but keeps the per-feature significance stars in the stats table rather than on
the figure. The `pca` PERMANOVA enumerates every labeling of the samples when there
are 5000 or fewer, so the P value is exact and the same in both engines. A larger
design uses 999 seeded permutations, and each engine draws its own, so the two P
values differ by the Monte Carlo error.

## Setup

pubready needs Podman and one built image.

```bash
figkit build          # builds localhost/pubready:1.0.0 from container/Containerfile
```

The image carries R (ggpubr, rstatix, ggprism, ggrepel, ggdist, survival,
survminer, ARTool, ComplexHeatmap, pROC, metafor) and Python (matplotlib, seaborn,
statannotations, pingouin, statsmodels, scikit-posthocs, lifelines,
PyComplexHeatmap, upsetplot, adjustText), all pinned. The
`core/` scripts are mounted at runtime, so editing a recipe does not need a
rebuild; only changing the dependencies does.

## Testing and versioning

Both engines carry a unit suite. `figkit test` runs them inside the image
against the mounted source, so the tested code is exactly what figkit runs: the
R suite (`core/r/tests/testthat`) loads the source package with `pkgload`; the
Python suite (`core/python/tests`) runs under `pytest`. The end-to-end
`tests/smoke_test.sh` drives every recipe on both engines.

The public benchmarks in `benchmarks/` check the numbers. Each recipe runs on
a published dataset, for example the `sleep` data of Student (1908), the NCCTG
lung cancer trial, the BCG vaccine trials or the PEFR data of Bland and Altman
(1986). Both suites compare the result with a reference value that the
canonical R function computed, such as `stats::t.test`, `survival::coxph`,
`pROC::roc` or `metafor::rma`. `tests/benchmark_test.sh` runs the same cases
end to end through `figkit plot`. See `benchmarks/README.md`.

The version lives in the root `VERSION` file and moves in lockstep across
`core/r/DESCRIPTION`, `core/python/pubready/version.py`, and the image tag.
`CHANGELOG.md` records every release. "Breaking" means a change to the figkit
command surface or to figure output, not just an internal refactor.

## Adding a recipe

A recipe is one function, `recipe_<name>`, in `core/r/R/<name>.R` and
`core/python/pubready/recipes/<name>.py`. It receives the raw data frame and the
spec, and returns:

- the figure (a ggplot in R, a matplotlib figure in Python),
- a tidy stats table,
- a test-metadata dict for the manifest,
- a methods paragraph,
- a builder for the standalone script,
- the cleaned data it used and the cleaning steps.

The bundle writer does the rest. Because recipes own their own methods text and
script, new ones plug in without touching the writer. Register the recipe by
adding an example dataset, a smoke case in `tests/smoke_test.sh`, a reference
bundle in `tools/make_examples/generate_expected.sh`, and unit tests under
`core/r/tests/testthat` and `core/python/tests`. Add a public dataset and its
reference values to `tools/make_benchmarks/make_benchmark_data.R`, and a case
to both benchmark suites and to `tests/benchmark_test.sh`.

## Design principles

- Thin agent, thick recipes. The model reasons and fills arguments; fixed code
  does the analysis.
- Measure the data before choosing a test.
- Put the statistic on the figure, not just in a table.
- Ship the record with the result. Every figure is reproducible or it does not
  ship.
- When an engine cannot do something faithfully, say so and route to the one that
  can.
