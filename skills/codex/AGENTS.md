# pubready (Codex adapter)

When the user has tidy data and wants a publication figure, use pubready. It makes
a publication-ready figure with the statistical test drawn on top. A small
library of tested, parameterized recipes
does the plotting and the statistics. Your job is to choose the right test with
the user and call the recipe. Do not write ggplot or matplotlib by hand. Write
original plotting or stats code only when no recipe fits, and say so when you do.

## The one command surface

Everything goes through `./cli/figkit`, which runs the pinned Podman image. Never
call R or Python directly.

In a clone of pubready, call `./cli/figkit` from the root of the repository.
After a uv install, call `figkit`, or `uv run figkit` in a uv project. A uv
install has no repository folder. `figkit where` prints the folder that holds
`reference/` and `docs/`. Read those files from that folder.
`figkit` mounts only the current folder into the container. The data file must
be in that folder or below it.

```
figkit build                                  # once, builds the image
figkit doctor                                 # checks Podman, the image and the mounts, prints each fix
figkit inspect --data FILE [--format json]    # columns, types, group sizes, candidate x/y
figkit advise --data FILE --format json ...  # recipe, command and unanswered questions
figkit diagnose --recipe NAME --data FILE ... # QQ + assumption check, recommend a test (no figure)
figkit plot --recipe NAME --data FILE ...     # draw + test, writes a bundle
figkit render --config plot_config.yaml       # redraw from an edited spec, recompute stats
```

If a command reports that the image is missing, run `figkit build` first. If a
command reports a Podman problem, run `figkit doctor` and show its output to the
user. Do not install Podman, and do not run `sudo`. The user runs each fix that
needs administrator rights.

Each recipe reads its own options. An option that the recipe does not read stops
the run, and the message names the recipes that read it. `core/recipe_options.txt`
lists the options of each recipe. Both engines read a CSV file with one rule: only
an empty cell or "NA" is missing, the spaces around a value are removed, and a
number with a comma, such as "85,45" or "4,380", stops the run. Change a decimal
comma to "." and remove each thousands mark before the run.

In a clone, `figkit plot` and `figkit render` stop when the repository holds an
untracked file that git does not ignore, because the bundle records the commit.
Write each helper file, such as a reshaped CSV, under `pubready_output/`, which
git ignores.

## Workflow

1. Run `figkit inspect` on the data. Read the real column types and group sizes.
   Do not guess them. `inspect` also names the candidate groups, numeric codes
   among them, and each column that looks like a subject id.
2. Run `figkit advise --data FILE --format json` with the known plot flags and
   `--answer KEY=VALUE` answers. Ask the scientist each unanswered question.
   Never supply a scientist answer yourself. Repeat advise until it returns
   `complete`, then follow `next_step`. Resolve a conflict before any plot.
   Exit code 5 means that no recipe supports the design.
3. If no test is named, run `figkit diagnose --recipe … --data … [--x --y …]`. It
   draws a QQ plot with the Shapiro-Wilk p, skewness, and kurtosis for the quantity
   the test depends on, and prints a recommended test. Show the scientist the
   recommendation and the QQ, and let them accept it or force another with `--test`.
   Shapiro is a default, not a verdict (underpowered at small n, over-rejects at
   large n), so the QQ plot is the tiebreaker.
4. Call `figkit plot` with the recipe and columns. Leave `--test auto` unless a test
   is forced. When the recipe chooses its test from a normality check, the bundle
   also includes the QC panel (`qc_normality_<ts>.png`).
   Pass `--llm-provider`, `--llm-model`, `--agent-harness` and
   `--agent-harness-version` with your own values, so that the report names the
   model and the harness that made the figure.
5. Point them at the bundle folder under `pubready_output/`: the figure, the
   standalone script, the stats table, the input copy, and a methods paragraph.
6. To restyle, edit `plot_config.yaml` in the bundle and run `figkit render`. It
   recomputes the statistics, so the figure and the numbers never drift apart.
   It uses the engine of the config unless you pass `--engine`, and it finds the
   bundle input in the folder of the config.

## Recipes

| Recipe | Tidy input | Columns you pass | Test |
|---|---|---|---|
| `two_group_compare` | group, value | `--x --y [--levels --paired]` | Welch t-test or Mann-Whitney, paired t-test or Wilcoxon signed-rank with `--paired` |
| `multi_group_compare` | group, value | `--x --y [--levels]` | one-way ANOVA / Welch / Kruskal-Wallis + post hoc |
| `factorial_anova` | f1, f2[, f3], value | `--x --y --fill [--facet --levels]` | two-way / three-way ANOVA (R also does the aligned rank transform) |
| `survival_km` | time, event, group | `--time --event --x` | Kaplan-Meier + log-rank |
| `cox_forest` | time, event, covariates | `--time --event --covariates "a,b,c"` | Cox proportional hazards |
| `heatmap` | matrix CSV + annotation CSV | `--data --annotation [--test welch_t]` | clustered heatmap; a Welch t-test for each feature only with `--test welch_t` and two annotation groups |
| `correlation` | x, y (two numeric columns) | `--x --y` | Pearson / Spearman / Kendall, chosen from normality; not for the agreement of two methods, use `bland_altman` |
| `correlation_heatmap` | table of numeric variables | `--data` | pairwise correlation, clustered, BH-adjusted stars |
| `upset` | membership matrix of 1/0, TRUE/FALSE or yes/no | `--data [--top_n]` | descriptive intersection sizes, the 15 largest drawn by default (no test) |
| `volcano` | DE table (log2FC + p) | `--data [--x --y --label]` | EnhancedVolcano-style scatter, top hits labeled |
| `enrichment_dot` | GSEA / ORA result table | `--data [--label]` | dot plot, size = count, color = -log10 adjusted p; a raw p is marked as not adjusted |
| `paired_compare` | condition, value, id | `--x --y --id [--levels]` | paired t-test or Wilcoxon, connecting lines |
| `proportions` | two categoricals | `--x --y [--levels]` | chi-square / Fisher, 100% stacked bar |
| `pca` | samples × features + group | `--data --group [--features "a,b" --id --label]` | PC1/PC2 scatter, 95% ellipses, PERMANOVA; the `--id` and `--label` columns are not features |
| `spider_response` | patient, time, value | `--x --y --id [--group]` | descriptive change from baseline, RECIST lines (no test) |
| `waterfall` | patient, time, value (or one best change per patient) | `--y --id [--x --group]` | descriptive best change from baseline, sorted bars, RECIST lines (no test) |
| `dose_response` | dose, response[, curve] | `--x --y [--group]` | four-parameter logistic fit, EC50 or IC50 with 95% CI, shared-EC50 F test |
| `roc_curve` | binary outcome, numeric markers | `--y --x "m1,m2" [--positive]` | ROC curves, AUC with DeLong CI, Youden threshold, DeLong test between markers |
| `bland_altman` | method A, method B (one row per subject) | `--x --y` | bias and 95% limits of agreement with CIs, proportional bias test |
| `forest_meta` | one row per study: estimate and 95% CI or SE | `--y --ci_low --ci_high [--se --label --measure]` | fixed-effect and DerSimonian-Laird pooling, Q, I2, tau2 |
| `timecourse` | subject, time, value[, group] | `--x --y --id [--group --levels]` | mixed (or repeated-measures) ANOVA, Greenhouse-Geisser when Mauchly rejects |
| `swimmer` | one row per patient: duration[, status, arm, event times] | `--id --time [--event --group --markers "a,b"]` | descriptive bar per patient, arrows for ongoing courses (no test) |
| `distribution` | one numeric variable[, group] | `--y [--x]` | descriptive histogram, kernel density and normal fit per group, Shapiro-Wilk (no comparison) |
| `composition` | sample, category[, count, condition] | `--id --x [--y --group]` | 100% stacked bar per sample, per-category Mann-Whitney or Kruskal-Wallis, BH-adjusted (exploratory) |
| `estimation_plot` | group, value[, subject] | `--x --y [--ref --paired --id --levels]` | mean difference against the reference with a 95% BCa bootstrap CI, Welch or paired t-test; Gardner-Altman or Cumming plot |
| `enzyme_kinetics` | substrate, rate[, condition] | `--x --y [--group]` | Michaelis-Menten fit, Vmax and Km with 95% CIs, F tests of one shared Km and one shared Vmax |
| `standard_curve` | concentration (empty for an unknown), signal[, sample] | `--x --y [--id --test auto\|linear\|4pl]` | 4PL or linear standard curve, unknowns read off the curve, no extrapolation, recovery of the standards |
| `cumulative_incidence` | time, status (0 = censored, other codes = causes)[, group] | `--time --event [--x --cause]` | Aalen-Johansen cumulative incidence of each cause, Gray's test across groups, number at risk |
| `calibration` | binary outcome, predicted probabilities | `--y --x "p1,p2" [--positive]` | Brier score, calibration intercept and slope, Hosmer-Lemeshow test, decile plot with a lowess smooth |
| `sankey` | patient, visit, state (one row per patient visit) | `--id --x --y [--levels --visit_levels]` | state transitions across visits; Bowker's or McNemar's test for two visits |
| `ma_plot` | DE table: log2FC, p, mean expression (`baseMean`, `AveExpr` or `logCPM`) | `--data [--x MEAN --y P --label]` | M against A, up and down genes coloured, top hits labeled (the test is upstream) |
| `manhattan` | one row per SNP: chromosome, position, p | `--group --x --y [--label]` | Manhattan and QQ plot, 5e-8 and 1e-5 lines, lead SNPs, lambda_GC |
| `feature_panel` | group, value, feature (long) or one column per feature (wide) | `--x [--y --facet --levels]` | one panel per feature; Welch t or Mann-Whitney, ANOVA or Kruskal-Wallis; BH across the features |
| `embedding` | UMAP or t-SNE coordinates and a cluster, computed upstream | `--x --y --group` | scatter by cluster, labels at the medians, silhouette width per cluster (no test) |

Common options: `--engine r|python` (default r), `--geom box|violin|bar|raincloud|interaction` (interaction: `factorial_anova` only),
`--paired`, `--test` (`two_group_compare`: `auto|welch_t|student_t|wilcoxon`;
`multi_group_compare`: `auto|anova|welch|kruskal`; `paired_compare`: `auto|t|wilcoxon`;
`factorial_anova`: `auto|anova|art`; `feature_panel`: `auto|welch_t|wilcoxon|anova|kruskal`;
`estimation_plot`: `auto|welch_t`; `correlation`: `auto|pearson|spearman|kendall`;
`correlation_heatmap`: `pearson|spearman|kendall`; `proportions`: `auto|chisq|fisher`;
`timecourse`: `auto|gg|uncorrected`; `standard_curve`: `auto|linear|4pl`; `heatmap`:
`auto|welch_t`; an unknown value stops the run), `--levels "a,b,c"` (the order of the
`--x` values; the states of `sankey`), `--visit_levels "a,b,c"` (the `sankey` visits),
`--features "a,b"` (`pca`), `--scale`/`--cluster` (heatmaps), `--fc_cutoff`/`--p_cutoff`/`--top_n`
(volcano, ma_plot; `--top_n` also sets the terms of enrichment_dot and the intersections
of upset), `--id --group --label`, `--positive LEVEL` (roc_curve case level, calibration event level),
`--cause CODE` (cumulative_incidence cause of interest, default 1),
`--ci_low --ci_high --se --measure ratio|difference` (forest_meta),
`--pd_threshold --pr_threshold` (spider_response, waterfall), `--markers "a,b"`
(swimmer event-time columns),
`--theme pubready_house|prism`, `--xlab --ylab --title --palette`, `--out DIR`.

For `roc_curve`, confirm the case level with the scientist: by default it is the
second level in sorted order, so a `healthy`/`disease` outcome needs
`--positive disease`. The recipe gives no warning when it takes the default.
`bland_altman` reports `--x` minus `--y`, so put the new method in `--x` and the
reference method in `--y`. For `forest_meta`, pass `--measure ratio` for an odds, risk
or hazard ratio when the CIs are not given on the ratio scale.

For `estimation_plot`, confirm the control group with the scientist and name
it with `--ref LEVEL`. By default the control is the first level in sorted order.
The figure shows the effect size and its bootstrap CI, and the stats table gives
the t-test. For `standard_curve`, leave the concentration empty on each unknown
row. The recipe does not extrapolate. It flags each unknown outside the range of
the standards and gives it no concentration.

For `cumulative_incidence`, code the status as 0 for a censored subject and as
one code for each event type. Name the cause of interest with `--cause`. For
`calibration`, give the predicted probabilities between 0 and 1, and confirm the
event level as for `roc_curve`. For `sankey`, give one row for each patient and
visit. A status of `cumulative_incidence` that holds text or a negative code stops
the run.

The order of the levels follows two rules in both engines. The groups of a
comparison sort as numbers when every label is a number, else in code-point
order, so "Drug" comes before "control". The first level is group 1 of the
contrast, and the first level of a text covariate is the `cox_forest` reference.
Paired conditions, the `timecourse` time points and the `sankey` visits follow
the time rule: numeric order, else the one number that each label holds, so
"week 8" comes before "week 16", else the order of first appearance. The
`proportions` groups and the `sankey` states keep the order of their first
appearance. To set an order, pass every value with `--levels`, for example
`--levels "CR,PR,SD,PD"` for the `sankey` states, or `--visit_levels` for the
`sankey` visits. When a value holds a comma, separate the values with "|".
`paired_compare` reads long rows. When the two conditions are in two columns,
such as `pre` and `post`, first make one row for each subject and condition,
with the subject id in its own column.

## Choosing the test

The full test-selection tree is in `reference/decision_tree.md`, and the assumption
checks behind it in `reference/assumptions.md`; walk those with the scientist rather
than deciding from memory. In short, `--test auto` measures before it picks. It
checks normality with Shapiro-Wilk: of each group, of the paired differences, or
of the model residuals. A constant group is not normal. Two independent groups get
Welch's t-test when both are normal, and the Mann-Whitney U test otherwise. For
three groups or more, Levene's test, centered on the median, chooses between the
ANOVA and Welch's ANOVA. When the rank test cannot reach P < 0.05 at the group
sizes, for example 3 against 3, auto keeps the parametric test, and the methods
text says why. `figkit diagnose` takes the options of `figkit plot` and shows that
evidence as a QQ plot, so the scientist can confirm the choice. Override any step
with `--test`. Event columns code `1 = event, 0 = censored`, or TRUE and FALSE.
Any other value, such as "dead" or 2, stops the run. `survival_km` needs `--x` for
the grouping arm.

## Designs that no recipe fits

Some designs need a model that pubready does not have. When the design is one of
these, tell the scientist, and do not force the nearest recipe. A pooled test on
such data can give the wrong answer. For example, the pooled admissions of UC
Berkeley give an odds ratio of 1.84 against women, but within the departments the
Mantel-Haenszel odds ratio is 0.90.

- A third variable splits the data into strata, such as a department or a
  center. The analysis needs a Mantel-Haenszel test or a regression model.
- The comparison must adjust for a covariate, such as age. The analysis needs
  ANCOVA or a regression model.
- The question is about a rate of change, or the repeated measures have
  missing visits. The analysis needs a mixed model.
- Each group is compared with one control only. `multi_group_compare` compares
  all pairs, and pubready has no Dunnett test.

## Engine limits

Pick the engine per figure with `--engine`. Both write the same bundle and, where
a matching method exists, the same numbers. The two-group and paired recipes agree
digit for digit on the t-test, the Mann-Whitney U test and the Wilcoxon signed-rank
test, and the smoke test asserts it. The public benchmarks in `benchmarks/` hold
32 of the 34 recipes in both engines to an independent reference value. No
public dataset in the image fits `spider_response` or `waterfall`, so only the
smoke test checks those two. Five gaps remain.
The aligned rank transform for non-parametric factorial ANOVA is R only, so
`--test art --engine python` stops and routes you to R. Under `--test auto`, when
the model residuals fail Shapiro-Wilk, R runs the aligned rank transform and
Python stops with the same route, so the two engines never report two tests.
The Python `survival_km` table has no confidence interval for the median survival.
In `multi_group_compare` a Tukey or Games-Howell P value below about 1e-5 can
differ between the engines in the fourth digit or later, because R's `ptukey`
loses relative precision in the far tail. Both engines write a P value below 1e-12
as "< 1e-12". The `pca`
PERMANOVA is exact, and identical in both engines, when the design has 5000
labelings or fewer. A larger design uses 999 seeded permutations, and the two
engines draw different permutations, so the two P values differ by the Monte Carlo
error.
For a table larger than 2 x 2, the Python Fisher's exact test enumerates every
table with the observed margins. Above one million tables it stops and routes you
to `--test chisq` or to the R engine, which uses the network algorithm.

## Read more

- `docs/how_it_works.md` — the full option table, an example for every recipe, and the bundle contents.
- `example/README.md` — worked commands and the committed reference figures.
  A uv install does not hold this file. Read it in a clone or on GitHub.
- `reference/` — the statistical decision tree and assumption-check notes.

## Install

Codex reads the nearest `AGENTS.md` from the working directory up. To activate
this adapter for the repo, link it at the repo root:

```bash
ln -s skills/codex/AGENTS.md AGENTS.md
```

If the repo already has a root `AGENTS.md`, append the contents of this file to it
instead of overwriting.
