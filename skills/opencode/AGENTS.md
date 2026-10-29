# pubready (opencode adapter)

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
figkit inspect --data FILE [--format json]    # columns, types, group sizes, candidate x/y
figkit diagnose --recipe NAME --data FILE ... # QQ + assumption check, recommend a test (no figure)
figkit plot --recipe NAME --data FILE ...     # draw + test, writes a bundle
figkit render --config plot_config.yaml       # redraw from an edited spec, recompute stats
```

If a command reports the image is missing, run `figkit build` first.

## Workflow

1. Run `figkit inspect` on the data. Read the real column types and group sizes.
   Do not guess them.
2. Talk through the design with the scientist: how many groups, paired or not, one
   factor or several, time to event, or a matrix. If they are unsure which test
   applies, walk the plain-language questions in `reference/decision_tree.md`.
   Confirm the design before drawing.
3. If no test is named, run `figkit diagnose --recipe … --data … [--x --y …]`. It
   draws a QQ plot with the Shapiro-Wilk p, skewness, and kurtosis for the quantity
   the test depends on, and prints a recommended test. Show the scientist the
   recommendation and the QQ, and let them accept it or force another with `--test`.
   Shapiro is a default, not a verdict (underpowered at small n, over-rejects at
   large n), so the QQ plot is the tiebreaker.
4. Call `figkit plot` with the recipe and columns. Leave `--test auto` unless a test
   is forced. The bundle also includes the QC panel (`qc_normality_<ts>.png`).
5. Point them at the bundle folder under `pubready_output/`: the figure, the
   standalone script, the stats table, the input copy, and a methods paragraph.
6. To restyle, edit `plot_config.yaml` in the bundle and run `figkit render`. It
   recomputes the statistics, so the figure and the numbers never drift apart.

## Recipes

| Recipe | Tidy input | Columns you pass | Test |
|---|---|---|---|
| `two_group_compare` | group, value | `--x --y` | t-test or Mann-Whitney, paired or unpaired |
| `multi_group_compare` | group, value | `--x --y` | one-way ANOVA / Welch / Kruskal-Wallis + post hoc |
| `factorial_anova` | f1, f2[, f3], value | `--x --y --fill [--facet]` | two-way / three-way ANOVA (R also does the aligned rank transform) |
| `survival_km` | time, event, group | `--time --event --x` | Kaplan-Meier + log-rank |
| `cox_forest` | time, event, covariates | `--time --event --covariates "a,b,c"` | Cox proportional hazards |
| `heatmap` | matrix CSV + annotation CSV | `--data --annotation` | per-feature Welch t-test when the annotation has two groups |
| `correlation` | x, y (two numeric columns) | `--x --y` | Pearson / Spearman / Kendall, chosen from normality |
| `correlation_heatmap` | table of numeric variables | `--data` | pairwise correlation, clustered, BH-adjusted stars |
| `upset` | binary membership matrix | `--data [--top_n]` | descriptive intersection sizes, the 15 largest drawn by default (no test) |
| `volcano` | DE table (log2FC + p) | `--data [--x --y --label]` | EnhancedVolcano-style scatter, top hits labeled |
| `enrichment_dot` | GSEA / ORA result table | `--data [--label]` | dot plot, size = count, color = -log10 adjusted p |
| `paired_compare` | condition, value, id | `--x --y --id` | paired t-test or Wilcoxon, connecting lines |
| `proportions` | two categoricals | `--x --y` | chi-square / Fisher, 100% stacked bar |
| `pca` | samples × features + group | `--data --group` | PC1/PC2 scatter, 95% ellipses, PERMANOVA |
| `spider_response` | patient, time, value | `--x --y --id [--group]` | descriptive change from baseline, RECIST lines (no test) |
| `waterfall` | patient, time, value (or one best change per patient) | `--y --id [--x --group]` | descriptive best change from baseline, sorted bars, RECIST lines (no test) |
| `dose_response` | dose, response[, curve] | `--x --y [--group]` | four-parameter logistic fit, EC50 or IC50 with 95% CI, shared-EC50 F test |
| `roc_curve` | binary outcome, numeric markers | `--y --x "m1,m2" [--positive]` | ROC curves, AUC with DeLong CI, Youden threshold, DeLong test between markers |
| `bland_altman` | method A, method B (one row per subject) | `--x --y` | bias and 95% limits of agreement with CIs, proportional bias test |
| `forest_meta` | one row per study: estimate and 95% CI or SE | `--y --ci_low --ci_high [--se --label --measure]` | fixed-effect and DerSimonian-Laird pooling, Q, I2, tau2 |
| `timecourse` | subject, time, value[, group] | `--x --y --id [--group]` | mixed (or repeated-measures) ANOVA, Greenhouse-Geisser when Mauchly rejects |
| `swimmer` | one row per patient: duration[, status, arm, event times] | `--id --time [--event --group --markers "a,b"]` | descriptive bar per patient, arrows for ongoing courses (no test) |
| `distribution` | one numeric variable[, group] | `--y [--x]` | descriptive histogram, kernel density and normal fit per group, Shapiro-Wilk (no comparison) |
| `composition` | sample, category[, count, condition] | `--id --x [--y --group]` | 100% stacked bar per sample, per-category Mann-Whitney or Kruskal-Wallis, BH-adjusted (exploratory) |
| `estimation_plot` | group, value[, subject] | `--x --y [--ref --paired --id]` | mean difference against the reference with a 95% BCa bootstrap CI, Welch or paired t-test; Gardner-Altman or Cumming plot |
| `enzyme_kinetics` | substrate, rate[, condition] | `--x --y [--group]` | Michaelis-Menten fit, Vmax and Km with 95% CIs, F tests of one shared Km and one shared Vmax |
| `standard_curve` | concentration (empty for an unknown), signal[, sample] | `--x --y [--id --test auto\|linear\|4pl]` | 4PL or linear standard curve, unknowns read off the curve, no extrapolation, recovery of the standards |
| `cumulative_incidence` | time, status (0 = censored, other codes = causes)[, group] | `--time --event [--x --cause]` | Aalen-Johansen cumulative incidence of each cause, Gray's test across groups, number at risk |
| `calibration` | binary outcome, predicted probabilities | `--y --x "p1,p2" [--positive]` | Brier score, calibration intercept and slope, Hosmer-Lemeshow test, decile plot with a lowess smooth |
| `sankey` | patient, visit, state (one row per patient visit) | `--id --x --y` | state transitions across visits; Bowker's or McNemar's test for two visits |
| `ma_plot` | DE table: log2FC, p, mean expression (`baseMean`, `AveExpr` or `logCPM`) | `--data [--x MEAN --y P --label]` | M against A, up and down genes coloured, top hits labeled (the test is upstream) |
| `manhattan` | one row per SNP: chromosome, position, p | `--group --x --y [--label]` | Manhattan and QQ plot, 5e-8 and 1e-5 lines, lead SNPs, lambda_GC |
| `feature_panel` | group, value, feature (long) or one column per feature (wide) | `--x [--y --facet]` | one panel per feature; Welch t or Mann-Whitney, ANOVA or Kruskal-Wallis; BH across the features |
| `embedding` | UMAP or t-SNE coordinates and a cluster, computed upstream | `--x --y --group` | scatter by cluster, labels at the medians, silhouette width per cluster (no test) |

Common options: `--engine r|python` (default r), `--geom box|violin|bar|raincloud|interaction` (interaction: `factorial_anova` only),
`--paired`, `--test auto|welch_t|wilcoxon|anova|kruskal|art` (correlation:
`pearson|spearman|kendall`; proportions: `chisq|fisher`; timecourse:
`gg|uncorrected`), `--scale`/`--cluster` (heatmaps), `--fc_cutoff`/`--p_cutoff`/`--top_n`
(volcano, ma_plot; `--top_n` also sets the terms of enrichment_dot and the intersections
of upset), `--id --group --label`, `--positive LEVEL` (roc_curve case level, calibration event level),
`--cause CODE` (cumulative_incidence cause of interest, default 1),
`--ci_low --ci_high --se --measure ratio|difference` (forest_meta),
`--pd_threshold --pr_threshold` (spider_response, waterfall), `--markers "a,b"`
(swimmer event-time columns),
`--theme pubready_house|prism`, `--xlab --ylab --title --palette`, `--out DIR`.

For `roc_curve`, confirm the case level with the scientist: by default it is the
second level in sorted order, so a `healthy`/`disease` outcome needs
`--positive disease`. For `forest_meta`, pass `--measure ratio` for an odds, risk
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
visit. Numeric states sort as numbers, and text states keep the order of their
first appearance in the data.

## Choosing the test

The full test-selection tree is in `reference/decision_tree.md`, and the assumption
checks behind it in `reference/assumptions.md`; walk those with the scientist rather
than deciding from memory. In short, `--test auto` measures before it picks: it checks
normality (Shapiro-Wilk) and, for group comparisons, equal variance (F test or
Levene), then goes parametric when the assumptions hold and non-parametric otherwise,
with the Welch correction when variances differ. `figkit diagnose` shows that evidence
as a QQ plot so the scientist can confirm the choice. Override any step with `--test`.
Event columns code `1 = event, 0 = censored`; `survival_km` needs `--x` for the
grouping arm.

## Engine limits

Pick the engine per figure with `--engine`. Both write the same bundle and, where
a matching method exists, the same numbers. The two-group and paired recipes agree
digit for digit on the t-test, the Mann-Whitney U test and the Wilcoxon signed-rank
test, and the smoke test asserts it. The public benchmarks in `benchmarks/` hold
every recipe in both engines to an independent reference value. Four gaps remain.
The aligned rank transform for non-parametric factorial ANOVA is R only, so
`--test art --engine python` errors and routes you to R.
The Python `survival_km` table has no confidence interval for the median survival.
In `multi_group_compare` a Tukey P value below about 1e-8 differs between the
engines in the second or third digit, because R's `ptukey` loses relative precision
in the far tail; both values stay far below any significance threshold. The `pca`
PERMANOVA is exact, and identical in both engines, when the design has 5000
labelings or fewer. A larger design uses 999 seeded permutations, and the two
engines draw different permutations, so the two P values differ by the Monte Carlo
error.

## Read more

- `docs/how_it_works.md` — the full option table, an example for every recipe, and the bundle contents.
- `example/README.md` — worked commands and the committed reference figures.
- `reference/` — the statistical decision tree and assumption-check notes.
