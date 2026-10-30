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
2. `figkit inspect` reports the columns, their types, group sizes, candidate
   x and y columns, and the columns that look like a subject id. The data shape is
   measured, not guessed.
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
figkit doctor         # check Podman, the image and the mounts, and print each fix
figkit shell          # open a shell in the image
figkit test           # run both engines' test suites in the image
figkit version        # print version, image tag/digest, podman, git commit
figkit where          # print the folder that holds core/, the adapters, reference/ and docs/
```

A clone runs `./cli/figkit` from the repository root. A uv install puts a
`figkit` command on the PATH. That command is a small launcher in
`src/pubready_cli/`, and it runs the same `cli/figkit` from the files that the
package carries. A uv install is not a Git checkout, so `cli/figkit` cannot ask
Git for the commit. The launcher reads the commit that uv recorded in the
`direct_url.json` file of the install, and it passes the commit in
`PUBREADY_SOURCE_COMMIT`. Each bundle then records that commit and a dirty state
of false. An install from a local folder has no recorded commit, so `figkit plot`
and `figkit render` stop.

Common options:

| Option | Meaning |
|---|---|
| `--engine r\|python` | which engine draws the figure (default r) |
| `--x`, `--y` | the axis columns for tidy data |
| `--fill`, `--facet` | second and third factors (factorial) |
| `--time`, `--event` | survival columns (event coded 1 = event, 0 = censored, or TRUE and FALSE; any other value stops the run) |
| `--covariates "a,b,c"` | Cox model covariates |
| `--annotation FILE` | heatmap column annotation CSV |
| `--test` | force a test instead of `auto`; the values of each recipe are in the table below, and an unknown value stops the run |
| `--paired` | paired design |
| `--geom box\|violin\|bar\|raincloud\|interaction` | figure geometry; `interaction` is the interaction plot of `factorial_anova` |
| `--scale`, `--cluster` | heatmap scaling and clustering |
| `--id`, `--group`, `--label` | the subject, the grouping and the label columns of the recipes that use them |
| `--features "a,b"` | the feature columns of `pca`; without it, every numeric column other than the group, `--id` and `--label` is a feature |
| `--pd_threshold`, `--pr_threshold` | the two RECIST reference lines of `spider_response` and `waterfall`, +20 and -30 by default |
| `--positive LEVEL` | the case level of the `roc_curve` outcome and the event level of the `calibration` outcome (default: the second level in sorted order) |
| `--ci_low`, `--ci_high`, `--se` | the 95% CI bounds, or the standard error, of each study in `forest_meta` |
| `--measure auto\|ratio\|difference` | the pooling scale of `forest_meta` (a ratio is pooled on the log scale) |
| `--markers "a,b"` | the columns of event times that `swimmer` draws on each bar |
| `--top_n N` | the number of labels in `volcano`, of terms in `enrichment_dot` and of intersections in `upset`, 15 by default. The `upset` title states how many intersections it hides. |
| `--xlab`, `--ylab`, `--title`, `--palette` | labels and colors |
| `--out DIR` | output root (default `pubready_output`) |
| `--llm-provider`, `--llm-model` | the provider and the model of the LLM that chose the figure, recorded in the report (default `none`) |
| `--agent-harness`, `--agent-harness-version` | the agent harness and its version, recorded in the report (default `manual` and `not_recorded`) |
| `--ref LEVEL` | the reference group of `estimation_plot` (default: the first level in sorted order) |
| `--test auto\|linear\|4pl` | the curve of `standard_curve`: auto fits a 4PL when the standards have 5 or more positive concentrations and the fit converges, and a line in all other conditions |
| `--cause CODE` | the event code of the cause of interest of `cumulative_incidence` (default 1) |
| `--levels "a,b,c"` | the order of the `--x` values: the groups of the comparison recipes and of `proportions`, the conditions of `paired_compare` and the time points of `timecourse`; in `sankey`, the order of the states. The list must name every value of the data. When a value holds a comma, separate the values with "\|" |
| `--visit_levels "a,b,c"` | the order of the `sankey` visits; the list must name every visit of the data |
| `--fc_cutoff`, `--p_cutoff`, `--top_n` | the fold-change cutoff, the p cutoff and the count of labels of `volcano` and `ma_plot` (default 1, 0.05 and 15); `--top_n` also sets the most lead SNPs that `manhattan` labels |
| `--facet COL` (`feature_panel`) | the feature column of long data; without it, every numeric column other than `--x` is a feature |

Write a value that starts with "-" in the `--option=value` form, for example
`"--ylab=-log10 adjusted p"`. The R option parser reads `--ylab "-log10"` as
two options and stops.

The `--test` values of each recipe:

| Recipe | `--test` values |
|---|---|
| `two_group_compare` | `auto`, `welch_t` (also `t`), `student_t`, `wilcoxon` (also `mann_whitney`) |
| `multi_group_compare` | `auto`, `anova`, `welch`, `kruskal` |
| `paired_compare` | `auto`, `t`, `wilcoxon` (also `signed_rank`) |
| `factorial_anova` | `auto`, `anova`, `art` (R only) |
| `feature_panel` | `auto`, `welch_t`, `wilcoxon`, `anova`, `kruskal` |
| `estimation_plot` | `auto`, `welch_t` |
| `correlation` | `auto`, `pearson`, `spearman`, `kendall` |
| `correlation_heatmap` | `auto` (Pearson), `pearson`, `spearman`, `kendall` |
| `proportions` | `auto`, `chisq`, `fisher` |
| `timecourse` | `auto`, `gg`, `uncorrected` |
| `standard_curve` | `auto`, `linear`, `4pl` |
| `heatmap` | `auto` (no test), `welch_t` |

### The options of each recipe

Each recipe reads its own options, and `core/recipe_options.txt` lists them. Both
engines read that one file. `figkit plot` and `figkit diagnose` stop when you give
an option that the recipe does not read, before any output, for example
"--title is not an option of survival_km; it is read by: ...". An option can
therefore never run with no effect. A new recipe needs a line in the file.

### The read rule

Both engines read a CSV file with one rule. Each cell is read as text, and the
spaces and tabs around each value and each column name are removed. Only an empty
cell or "NA" is missing, so a response named "None" stays a response. A column is
numeric when each present cell is a number with "." as the decimal mark, and
logical when each present cell is TRUE or FALSE. A cell that is a number with a
comma, such as "85,45" or "4,380", stops the run and names the column, because
the comma can be a decimal mark or a thousands mark. An Excel cell keeps the type
that the sheet gives it, and its text cells get the same rule. The data log has a
"Read rule" section that counts each value that the read changed.

### What `inspect` reports

`figkit inspect` prints the same words in both engines. For each column it gives
the type (numeric, text, logical, date or other), the missing count, the count of
unique values, and each level in quotes with its count of rows. A numeric column
with 2 to 6 distinct values is a candidate group, marked as a numeric code. A
column whose values repeat across the levels of another column, at most once in
each level, is a possible subject id. With `--format json`, the output adds
`level_counts`, `group_candidates`, `outcome_candidates` and `id_candidates`.

### Render

`figkit render --config plot_config.yaml` looks for the bundle input in the
folder of the config first and then in the current folder, so it runs from the
root of the repository. `--data` names another file. Without `--engine`, render
uses the engine that the config names, and the manifest records the engine that
drew the figure. The Python engine prints a data error as one line,
`Error: <message>`, with no traceback.

### Check the setup

`figkit doctor` checks the host setup and prints one line for each check. A
failed check also prints its fix. The command never installs a package and never
runs `sudo`. The user runs each fix that needs administrator rights. The exit
code is 1 when a check fails. The exit code is 0 when all checks pass or when a
check gives only a warning.

| Check | Result | Fix that doctor prints |
|---|---|---|
| Platform | Fails in Git Bash, MSYS or Cygwin on Windows. | Install WSL2 and run figkit in the Ubuntu shell. |
| bash | Fails before bash 3.2. `cli/figkit` runs on the bash 3.2 of macOS. | Install a newer bash. |
| podman | Fails when `podman` is not on the PATH. | The install command of the platform. |
| Podman machine | On macOS, fails when no machine exists or when the machine is stopped. | `podman machine init` and `podman machine start`. |
| Podman service | Fails when `podman info` fails. Doctor shows the first line of the error. | Read the full error of `podman info`. |
| Rootless mode | Warns when Podman runs rootful. figkit works, but the files that it writes can belong to root. | Run figkit as your own user. On macOS, set the machine back to rootless. |
| Subordinate IDs | On Linux, fails when rootless Podman has no subordinate UID range. Podman then cannot unpack the image. | `sudo usermod --add-subuids 100000-165535 --add-subgids 100000-165535 USER`, then `podman system migrate`. |
| Folders | On macOS, warns when the current folder or the pubready files are outside the home folder. On WSL2, warns when the current folder is on a Windows drive. | Move the files. |
| Image | Fails when the image is not built. | `figkit build`. |
| Test container | Fails when a container with the mounts of `figkit plot` cannot list the two mounted folders. An agent sandbox that blocks user namespaces gives this failure. | Run figkit outside the sandbox. |

The install command depends on the platform. Doctor uses `apt` on Ubuntu and
Debian, `dnf` on Fedora, RHEL and CentOS, `pacman` on Arch and `zypper` on
openSUSE. On macOS it uses `brew install podman`, then `podman machine init` and
`podman machine start`. On other platforms it points to the Podman install guide.

The other commands check Podman before they use it. If `podman` is not on the
PATH, they print the install command of the platform. If Podman cannot reach its
service, for example when the macOS machine is stopped, they tell the user to run
`figkit doctor`. Only a missing image gives the hint `figkit build`.

Each `podman run` of figkit has the option `--security-opt label=disable`. On a
host with SELinux in enforcing mode, such as Fedora, the container can then read
the mounted folders, and Podman does not relabel your files. The option does
nothing on a host without SELinux. The `REPRODUCE.md` command of each bundle has
the same option.

## The recipes

Each recipe is a named function with a fixed argument surface, present in both
engines. Pick the engine per figure with `--engine`.

| Recipe | Tidy input | Test | Figure |
|---|---|---|---|
| `two_group_compare` | `group, value` | Welch t-test or Mann-Whitney; with `--paired`, the paired t-test or the Wilcoxon signed-rank test from the normality of the differences | box or violin with a bracket |
| `multi_group_compare` | `group, value` | one-way ANOVA / Welch / Kruskal-Wallis + post hoc | box with brackets for significant pairs, or a compact letter display when more than six pairs differ; omnibus as a subtitle |
| `factorial_anova` | `f1, f2[, f3], value` | two-way and three-way ANOVA (R also does the aligned rank transform) | grouped box or bar, faceted, effects on top |
| `survival_km` | `time, event, group` | Kaplan-Meier with a log-rank test | step curves, log-rank p, number-at-risk table (time >= t), legend right of the panel |
| `cox_forest` | `time, event, covariate…` | Cox proportional hazards; a term that separates the events is marked "not estimable, separation" in the `note` column | forest plot of hazard ratios with 95% CI |
| `heatmap` | matrix CSV + annotation CSV | a Welch t-test of each feature only with `--test welch_t` and two annotation groups, BH-adjusted | clustered heatmap, z-scored, group annotation bar, stars with the test |
| `correlation` | `x, y` (two numeric columns) | Pearson / Spearman / Kendall, chosen from normality | scatter with method, coefficient and P in the subtitle; a Pearson scatter also has the least-squares line and its 95% CI band, and a rank correlation has no line |
| `correlation_heatmap` | table of numeric variables | pairwise correlation, BH-adjusted | clustered correlation matrix, coefficients + significance stars |
| `upset` | membership matrix of 1/0, TRUE/FALSE, T/F, yes/no or y/n; an empty cell is not a member, and any other value stops the run | descriptive (intersection sizes, no test) | UpSet plot of the largest set intersections, 15 by default |
| `volcano` | DE table (log2FC + p) | thresholds on fold change and p (test is upstream) | EnhancedVolcano-style four-color scatter, top hits labeled |
| `enrichment_dot` | GSEA / ORA result table | descriptive (enrichment is upstream); a table with only a raw p is marked as not adjusted | dot plot, size = count, color = -log10 adjusted p; an adjusted p of 0 is a triangle |
| `paired_compare` | `condition, value, id` | paired t-test or Wilcoxon signed-rank | before/after points with connecting lines and a bracket |
| `proportions` | two categoricals | chi-square, or Fisher when a cell is sparse | 100% stacked bar of proportions with the p |
| `pca` | samples × features + group; `--features` names the features, and the `--id` and `--label` columns are never features | PERMANOVA for group separation | PC1/PC2 scatter, 95% ellipses, variance on the axes |
| `spider_response` | `patient, time, value` | descriptive (change from baseline, no test) | one line per patient, RECIST reference lines at +20% and -30% |
| `waterfall` | `patient, time, value`, or one best change per patient | descriptive (best change from baseline, no test) | one sorted bar per patient, RECIST reference lines |
| `dose_response` | `dose, response[, curve]` | four-parameter logistic fit; extra-sum-of-squares F test of one shared EC50 over the determined curves | fitted curves on a log dose axis, mean +/- SEM at each dose, EC50 or IC50 in the legend; a curve that is not determined is drawn as points only |
| `roc_curve` | binary outcome, numeric markers | AUC with a DeLong 95% CI; DeLong test of each marker against the first | ROC curves with the Youden point, AUCs in the legend |
| `bland_altman` | `method_a, method_b`, one row per subject | bias and 95% limits of agreement with CIs; regression test of proportional bias | difference against mean, shaded CIs of the bias and the limits |
| `forest_meta` | one row per study: estimate and 95% CI or SE | fixed-effect and DerSimonian-Laird pooling; Cochran's Q, I2, tau2 | forest plot with weighted squares, pooled diamonds and an estimate column |
| `swimmer` | one row per patient: `duration[, status, arm, event times]` | descriptive (no test) | one bar per patient sorted by duration, arrows for ongoing courses, symbols for events |
| `timecourse` | `subject, time, value[, group]` | mixed ANOVA (group between, time within) or repeated-measures ANOVA; Greenhouse-Geisser when Mauchly's test rejects sphericity | mean +/- SEM per group over time, one thin line per subject |
| `distribution` | `value[, group]` | descriptive: moments, quantiles and the Shapiro-Wilk test of each group, no comparison | histogram with Freedman-Diaconis bins, kernel density and normal fit, one panel per group |
| `composition` | `sample, category[, count, condition]` | Mann-Whitney U or Kruskal-Wallis test of each category on the per-sample proportions, BH-adjusted | one 100% stacked bar per sample, one strip per condition |
| `estimation_plot` | `group, value[, subject]` | mean difference against the reference with a 95% BCa bootstrap CI from 5000 Park-Miller resamples, the same in both engines; Welch or paired t-test, Hedges' g | Gardner-Altman plot for two groups, Cumming plot for three or more; a slopegraph with `--paired` |
| `enzyme_kinetics` | `substrate, rate[, condition]` | Michaelis-Menten fit for each curve, Vmax and Km with 95% CIs; extra-sum-of-squares F tests of one shared Km and one shared Vmax | mean +/- SEM at each concentration, the fitted curves on linear axes from 0, Vmax and Km in the legend |
| `standard_curve` | `concentration, signal[, sample]`, an empty concentration for an unknown | 4PL or linear fit to the standards; each unknown read off the inverse curve, with no extrapolation; recovery of each standard level | standards as mean +/- SD, the curve, the unknowns as open diamonds on the curve with drop lines |
| `cumulative_incidence` | `time, status[, group]`, status 0 = censored, other codes = causes | Aalen-Johansen cumulative incidence of each cause; Gray's test across groups for each cause | one panel per cause, step curves per group, Gray's P in the panel title, number-at-risk table |
| `calibration` | binary outcome, one column of predicted probabilities per model | Brier score, calibration intercept and slope with Wald CIs, Hosmer-Lemeshow test, observed/expected ratio | decile points with Wilson CIs, lowess smooth, ideal line, spike histogram of the predictions |
| `sankey` | `patient, visit, state`, one row per patient visit | transitions table; Bowker's test of symmetry, or McNemar's test, when there are two visits | stacked nodes per visit, ribbons coloured by the state of origin, labels with counts |
| `ma_plot` | DE table: log2FC, p and mean expression (`baseMean`, `AveExpr` or `logCPM`) | thresholds on fold change and p (the test is upstream) | M against A, up genes red, down genes blue, top hits labeled, counts in the subtitle |
| `manhattan` | one row per SNP: `chr, pos, p[, snp]` | descriptive; genomic inflation factor lambda_GC | Manhattan plot with the 5e-8 and 1e-5 lines and the lead SNPs, QQ plot with a 95% band |
| `feature_panel` | `group, value, feature` (long) or one column per feature (wide) | Welch t-test or Mann-Whitney U (2 groups), one-way ANOVA or Kruskal-Wallis (3 or more); one family for the panel; BH across the features | one panel per feature, box and points, the adjusted p in each title |
| `embedding` | `UMAP_1, UMAP_2, cluster`, computed upstream | descriptive; mean silhouette width of each cluster on the 2D coordinates | scatter coloured by cluster, labels at the cluster medians, equal aspect ratio |

### Test selection

`--test auto` walks a short decision tree rather than picking a test from a single
prompt. It checks normality with the Shapiro-Wilk test: of each group, of the
paired differences, or of the model residuals. A constant sample is not normal.
Two independent groups get Welch's t-test when both are normal, and the
Mann-Whitney U test otherwise. Levene's test runs too, and its P goes to the stats
table. For three groups or more, Levene's test, centered on the median, chooses
between the ANOVA and Welch's ANOVA. The `correlation` recipe uses the same
normality check to pick Pearson (both variables normal) or Spearman. The scientist
can override any step by passing `--test` explicitly.

At a small n a rank test cannot reach P < 0.05: 3 against 3 gives at best P = 0.10.
Before `--test auto` takes a rank test, the recipe runs the rank test on data with
the same group sizes, no tie and complete separation. When that smallest P is above
0.05, the parametric test stays, `test_meta` records `rank_test_attainable_p`, and
the methods text says why. The rule applies to `two_group_compare`,
`paired_compare`, `multi_group_compare` and `feature_panel`.

A test that divides by a zero variance stops and names the groups, and an outcome
with zero variance stops before any test. Cohen's d uses the pooled SD, and the
paired effect is Cohen's dz.

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
and over-rejects at large n, so read the plot and tolerate mild deviation. Each
`figkit plot` bundle of a recipe that chooses its test from normality carries the
same panel as `qc_normality_*`. `figkit diagnose` takes every option of
`figkit plot`.

### The bench assay recipes

`estimation_plot` gives each group against the reference as a mean difference
with a 95% bias-corrected and accelerated (BCa) bootstrap interval. The two
engines draw the same 5000 resamples, because neither engine uses its own random
number generator. Both engines run the Park-Miller minimal standard generator
with the multiplier 48271, the modulus 2^31 - 1 and the seed 20261029. This
calculation is exact in double precision, so the two intervals agree to the last
digits.

Each comparison starts again from the seed. A resample draws the reference
indices first and then the group indices, as floor(u n). A paired design draws
whole pairs. The bias correction is the normal quantile of the proportion of
bootstrap differences strictly below the estimate. The acceleration comes from
the jackknife of the mean difference, and the interval ends are type-7
quantiles.

`enzyme_kinetics` fits v = Vmax S / (Km + S) to each curve. Each 95% CI is the
estimate plus or minus t(0.975, n - 2) times the asymptotic standard error, as
GraphPad Prism reports it. With two or more curves, two extra-sum-of-squares F
tests compare one shared Km and one shared Vmax with the separate fits.

`standard_curve` reads each unknown off the inverse of the fitted curve. It fits
the four-parameter logistic of `dose_response` when the standards have 5 or more
positive concentrations and the fit converges. In all other conditions it fits a
straight line. The signal range is the range of the curve between the lowest and
the highest standard. An unknown outside it gets the flag `below range` or
`above range` and no concentration, so the recipe never extrapolates. The
recovery of each standard level, read back through the curve, is the QC of the
fit.

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

- The visits follow the time rule: `--visit_levels` when it is given, else
  numeric order when every visit is a number, else the order of the one number
  that each visit holds, so "week 8" comes before "week 16", else the order of
  first appearance. The data log records the order and the rule. `timecourse`
  orders its time points and `paired_compare` its two conditions with the same
  rule, and `--levels` sets their order.
- `--levels "a,b,c"` sets the state order. The list must name every state of
  the data, else the recipe stops and names the missing states. A listed state
  that no row holds is left out of the diagram.
- Without `--levels`, states sort as numbers when every state is a number.
  Other states keep the order of their first appearance in the data.
- A patient who is missing at a visit takes no part in the transitions into or
  out of that visit.
- With exactly two visits that carry the same set of states, Bowker's test of
  symmetry runs on the patients present at both visits. With two states, this
  test is McNemar's test. Neither engine applies a continuity correction.
- A pair of cells that holds no patient adds nothing to the statistic and no
  degree of freedom. `stats::mcnemar.test` returns NaN for such a table, so both
  engines then compute the statistic from its formula.

### Shared figure rules

Both engines apply the same rules to the layout and to the reported numbers.

| Rule | Recipes | What it does |
|---|---|---|
| Brackets or letters | `multi_group_compare` | Up to six significant pairs get one bracket each, on its own level. More than six pairs give a compact letter display. Groups that share a letter do not differ. `test_meta` records the display. |
| Leaf order | `heatmap`, `correlation_heatmap` | Both engines build a complete-linkage tree on the Euclidean distance of the drawn values. At each merge, the branch with the higher mean of the unscaled values comes first. For `correlation_heatmap`, the mean is the mean correlation of each variable. On a tie, the branch with the earlier input row or column comes first. `test_meta` records `row_order` and `column_order`, and the smoke test checks that the two engines agree. The emitted scripts carry the same rule. |
| Not determined | `dose_response` | Two checks before the fit find a flat curve. The first finds a curve whose mean responses vary by less than 10% of the range of the mean responses of all curves. The second finds a curve whose responses do not differ between the doses in a one-way ANOVA (P >= 0.05). The second check also covers a single curve, and it needs a replicate at one dose or more. A flat curve gets no fit. The checks use no fitted value, so both engines give the same reason for a flat curve. A curve is also not determined when its fit fails, when its fitted span is below 10% of the range of the mean responses, when its EC50 is outside the tested doses, or when its 95% CI is not finite. Its EC50 is NA, `fit_status` and `status_reason` give the cause, and the shared-midpoint F test leaves it out. |
| Number at risk | `survival_km` | The table counts the subjects whose time is at or after each break. The breaks use the 1-2-5 step that gives at most six intervals, and the x axis runs to the largest time. |
| Heterogeneity | `forest_meta` | The DerSimonian-Laird tau2 is truncated at zero, so I2 is 0% when Q is below its degrees of freedom. |
| Labels | the comparison recipes | A category label wraps at 12 characters. The y axis shows plain numbers with a comma between the thousands, and scientific notation when the largest break is below 1e-4. |
| Level order | the comparison recipes, `survival_km`, `cox_forest`, `cumulative_incidence` | The groups sort as numbers when every label is a number, else in code-point order, so the order does not depend on the locale. `--levels` sets another order. The first level is group 1 of a contrast and the reference level of a Cox covariate. A whole number reads "80", not "80.0", in both engines. |
| Not estimable | `cox_forest` | A covariate that separates the events, so that at each event time the subject with the event holds the highest or the lowest value in the risk set, has no finite hazard ratio. Its HR, CI and P are NA, and the `note` column, the figure and the methods text say "not estimable, separation". |
| Not determined | `enzyme_kinetics` | A curve whose Km is above the highest substrate, whose Km CI holds 0, whose Vmax is not positive, or whose Vmax CI holds 0 is not determined. Its estimates are NA, `fit_status` and `status_reason` give the cause, and the shared F tests leave it out. |
| p = 0 | `volcano`, `enrichment_dot`, `manhattan` | A p of exactly 0 becomes the smallest positive double, and the figure draws it as a triangle at the top. A p outside 0 to 1 stops `volcano` and `manhattan`. |
| Small P | `multi_group_compare` | A Tukey or Games-Howell P below 1e-12 is written "< 1e-12" in the stats table and on the figure. |

The compact letter display uses the insert-and-absorb algorithm of Piepho (2004),
which is also the algorithm of `multcompView::multcompLetters`. Each engine has
its own copy of the steps, `cld_letters()` in `util.R` and in `util.py`, and the
two unit suites check that the copies give the same letters.

### Checks on the data

Each recipe checks its data before the test, in both engines. A check that
fails stops the run with a message that names the column, the group or the row,
or it writes the cause to the data log and the stats table.

| Recipes | Check | Result |
|---|---|---|
| comparison recipes, `correlation`, `timecourse` | An outcome with zero variance, or a variance that is not finite | The run stops. |
| comparison recipes | A forced t-test or Welch ANOVA that divides by a zero variance | The run stops and names the groups. |
| `feature_panel` | A constant feature | The feature is "not tested: constant", and the BH adjustment leaves it out. |
| `factorial_anova` | An empty cell, a factor with one level, or no residual degrees of freedom | The run stops and names the problem. |
| `paired_compare`, `two_group_compare --paired` | Two rows of one subject in one condition, or groups of unequal size under `--paired` | The run stops and names the subjects. |
| `survival_km`, `cox_forest`, `swimmer` | An event value other than 0, 1, TRUE or FALSE | The run stops and lists the values. |
| `cumulative_incidence` | A status that holds text or a negative code | The run stops and lists the codes. |
| `cox_forest` | A covariate with one value, a collinear covariate, or no event | The run stops and names the column. A covariate is numeric when 80% or more of its cells read as numbers, and the data log names each cell that does not. |
| `waterfall`, `spider_response` | Two rows of one patient at one time point | The run stops and names the patients. A patient whose first row is after the first visit gets that row as the baseline, and the data log records it. |
| `swimmer` | A marker time after the end of the bar or below 0 | The data log records it. |
| `timecourse` | Too few subjects for Mauchly's test | Mauchly's test is not available, and auto applies the Greenhouse-Geisser correction. |
| `forest_meta` | A CI that does not hold its estimate, reversed bounds, a CI of zero width, or an SE of 0 or below | The run stops and names the study. A study with no label is "Study <row>". |
| `heatmap` | A sample that the annotation does not list | The run stops and names the samples. A degenerate test gets P = NA and stays out of the BH count. |
| `correlation_heatmap`, `pca` | A numeric column with a text cell, such as "<0.1" | The cell becomes NA when 80% or more of the column reads as numbers, and the data log names each excluded column. |
| `correlation_heatmap` | A pair with fewer than 3 shared complete rows, or a constant in them | The run stops and names the pair. |
| `pca` | A sample with a missing group or a missing feature value, or a feature with zero variance | The sample or the feature is dropped, and the data log records each drop. |
| `upset` | A cell other than 1/0, TRUE/FALSE, T/F, yes/no or y/n | The run stops. An empty cell is not a member, and the data log counts it. |
| `enrichment_dot` | A repeated term name | The second name gets " (2)", and the data log records it. |
| `calibration` | Fewer than 3 distinct predictions in a model | The calibration slope and the Hosmer-Lemeshow test are not estimable. |
| `roc_curve` | A constant marker or complete separation | The DeLong CI is NA, and the `constant_marker` and `complete_separation` columns say why. |
| `bland_altman` | Constant differences | The proportional-bias regression does not run, and `differences_constant` is TRUE. |
| `standard_curve` | Standards whose mean signals are not monotonic | The subtitle and the methods text warn of a hook effect. |
| `dose_response` | A response of the order 1e3 or more, or 1e-3 or less | The fit runs on the response divided by a power of ten, and the plateaus are scaled back. The methods text counts the rows with a dose of zero. |
| `composition` | A category with a total count of zero | The category is dropped, and the data log records it. |
| `distribution` | A group with no spread | The group gets no density curve. |

### Examples

Each command runs from the root of a clone on a file in `example/`.
`example/README.md` shows the figure of each command.

```bash
# two groups, test chosen automatically
figkit plot --recipe two_group_compare --data example/tumor_volume.csv --x group --y volume

# four groups, one-way ANOVA with Tukey post hoc, significant pairs only
figkit plot --recipe multi_group_compare --data example/gene_expression.csv --x genotype --y expression

# two-way ANOVA, grouped box coloured by the second factor
figkit plot --recipe factorial_anova --data example/twoway_response.csv --y response \
  --x genotype --fill treatment

# Kaplan-Meier, drawn with the Python engine
figkit plot --engine python --recipe survival_km --data example/survival_trial.csv \
  --time time --event event --x arm

# Cox hazard ratios
figkit plot --recipe cox_forest --data example/survival_trial.csv --time time --event event \
  --covariates "arm,age,sex,stage"

# clustered heatmap with a sample annotation and a Welch t-test of each gene
figkit plot --recipe heatmap --data example/expression_matrix.csv \
  --annotation example/expression_annotation.csv --test welch_t

# correlation scatter, Pearson or Spearman chosen from the data
figkit plot --recipe correlation --data example/correlation_xy.csv --x gene_a --y gene_b

# clustered correlation heatmap across a panel of numeric variables
figkit plot --recipe correlation_heatmap --data example/correlation_vars.csv

# UpSet plot from a membership matrix
figkit plot --recipe upset --data example/set_membership.csv

# any comparison recipe in the GraphPad Prism style, or as a raincloud or bar
figkit plot --recipe two_group_compare --data example/tumor_volume.csv --x group --y volume --theme prism
figkit plot --recipe two_group_compare --data example/tumor_volume.csv --x group --y volume --geom raincloud

# volcano plot from a DESeq2 / limma results table, on the adjusted p
figkit plot --recipe volcano --data example/de_results.csv --y padj

# GSEA / ORA enrichment dot plot
figkit plot --recipe enrichment_dot --data example/enrichment_results.csv

# paired before/after with connecting lines
figkit plot --recipe paired_compare --data example/paired_response.csv --x condition --y value --id subject

# proportions with a chi-square or Fisher test
figkit plot --recipe proportions --data example/response_by_arm.csv --x arm --y response

# PCA scatter with ellipses and a PERMANOVA p
figkit plot --recipe pca --data example/pca_samples.csv --group group

# four-parameter logistic dose-response curves, one per compound
figkit plot --recipe dose_response --data example/dose_response.csv --x concentration_um \
  --y viability_pct --group compound

# ROC curves of two markers, with the case level stated
figkit plot --recipe roc_curve --data example/biomarker_roc.csv --y status --x marker_a,marker_b \
  --positive disease

# Bland-Altman agreement between two assays
figkit plot --recipe bland_altman --data example/method_agreement.csv --x new_assay --y reference_assay

# meta-analysis forest plot from hazard ratios and their CIs
figkit plot --recipe forest_meta --data example/meta_trials.csv --y hazard_ratio --ci_low ci_low \
  --ci_high ci_high --label trial

# repeated-measures time course with a mixed ANOVA
figkit plot --recipe timecourse --data example/tumor_growth.csv --x day --y volume_mm3 --id mouse \
  --group arm

# spider plot and waterfall of the change from baseline
figkit plot --recipe spider_response --data example/spider_response.csv --x week \
  --y target_lesion_mm --id patient --group arm
figkit plot --recipe waterfall --data example/spider_response.csv --x week \
  --y target_lesion_mm --id patient --group arm

# swimmer plot of the time on treatment, with response and progression marks
figkit plot --recipe swimmer --data example/swimmer_patients.csv --id patient --time months \
  --event ended --group arm --markers response_month,progression_month

# the shape of one variable in each group: histogram, kernel density and normal fit
figkit plot --recipe distribution --data example/distribution.csv --y cell_area_um2 --x condition

# the cell-type composition of each sample, compared between two conditions
figkit plot --recipe composition --data example/composition.csv --id sample --x cell_type \
  --y n_cells --group condition

# the interaction plot of a two-way ANOVA
figkit plot --recipe factorial_anova --data example/twoway_response.csv --y response \
  --x genotype --fill treatment --geom interaction

# estimation plot: each treatment against the vehicle, with a bootstrap CI
figkit plot --recipe estimation_plot --data example/estimation_plot.csv --x treatment \
  --y speed_um_h --ref vehicle

# Michaelis-Menten kinetics, one curve for each inhibitor concentration
figkit plot --recipe enzyme_kinetics --data example/enzyme_kinetics.csv --x substrate_um \
  --y rate --group inhibitor

# ELISA standard curve, the unknowns read off the curve
figkit plot --recipe standard_curve --data example/standard_curve.csv --x conc_ng_ml \
  --y absorbance --id sample

# cumulative incidence of each cause, compared between the arms
figkit plot --recipe cumulative_incidence --data example/cumulative_incidence.csv \
  --time months --event status --x arm --cause 1

# calibration of two risk models in a validation cohort
figkit plot --recipe calibration --data example/calibration.csv --y event_1y \
  --x risk_clinical,risk_extended --positive yes

# transitions of the response category across assessments
figkit plot --recipe sankey --data example/sankey.csv --id patient --x week --y response \
  --levels "complete response,partial response,stable disease,progressive disease"

# MA plot of a DESeq2 table; baseMean is drawn on the log10 scale
figkit plot --recipe ma_plot --data example/de_basemean.csv

# Manhattan and QQ plot of GWAS summary statistics
figkit plot --recipe manhattan --data example/gwas_summary.csv --group chr --x pos --y p --label snp

# one panel per cytokine, with the test of each and BH across the cytokines
figkit plot --recipe feature_panel --data example/cytokine_panel.csv --x group \
  --y concentration_pg_ml --facet cytokine

# UMAP coloured by cluster, with silhouette widths
figkit plot --recipe embedding --data example/umap_clusters.csv --x UMAP_1 --y UMAP_2 --group cell_type
```

### The distribution recipe

The `distribution` recipe uses one set of bins for all the groups. The bin width
comes from the Freedman-Diaconis rule on the pooled values, 2 x IQR x n^(-1/3),
with type-7 quantiles. When the IQR is zero, the recipe uses the Sturges rule.
When the rule gives more than 200 bins, the recipe makes the bins wider so that
200 bins cover the range. The first bin starts at floor(min / width) x width, and
each bin includes its left edge.

The solid curve is a Gaussian kernel density with the bandwidth of
`stats::bw.nrd0`. Both engines sum the kernels directly on 512 points from
min - 3 bw to max + 3 bw. The dashed curve is a normal distribution with the
sample mean and SD. Both curves show counts, that is the density x n x bin width.
A panel with 500 values or fewer also shows a rug.

The horizontal axis shows the variable. `--xlab` names it, and `--ylab` names it
when `--xlab` is absent. The Shapiro-Wilk test runs on 3 to 5000 values. The
recipe does not compare the groups. Use `two_group_compare` or
`multi_group_compare` to test a difference.

### The composition recipe

The `composition` recipe sums the counts of each sample and category. Then it
divides each count by the sample total. Without `--y`, each row counts as one
observation. A sample must belong to one group. The figure sorts the samples by
group and then by name, and a trailing number sorts as a number.

With two groups, the recipe compares each category with the Mann-Whitney U test
on the per-sample proportions. With three or more groups, it uses the
Kruskal-Wallis test. The Benjamini-Hochberg method adjusts the P values across
the categories. Stars in the legend mark each category with an adjusted P below
0.05.

Proportions are compositional, because they sum to one in each sample. A rise in
one category therefore forces the other categories down. The tests are thus not
independent, and the methods paragraph calls them exploratory.

### The interaction geom

`--geom interaction` is a geom of `factorial_anova`. It draws the mean +/- SEM of
each cell, with one line for each `--fill` level across the `--x` levels. The
lines are dodged by 0.25, so the error bars do not overlap. A `--facet` factor
gives one panel for each level. The ANOVA, the subtitle and the stats table do not
change.

For every geom of `factorial_anova`, an `--x` level longer than 12 characters
wraps on the spaces. The figure width grows with the number of `--x` levels
multiplied by the number of `--fill` levels.

The four omics recipes follow these rules:

- `ma_plot` finds the columns of DESeq2, limma and edgeR tables. `--x` sets the
  mean column and `--y` sets the p column. A gene without a p counts as not
  significant. Above 20,000 grey genes, the recipe draws one grey point for
  each cell of a 600 by 600 grid, so the vector files stay small.
- `manhattan` puts the chromosomes in the order 1 to 22, X, Y, XY, MT and
  removes a "chr" prefix. A p outside 0 to 1 stops the recipe. A p of exactly 0
  becomes the smallest positive double for the statistics, and the figure
  draws it as a triangle at the top. Above 50,000 SNPs, the figure draws every
  SNP with p < 0.01 and one in every k of the others. Every statistic uses all
  SNPs.
- `feature_panel` uses one test family for all the features. With `--test
  auto`, the parametric test is used only when the Shapiro-Wilk test accepts
  every group of every feature. The recipe stops above 36 features and points
  to the `heatmap` recipe.
- `embedding` does not compute the embedding. A numeric `--group` with more
  than 12 values is a continuous score on the viridis scale, unless it holds 50
  or fewer whole numbers. Above 5,000 cells, the silhouette uses one in every k
  cells in input order. A silhouette on two embedding coordinates describes the
  separation in the embedding and not in the original feature space.

## The artifact bundle

Every `figkit plot` run writes one folder. That folder is the reproducible unit.

- `figure_*.pdf`, `.png`, `.svg` — the figure in each format.
- `qc_normality_*.png`, `.pdf` — a QQ panel per quantity the normality check ran on,
  for the recipes that choose a test from normality (absent for the rest).
- `script_*.R` or `script_*.py` — a standalone script that redraws the figure.
- `stats_*.csv` — the test result: statistic, p-value, adjusted p, effect size,
  group sizes. Both engines write a missing value as NA, a flag as TRUE or FALSE,
  an infinite value as Inf or -Inf, and a whole number with no decimal part.
- `input_*` — a copy of the exact data used, with the extension of the input
  (heatmaps also copy the annotation), so the bundle does not depend on a file that
  may change.
- `data_log.md` — the original filename, md5 and sha256 checksums, the read rule
  and each value that the read changed, the column types, and each cleaning step
  applied.
- `plot_config.yaml` — the editable spec that drives `figkit render`.
- `session_info.txt` — the full environment: R `sessionInfo()` or a Python
  package list.
- `packages_*.txt` — the package versions for this figure: the ones the script loads
  directly first, then the complete installed set, headed with the R or Python version
  and the container OS.
- `manifest_*.json` — structured provenance: the recipe, every plot option with
  `null` for an option that the run did not set, the engine that drew the figure,
  the container pinned by digest, the pubready git commit, and the versions of the
  main packages. `packages_*.txt` holds the full list. Numbers keep their full
  precision, and NaN and infinity are written as `null`, so the file is valid JSON.
- `methods_*.md` — a manuscript-ready methods paragraph, templated from the same
  manifest so the reported numbers always match the figure.
- `REPRODUCE.md` — the one Podman command to regenerate the figure, with the
  image digest.
- `report.qmd`, `report.html`, `report_provenance.json` — a Quarto report that
  embeds the figure, the stats table and the methods text, and its provenance. The
  report names the LLM and the agent harness that `--llm-provider`, `--llm-model`,
  `--agent-harness` and `--agent-harness-version` give.

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
| brackets | rstatix + ggpubr | statannotations (two groups), drawn by pubready (multi-group) |
| ANOVA / post hoc | rstatix | pingouin + scikit-posthocs + statsmodels |
| survival | survival / survminer; the at-risk table is drawn by pubready | lifelines; the at-risk table is drawn by pubready |
| heatmap | ComplexHeatmap | PyComplexHeatmap (seaborn clustermap fallback) |
| correlation | ggplot2 | seaborn `regplot` |
| correlation heatmap | ComplexHeatmap | seaborn `clustermap` |
| UpSet | ComplexHeatmap `UpSet` | matplotlib, in the ComplexHeatmap layout |
| Prism style (`--theme prism`) | ggprism `theme_prism` | matched matplotlib style |
| volcano / labels | ggplot2 + ggrepel | matplotlib + adjustText |
| raincloud geom | ggdist half-eye | matplotlib half-violin |
| PCA / PERMANOVA | prcomp + exact or seeded permutation | numpy SVD + exact or seeded permutation |
| dose-response (4PL) | `stats::nls` | `scipy.optimize.curve_fit` |
| ROC / DeLong | pROC | the rules of pROC, restated with numpy |
| meta-analysis | metafor `rma` | statsmodels `combine_effects` for the fixed effect; DerSimonian-Laird with numpy, tau2 truncated at zero |
| mixed ANOVA | split-plot sums of squares, the GG epsilon and Mauchly's test of car | the same formulas with numpy |
| histogram and density | Freedman-Diaconis bins, `stats::bw.nrd0`, a direct Gaussian kernel sum | the same bins, bandwidth and kernel sum with numpy |
| composition tests | `stats::wilcox.test`, `stats::kruskal.test`, `stats::p.adjust` | the `wilcox.test` p-value rule restated, `scipy.stats.kruskal`, the BH rule restated |
| estimation bootstrap | the Park-Miller generator and the BCa interval, written out | the same generator and formulas with numpy, so the resamples are identical |
| Michaelis-Menten | `stats::nls` | `scipy.optimize.curve_fit` |
| standard curve | the 4PL of `dose_response`, or `stats::lm` | the 4PL of `dose_response`, or least squares with numpy |
| competing risks | cmprsk `cuminc` | the `cinc` and `crstm` routines of cmprsk, restated with numpy |
| calibration | `stats::glm`, `stats::lowess` | statsmodels `GLM` and `lowess`, the groups of `hoslem.test` restated |
| Sankey diagram | ggplot2 polygons, `stats::mcnemar.test` | matplotlib polygons, statsmodels `SquareTable` |
| MA and Manhattan labels | ggplot2 + ggrepel, ggpubr `ggarrange` for the two panels | matplotlib + adjustText |
| QQ band, lambda_GC | `stats::qbeta`, `stats::qchisq` | `scipy.stats.beta`, `scipy.stats.chi2` |
| feature panel tests | `stats::t.test`, `wilcox.test`, `oneway.test`, `kruskal.test`, `p.adjust` | scipy with R's wilcox rule, BH restated with numpy |
| embedding silhouette | `cluster::silhouette` | the rules of `cluster::silhouette`, restated with numpy |

Where Python has no faithful equivalent, pubready routes you to R rather than
running something different. The one current case is the aligned rank transform
for non-parametric factorial ANOVA: `--test art --engine python` stops and points
to `--engine r`. Under `--test auto`, when the model residuals fail Shapiro-Wilk,
the R engine runs the aligned rank transform and the Python engine stops with the
same route. The `pca` PERMANOVA enumerates every labeling of the samples when there
are 5000 or fewer, so the P value is exact and the same in both engines. A larger
design uses 999 seeded permutations, and each engine draws its own, so the two P
values differ by the Monte Carlo error.

## Setup

pubready needs Podman and one built image.

```bash
figkit build          # builds localhost/pubready:3.0.0 from container/Containerfile
```

The image carries R (ggpubr, rstatix, ggprism, ggrepel, ggdist, survival,
survminer, ARTool, ComplexHeatmap, pROC, metafor) and Python (matplotlib, seaborn,
statannotations, pingouin, statsmodels, scikit-posthocs, lifelines,
PyComplexHeatmap, upsetplot, adjustText). The R packages come from the fixed
CRAN snapshot of the `rocker/r-ver:4.4.1` base image. `container/requirements.txt`
gives no Python versions, so the build writes the versions that it resolved to
`/opt/pubready/requirements.lock` in the image. A built image is therefore
pinned, and each manifest records its digest, but a new build can resolve newer
Python packages. The `core/` scripts are mounted at runtime, so editing a recipe
does not need a rebuild; only changing the dependencies does.

## Testing and versioning

Both engines carry a unit suite. `figkit test` runs them inside the image
against the mounted source, so the tested code is exactly what figkit runs: the
R suite (`core/r/tests/testthat`) loads the source package with `pkgload`; the
Python suite (`core/python/tests`) runs under `pytest`. The end-to-end
`tests/smoke_test.sh` drives every recipe on both engines.

The public benchmarks in `benchmarks/` check the numbers. Each recipe other than
`spider_response` and `waterfall` runs on a published dataset, for example the
`sleep` data of Student (1908), the NCCTG lung cancer trial, the BCG vaccine
trials or the PEFR data of Bland and Altman (1986). Both suites compare the result with a reference value that the
canonical R function computed, such as `stats::t.test`, `survival::coxph`,
`pROC::roc` or `metafor::rma`. `tests/benchmark_test.sh` runs the same cases
end to end through `figkit plot`. See `benchmarks/README.md`.

The version lives in the root `VERSION` file and moves in lockstep across
`core/r/DESCRIPTION`, `core/r/R/version.R`, `core/python/pubready/version.py`,
the default image tag in `core/r/R/run.R`, `core/r/R/bundle.R`,
`core/python/pubready/run.py` and `core/python/pubready/bundle.py`, and the
README install command.
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
adding a line with its options to `core/recipe_options.txt`, an example dataset, a smoke case in `tests/smoke_test.sh`, a reference
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
