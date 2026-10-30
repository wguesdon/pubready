# Changelog

All notable changes to pubready are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because reproducibility is the point of this tool, "breaking" here means either a
change to the `figkit` command surface or a change that alters figure output
(a default test, a post-hoc rule, or a plotting detail). Those bump the major
version; new recipes bump the minor; fixes that leave outputs unchanged bump the
patch. The repo version, both engine package versions, and the container image
tag move together.

## [Unreleased]

### Added

### Changed
- Before `--test auto` takes a rank test because the normality check failed,
  `two_group_compare`, `paired_compare`, `multi_group_compare` and
  `feature_panel` compute the smallest P that the rank test can reach at the
  group sizes. When that P is above 0.05, the recipe keeps the parametric test,
  and the methods text gives the attainable P. A forced rank test still runs,
  and the text gives the attainable P. `test_meta` records it as
  `rank_test_attainable_p` (F61).
- Under `--test auto`, `two_group_compare` runs Welch's t-test when both groups
  are normal. Levene's test still runs and goes to the stats table, but it no
  longer chooses Student's t-test. `--test t` also runs Welch's t-test, and
  `--test student_t` forces Student's t-test. The methods text says that the
  Welch t-test does not assume equal variances (F109).
- The comparison recipes order their groups with one rule in both engines:
  the `--levels` order, else numeric order when every label is a number, else
  code-point order. A numeric group reads "80", not "80.0", and a dose axis
  reads 0, 80, 160 in Python too. `--levels` now works in `two_group_compare`,
  `multi_group_compare`, `paired_compare`, `factorial_anova` (the `--x`
  levels), `estimation_plot` and `feature_panel`. The first level is group 1 of
  a contrast. `paired_compare` no longer takes the order of first appearance,
  so pass `--levels pre,post` to keep pre on the left. `estimation_plot` keeps
  `--ref` for the reference, and `proportions` keeps its own order rule with
  the new labels (F38, F64).
- Under `--test auto`, the Python engine of `factorial_anova` runs the
  Shapiro-Wilk test on the residuals of the full model, as R does. When
  P < 0.05, R takes the aligned rank transform, which Python does not have, so
  Python stops with a message that gives the P and routes to `--engine r` or
  to `--test anova`. Python ran the ANOVA before, so the engines reported
  different tests. Both engines write the residual P to the stats table as
  `residual_shapiro_p` and to `test_meta` (F37).

### Fixed
- `--levels`, `--markers` and `--covariates` accept "|" as the separator of
  the items, so an item can hold a comma, for example
  `--levels "none|edema, no diuretics"`. A value with no "|" splits on commas
  as before (F96).
- A sample with zero range is not normal in both engines, with a Shapiro P of
  NA, so `--test auto` chooses the same test. A t-test, a one-way ANOVA or
  Welch's ANOVA whose variance is zero stops with a message that names the
  groups. `correlation` stops on a constant column and names it.
  `feature_panel` reports a constant feature as "not tested: constant" and
  leaves it out of the normality check and of the BH count (F43).
- The paired t-test of `two_group_compare --paired` writes no sentence about
  equal variances and names its effect size "Cohen's dz", the mean of the
  within-pair differences over their SD, group1 minus group2 (F66).
- `two_group_compare --paired` chooses between the paired t-test and the
  Wilcoxon signed-rank test from the Shapiro-Wilk test of the within-pair
  differences, not of each group. The rows still pair by their order in each
  group, and groups of unequal size stop the run. The stats table gets
  `shapiro_p_diff`, the QC panel shows the differences, and the methods text
  names them (F36).
- The R engine of `two_group_compare` computes Cohen's d over the pooled SD, as
  the Python engine and `feature_panel` do, and not over the average variance
  of `rstatix::cohens_d`. On the cats data both engines now give -0.95. The
  methods text states the formula (F65).
- The R engine of `factorial_anova` makes its factors after it drops the rows
  with a missing value, so a level whose every outcome is missing leaves no
  empty slot, and the points sit on their boxes (F54).
- The R engine of `two_group_compare`, `multi_group_compare` and
  `factorial_anova`, also with `--test art`, accepts a column name with a space
  or a symbol. The tests and the emitted scripts run on a copy with syntactic
  names, and the tables keep the original names (F44).
- `two_group_compare`, `multi_group_compare`, `paired_compare`,
  `factorial_anova`, `estimation_plot`, `feature_panel`, `correlation` and
  `proportions` stop on a `--test` value that they do not know, with a message
  that lists the values that they know. `proportions` and `factorial_anova` ran
  their default test before (F63).
- `paired_compare` stops and names the subjects when a subject has two rows in
  one condition, as `estimation_plot --paired` does. It kept the first row
  before, so the result depended on the row order (F55).
- `two_group_compare`, `multi_group_compare`, `paired_compare`,
  `estimation_plot` and `factorial_anova` stop with "the outcome has zero
  variance" when every value of the outcome is the same. R stopped with an
  unclear error, and Python wrote "nan" or F values from rounding noise (F90).
- The comparison recipes and `correlation` stop with a clear message when the
  variance of the outcome or of a column is not a finite positive number, as
  for values near 1e300 or near 1e-300. They gave a silent NA, P = 1 or a
  crash before (F83).
- The data log of `factorial_anova` counts each reason for a dropped row on its
  own: the outcome values that became NA as numbers, and the rows with no value
  in each column. It named every column as the reason before (F124).
- `factorial_anova` counts each cell of the design before the fit and stops
  when a cell holds no value, and names the cell. The engines tested different
  hypotheses on such a design: the interaction gave P = 2.5e-7 in Python and
  P = 0.82 in R (F53).
- `factorial_anova` stops with a message that names the column and the count
  when a factor has only one level, or when every cell holds one value, so no
  residual degree of freedom is left (F92).
- The methods text of `factorial_anova --geom bar` says that the bars show the
  mean +/- SEM of each cell, and the text of `--geom interaction` says the same
  of the points (F122).
- On the chi-squared branch of `proportions`, a 2 x 2 table with a zero cell
  gives no odds ratio in both engines (NA), and the methods text says that a
  cell is zero. R wrote Inf and Python an empty cell before (F93).
- `multi_group_compare` writes a Tukey or Games-Howell P below 1e-12 as
  "< 1e-12" in the stats table and in a `p.format` bracket. The studentized
  range loses its precision there in both engines: R wrote 0 or about 2.4e-13,
  and Python wrote 0 (F118).
- The y axis of `two_group_compare`, `multi_group_compare` and `paired_compare`
  reads in scientific notation, for example "1.5e-09", when its largest break
  is below 1e-4, as for molar concentrations. It read "0.000000" before. The R
  bracket of `two_group_compare` now sits 6% of the data span above the drawn
  geom, so it no longer falls to 0 on such values (F84).
- A wide `feature_panel` keeps a column when 80% or more of its cells read as
  numbers. Each other cell, such as "<0.1", becomes NA and the data log
  records it, and the data log names each text column that is left out. The
  whole column was dropped with no record before (F42).
- `estimation_plot` wraps its title and the subtitle of the Cumming layout to
  the figure width, so a long reference name no longer runs past the edge of
  the R figure. `factorial_anova` wraps its legend labels and title at 20
  characters in both engines, and the figure width follows the wrapped labels
  (F81).

## [3.0.0] - 2026-10-29

The fixes of the public data audit change figure output and stats tables in
both engines, and the rule at the top of this file calls for a major version for
such changes. The release also adds six benchmark cases, the
`proportions --levels` option and the `pubready.fisher` module of the Python
engine.

- The default image tag is `localhost/pubready:3.0.0`. Run `figkit build` once
  after the upgrade. The README installs the `v3.0.0` tag with uv.

### Added
- Six benchmark cases cover test paths that no case covered before: Fisher's
  exact test of a 2 x 2 table (`mtcars_fisher`) and of a 3 x 2 table
  (`infert_fisher`), Kendall's tau with ties (`swiss_kendall`), Welch's one-way
  ANOVA (`chickwts_welch`), the paired Wilcoxon signed-rank test
  (`sleep_wilcoxon`) and pooling on standard errors (`normand_meta`). The data
  come from the datasets and metadat packages, and the exporter writes 42
  datasets and 257 reference values.
- `proportions` accepts `--levels "a,b,c"` to set the order of the groups, as
  `sankey` does for its states. The list must name every group of the data.
  Without it, the groups keep the order of their first appearance, which drew
  the `titanic` classes as 3rd, 1st, 2nd, Crew.

### Changed
- The stats tables of `two_group_compare`, `paired_compare` and
  `multi_group_compare` have a new `contrast` column, for example "1 - 2", and
  the methods paragraph names the direction of the effect size, for example
  "Cohen's d = -0.83 (1 minus 2)". The recipes keep their sign conventions:
  `two_group_compare` reports group 1 minus group 2, `paired_compare` reports
  condition 2 minus condition 1, and the Tukey estimate of
  `multi_group_compare` is group 2 minus group 1. On the same `sleep` data,
  Cohen's d is -0.83 and Cohen's dz is +1.28, and before this change no text
  said why the signs differ.
- `feature_panel` reports eta-squared (H), (H - k + 1) / (n - k), for the
  Kruskal-Wallis test in both engines, where it reported epsilon squared,
  H / (n - 1). `multi_group_compare` already reports eta-squared (H), so the two
  recipes now give one value for the same test. The values change, for example
  from 0.651 to 0.646 for the sepal length of `iris`.
- A Spearman or a Kendall `correlation` scatter has no least-squares line and
  no band now. A rank correlation measures a monotonic association, and on the
  `airquality` benchmark the line predicted a negative ozone level below 60 F.
  The methods paragraph says that no line is drawn. A Pearson scatter keeps the
  line and its 95% band.
- The stats tables of `multi_group_compare`, `pca`, `proportions` and `sankey`
  now hold the omnibus test on each row: `omnibus_test`, `omnibus_statistic`,
  `omnibus_df1`, `omnibus_df2` and `omnibus_p`. `multi_group_compare` adds the
  omnibus effect size, `pca` adds the PERMANOVA method and the count of
  permutations, and `proportions` adds the odds ratio of a 2 x 2 table. Before,
  only the manifest and the figure title held these numbers.
- The `proportions` methods paragraph now gives the chi-square statistic and
  its df, and it names the direction of the odds ratio, for example "The odds
  ratio of Admitted for Male versus Female was 1.84".
- The Python stats tables wrote a count or a df as a float, for example 16470.0
  and 2.0, when the column also held a missing value. They now write 16470 and
  2, as the R tables do.
- The Python `timecourse` spread the group means at each time point k times
  wider than R, because each group took a full dodge step. The groups now
  share one dodge width, as `position_dodge` does in R.
- The `composition` facet strips now name the group column, for example "Cont:
  High". On the `housing` benchmark the strips read "High" and "Low", the same
  words as two levels of the stacked categories.
- The `cox_forest` subtitle now warns when the proportional-hazards check fails,
  for example "137 subjects, 128 events; PH assumption violated (P < 0.001)" on
  the `veteran` benchmark. Before, only the methods paragraph said so.
- The Python `heatmap` cut a legend at its right edge when a level name was
  longer than the width estimate of PyComplexHeatmap, for example the `cell`
  legend of the `airway` benchmark. The legends now draw in full, and the saved
  figure includes them.
- `correlation_heatmap` ordered the variables differently in the two engines,
  for example carb, hp, wt in R and wt, cyl, disp in Python on `mtcars`. It now
  uses the leaf rule of `heatmap` in both engines, with the mean correlation of
  each variable as the weight. `test_meta` records the order, the emitted
  scripts carry the rule, and the smoke test checks that the engines agree.
- The Python `cox_forest` checked the proportional-hazards assumption with the
  smallest per-covariate P of lifelines on the rank transform, where R uses the
  global `cox.zph` test. Python now computes the global `cox.zph` test with the
  "km" transform, as the C routine of the survival package does. On the `lung`
  benchmark the Python P was 0.13 and is now 0.2156, the R value. The benchmarks
  check the new `ph_global_p` reference of `lung_cox` and `veteran_cox`.

### Fixed
- The `multi_group_compare` methods paragraph wrote the Welch df unformatted in R,
  "F(5, 29.9520363861042)", and rounded to a whole number in Python, "F(5, 30)".
  Both engines now write a whole df as an integer and a fractional df with two
  decimals, "F(5, 29.95)". The new `chickwts_welch` case showed it.
- The Python `proportions` ran Fisher's exact test only on a 2 x 2 table and
  fell back to the chi-square test for a larger table, where R runs Fisher's
  test on any table. For a 2 x 2 table it reported the sample odds ratio with a
  Wald CI, where R reports the conditional maximum likelihood estimate with its
  exact CI. A new `pubready.fisher` module follows `stats::fisher.test`: it
  carries a port of the `zeroin` root finder of R, so the estimate and the CI
  match R, and it enumerates the tables of an r x c table. Above one million
  tables it stops, and the "Engine limits" section of the adapters records the
  gap.
- `paired_compare` wrote "The test was chosen from the Shapiro-Wilk normality
  of the within-subject differences" also when `--test` set the test. The
  methods paragraph now says "The test was set by the user" for a set test. On
  the `sleep` benchmark, the old text contradicted the paired t-test, because
  the Shapiro-Wilk P of the differences is 0.033.
- `two_group_compare` wrote "A non-parametric test was used because the normality
  assumption was not met" also when `--test wilcoxon` set the test. For a set
  test, the methods paragraph now says "The test was set by the user".
- `forest_meta` wrote "the DerSimonian-Laird tau2 was truncated at zero" in every
  methods paragraph. It now writes the clause only when Q is below its degrees
  of freedom, which is the only condition that truncates tau2. On the `bcg`
  benchmark, tau2 is 0.309, and the old text read as a truncation.
- The `volcano` methods paragraph wrote "p <= 0.05" also when the cutoff read
  an adjusted p column such as `padj`. It now names the column, as the axis and
  the legend do. It also states the label rule: the significant features with
  the largest product of the absolute log2 fold change and -log10 p are
  labelled. The R and the Python text now print the cutoffs in the same form.
- The `dose_response` and `enzyme_kinetics` methods paragraphs did not say that
  the points and the error bars show the mean +/- SEM of the replicates. They
  now say what each point shows, also for data without replicates.
- The `estimation_plot` methods paragraph now says that the P values against the
  reference carry no adjustment for multiple comparisons when two or more groups
  are compared.
- The `heatmap` methods sentence on the per-feature test was broken ("A
  per-feature Welch t-test (trt vs untrt), BH-adjusted was applied per feature;
  38 feature(s)"). It is rewritten, and it now says that the Welch test treats
  the samples as independent, so a reader sees the assumption for a paired
  design.
- The `cox_forest` and `survival_km` methods paragraphs wrote "P = < 0.001".
  They now write "P < 0.001". A new `p_phrase()` helper in each engine writes
  the P value of a methods sentence.
- The R `survival_km` methods paragraph padded each median with spaces, for
  example "adeno (median   51)". The padding is removed, so the text matches the
  Python engine.
- `dose_response`, `enzyme_kinetics`, `bland_altman`, `timecourse` and
  `estimation_plot` dropped trailing zeros in their labels and methods text, for
  example "EC50 = 3.7 (95% CI 3 to 4.58)" and "P = 0.01". A new `sig_text()`
  helper in each engine keeps a fixed count of significant digits, so the text
  now reads "EC50 = 3.70 (95% CI 3.00 to 4.58)" and "P = 0.010".
- Three methods templates had wording errors: `correlation` wrote "The
  association between Wt and mpg" with a capital in mid-sentence,
  `correlation_heatmap` wrote "pearson" in lower case, and `enrichment_dot` wrote
  "a over-representation".
- The `standard_curve` methods paragraph wrote "6 of 8 levels fell within 80% to
  120%" when 2 levels were not read back at all, because their mean signals lie
  outside the signal range of the curve. The text read as two failed standards.
  It now counts only the levels that were read back and names the others.
- The Python `correlation` stats table now has the statistic and the interval of
  `stats::cor.test`, as the R table has. For Spearman it writes S and for
  Kendall T or z, where it wrote nothing before. For Pearson it writes the
  Fisher z CI at full precision, where it wrote the 2-digit CI of pingouin. For
  Spearman and Kendall it writes no CI, as R does. On the `airquality` and
  `mtcars` benchmarks, the two tables now agree to 1e-15.
- The Python `multi_group_compare` gave no effect size for the Kruskal-Wallis
  test, where R gives eta-squared (H). It now gives the same value.

## [2.0.0] - 2026-10-29

The fixes below change figure output in both engines, and the rule at the top of
this file calls for a major version for such changes. The release also adds 12
recipes, the interaction geom, the uv install, `figkit doctor` and the
`sankey --levels` option.

### Added
- pubready installs with uv from its Git URL, as a uv tool or as a dependency
  of a uv project. The root `pyproject.toml` builds a package whose `figkit`
  launcher runs the same `cli/figkit` and the same pinned image. The launcher
  reads the commit that uv recorded at install time, so each bundle still names
  the exact code. `figkit where` prints the folder that holds the adapters,
  `reference/` and `docs/` of the install. `tests/install_test.sh` installs HEAD
  with uv and checks a figure from the installed command.
- `sankey --levels "a,b,c"` sets the state order, for example a clinical order
  of the response categories, in both engines. The list must name every state
  of the data, else the recipe stops and names the missing states. A listed
  state that no row holds is left out. `plot_config.yaml` keeps the list, so
  `figkit render` draws the same order, and the methods paragraph states it.
- `figkit doctor` checks the host setup. It checks the platform, bash, `podman`
  on the PATH, the macOS Podman machine, the Podman service, rootless mode, the
  subordinate UID ranges on Linux, the mounted folders, the image and a test
  container with the mounts of `figkit plot`. Each failed check prints its fix.
  The command never installs anything and never runs `sudo`, and it exits 1 when
  a check fails. `tests/doctor_test.sh` tests each result with stub `podman`
  and `uname` scripts, so it needs no image.
- `tools/render_audit/` writes stress datasets and renders every recipe with
  both engines, and it puts the two figures of each case side by side.
- The `distribution` recipe, in both engines, shows the shape of one numeric
  variable, overall or for each group. The figure is a histogram with a kernel
  density and a normal fit, with one panel for each group on a shared x axis.
  The bins come from the Freedman-Diaconis rule on the pooled values, and the
  bandwidth is that of `stats::bw.nrd0`. The stats table gives the moments, the
  quantiles and the Shapiro-Wilk test of each group. The recipe does not compare
  the groups.
- The `composition` recipe, in both engines, shows one 100% stacked bar for each
  sample, for example the cell-type proportions of a single-cell experiment. With
  `--group`, each category is compared with the Mann-Whitney U test for two
  groups or the Kruskal-Wallis test for three or more. The P values are adjusted
  with the Benjamini-Hochberg method. The methods paragraph states that
  proportions are compositional, so the tests are exploratory.
- `factorial_anova --geom interaction`, in both engines, draws the interaction
  plot of the same ANOVA. It shows the mean +/- SEM of each cell, with one line
  for each `--fill` level. The statistics do not change.
- Three public benchmark cases: `faithful_distribution`,
  `tooth_growth_distribution` and `housing_composition`.
- Three bench assay recipes, in both engines:
  - `estimation_plot` shows the mean difference of each group against a
    reference with a 95% BCa bootstrap CI. Two groups give the Gardner-Altman
    layout, and three or more give the Cumming layout. The two engines draw the
    same resamples from the Park-Miller generator, so the intervals agree. The
    `--ref` option names the reference. `--paired --id` resamples pairs and
    draws a slopegraph.
  - `enzyme_kinetics` fits the Michaelis-Menten equation to each curve and
    reports Vmax and Km with 95% CIs. With two or more curves, extra-sum-of-squares
    F tests compare one shared Km and one shared Vmax with the separate fits.
  - `standard_curve` fits a 4PL or a straight line to the standards and reads
    each unknown off the curve. An unknown outside the range of the standards
    gets a flag and no concentration. The recovery of each standard level is
    the QC of the fit.
- Three public benchmark cases for these recipes: `chickwts_estimation`,
  `puromycin_kinetics` and `dnase_standard_curve`.
- Three clinical recipes, in both engines:
  - `cumulative_incidence` estimates the cumulative incidence of each cause of a
    competing-risks outcome with the Aalen-Johansen estimator. Gray's test
    compares the groups for each cause. The R engine calls `cmprsk::cuminc`, and
    the Python engine restates the `cinc` and `crstm` routines of cmprsk, so tied
    times give the same numbers. The `--cause` option names the cause of
    interest.
  - `calibration` gives the Brier score, the calibration intercept and slope
    with Wald CIs, the Hosmer-Lemeshow test and the observed/expected ratio of
    one or more risk models. The groups of the test follow
    `ResourceSelection::hoslem.test`. The figure shows the decile points, a
    lowess smooth and a spike histogram of the predictions.
  - `sankey` draws the state of each patient at each visit and the transitions
    between visits. With two visits, Bowker's test of symmetry, or McNemar's
    test for two states, asks whether the patients move more in one direction.
- Example datasets for the new recipes: `cumulative_incidence.csv`,
  `calibration.csv` and `sankey.csv`, each from a seeded generator in
  `tools/make_examples/`.
- Three public benchmark cases: `mgus2_cif` from `survival::mgus2`,
  `pima_calibration` from `MASS::Pima.te`, and `pbcseq_sankey` from
  `survival::pbcseq`. They add 3 datasets and 24 reference values.
- Four omics recipes, in both engines:
  - `ma_plot` draws the log2 fold change against the mean expression of a
    DESeq2, limma or edgeR table. Up genes are red, down genes are blue and the
    top genes by p carry labels. The stats table lists the significant genes
    sorted by p, and the subtitle gives the counts.
  - `manhattan` draws a Manhattan plot and a QQ plot of GWAS summary
    statistics side by side. The recipe reports lambda_GC, the counts below
    5e-8 and 1e-5 and the lead SNP of each chromosome. It reads a p of exactly
    0 as the smallest positive double and thins the SNPs with p of 0.01 or more
    for the vector outputs.
  - `feature_panel` draws one panel per feature with one test family for the
    whole panel: the Welch t-test or the Mann-Whitney U test for two groups,
    and the one-way ANOVA or the Kruskal-Wallis test for more. The p values are
    adjusted across the features with the Benjamini-Hochberg method.
  - `embedding` draws UMAP or t-SNE coordinates from upstream, coloured by
    cluster, with a label at each cluster median. It reports the mean
    silhouette width of each cluster on the two coordinates.
- Example datasets for the new recipes: `de_basemean.csv`, `gwas_summary.csv`,
  `cytokine_panel.csv` and `umap_clusters.csv`, each from a seeded generator in
  `tools/make_examples/`.
- Benchmark cases `airway_ma`, `gwas_manhattan` (the `gwasResults` data of
  qqman), `iris_feature_panel_kruskal`, `iris_feature_panel_anova` and
  `iris_embedding`. The `airway_de.csv` file now holds the `baseMean` column.

### Changed
- The default image tag is `localhost/pubready:2.0.0`. Run `figkit build` once
  after the upgrade. The README installs the `v2.0.0` tag with uv.
- `factorial_anova` wraps an `--x` level longer than 12 characters, and the
  figure width grows with the number of `--x` levels multiplied by the number of
  `--fill` levels. A design that fitted the old size keeps it.
- `upset` drew every intersection. Seven sets gave 89 bars in a figure too wide
  to read. Both engines now draw the 15 largest intersections, and `--top_n`
  sets the number. The title states how many intersections are not shown. The
  stats table still lists every intersection, and `n_intersections_shown` in the
  manifest gives the number drawn.
- A `dose_response` curve is "not determined" in six cases. Two checks before
  the fit find a flat curve and give it no fit. The first finds a curve whose
  mean responses vary by less than 10% of the range of the mean responses of
  all curves. The second finds a curve whose responses do not differ between
  the doses in a one-way ANOVA (P >= 0.05), and it also covers a single curve.
  It needs a replicate at one dose or more. The checks use no fitted value, so
  both engines give the same reason for a flat curve. Before
  the check, R `nls` failed on a flat curve and Python converged, so the two
  reasons differed. A curve is also not determined when its fit fails, when its
  fitted span is below 10% of the range of the mean responses, when its EC50 is
  outside the tested doses, or when its 95% CI is not finite. The stats table
  gives EC50 NA, a `fit_status` and a `status_reason` for that curve. The
  figure draws its mean and SEM without a curve, and the legend labels it
  "(not determined)". The shared-midpoint F test leaves it out, and the
  subtitle and the methods paragraph say so.
- `multi_group_compare` shows a compact letter display when more than six
  pairs differ. Groups that share a letter do not differ. Both engines compute
  the letters with the same insert-and-absorb algorithm of Piepho (2004).
  `test_meta` records the display and the letters.
- `correlation` shows one subtitle in both engines, with the method, the
  coefficient symbol and the P value, for example "Pearson correlation,
  r = 0.71, P = 2.1e-10". R no longer uses the ggpubr `stat_cor` label.
- `survival_km` puts the legend to the right of the panel in both engines, and
  wraps the group names at 24 characters. The risk table has a 10 pt title over
  9 pt numbers, and a short line in the colour of each group. `test_meta`
  records the at-risk counts.
- The comparison recipes print the y axis as plain numbers with a comma between
  the thousands, for example 300,000. The x labels wrap at 12 characters. The R
  figure of `multi_group_compare` grows in width with the number of groups, as
  the Python figure does.

### Fixed
- `figkit` told the user to build the image when Podman was not installed, and
  also when Podman could not reach its service, such as a stopped macOS Podman
  machine. A missing `podman` now gives the install command of the platform. A
  Podman that does not answer now points to `figkit doctor`. Only a missing image
  gives the hint `figkit build`. `figkit version` prints "podman: not found" in
  place of an empty version.
- The points of `two_group_compare`, `multi_group_compare` and
  `factorial_anova` moved on each run, because R `geom_jitter`,
  `position_jitterdodge` and seaborn `stripplot` drew random offsets with no
  seed. Each point now sits at a fixed golden-ratio offset by its input row in
  its group or cell, the rule of the raincloud points. Both engines use the same
  offsets, so the two-group and multi-group points sit at the same places in R
  and Python. The emitted scripts use the same rule.
- On a host with SELinux in enforcing mode, such as Fedora, the containers of
  `figkit` could not read the mounted folders. Each `podman run` of `figkit`,
  the `REPRODUCE.md` command of each bundle and the maintainer commands now
  pass `--security-opt label=disable`. Podman does not relabel the user's
  files, and the option does nothing on a host without SELinux.
- The Python `correlation` band came from the bootstrap of seaborn `regplot`,
  so it moved on each run and differed from the R band. The Python engine now
  draws the t-based 95% CI of the fitted mean on 80 points, the band of R
  `geom_smooth(method = "lm")`, and the emitted script draws the same band.
- The `heatmap` row order differed between the engines. ComplexHeatmap
  reordered each branch by the row means of the drawn values, which are
  rounding noise after a row z-score, and PyComplexHeatmap kept the scipy
  order. Both engines now build the same complete-linkage tree and put the
  branch with the higher mean of the unscaled values first, with the earlier
  input row or column first on a tie. `test_meta` records `row_order` and
  `column_order`, and the smoke test checks that they agree. The Python emitted
  script used the average linkage of seaborn, and it now uses the same tree.
- The emitted scripts of `two_group_compare` and `multi_group_compare` drew a
  box for the bar and the raincloud geoms, in both engines. Each script now
  draws the geom of the figure, with the same layer settings, and it leaves out
  the points when the figure has none. A new test in each engine runs the
  emitted script of each geom and checks that it writes the figure.
- The Python `dose_response` legend drew a fitted curve as a line only. The key
  now shows the point on the line, as the R key does, and the mean points have
  the size of the R points.
- The Python `embedding` points were 27% wider than the R points, because the
  size conversion used 2.7 points for each R size unit. R draws a filled point
  of size 1 with a diameter of 0.75 x 72.27 / 25.4 = 2.13 points, and the
  Python engine now uses the same diameter. The legend keys and the emitted
  script use it too.
- The Python `embedding` labels moved between runs. adjustText stops after one
  second when it has no iteration limit, so the layout depended on the speed of
  the machine. The recipe now sets 2000 iterations, as `ma_plot` and
  `manhattan` do.
- `figkit inspect` and the other data commands stopped with "unbound variable"
  under the bash 3.2 of macOS when they had no arguments. `cli/figkit` now runs
  on bash 3.2.
- The R `calibration` figure dropped the 95% CI bar of a decile point whose
  observed proportion was 1, because the Wilson bound passed 1 by a rounding
  error. The bounds now stay inside 0 and 1, and the y scale keeps them.
- `estimation_plot` with two groups had a fixed width, so long group labels ran
  into each other in the Python figure. The Gardner-Altman figure now grows
  with the longest label line, by the same rule in both engines.
- The Python `factorial_anova` ignored `--geom bar` and drew a box plot. It now
  draws the mean of each cell as a bar with its SEM and the points over it, as
  the R engine does.
- The R `cox_forest` drew its figure with `survminer::ggforest`. The text was
  about 6 pt, grey bands crossed the rows, and a long level name ran into its
  hazard ratio. Both engines now draw the same house-style forest plot. Each
  model term has one row, and the first level of each factor has a "reference"
  row. The "HR (95% CI)" and "P" columns are on the right, and the subtitle
  gives the number of subjects and events. The numbers do not change. The
  generated R script no longer loads survminer.
- The hazard-ratio axis of `cox_forest` showed only 0.5 and 1 when a CI reached
  1.94. The axis now has 1-2-5 ticks that cover the CIs on both sides of 1, in
  both engines.
- The Python `heatmap` put the group legend and the z-score colour bar on top of
  the row names. PyComplexHeatmap measured the names before the first draw, when
  they were still empty. The legends now go to the right of the row names.
- The Python `heatmap` did not draw the significance stars that R draws when the
  annotation has two groups. It now writes them between the heatmap and the row
  names, with the same BH-adjusted Welch t thresholds as R.
- The two `upset` engines did not agree on the layout. R rotated the counts over
  the bars by 45 degrees and the set-size axis labels by 90 degrees, and it gave
  no set-size numbers. Python put the set sizes on the left, where the numbers
  overlapped the bars, and it sorted tied intersections in a different order.
  Python now draws the plot with matplotlib in the layout of ComplexHeatmap,
  with the set sizes on the right. Both engines write horizontal counts and
  set-size numbers. Both sort the intersections by size, then by degree, then by
  set name, and the stats table uses the same order.
- The R `enrichment_dot` gave the whole width to long GO term names, and the dot
  panel shrank to nothing with 25 terms of up to 90 characters. Both engines now
  wrap each label into balanced lines of at most 45 characters. The figure
  height grows with the number of label lines, and the width grows with the
  longest line.
- The Python `enrichment_dot` put the count legend inside the plot area, where
  it could cover dots. The count legend and the colour bar now sit right of the
  panel, as in R.
- The R `volcano` legend showed an "FC" key with no point when no feature passed
  the fold change alone. Both engines now give a key only to a class that holds
  at least one feature, and the key glyphs are opaque.
- The Python `volcano` labels collided, for example Gene15014 with Gene12099.
  adjustText ran before the final layout, so the later resize undid its
  spacing, and it stopped after a time limit. The labels are now placed after
  the layout is final, they repel from every point as ggrepel does, and a fixed
  iteration limit gives the same result on every run.
- The `bland_altman` labels "Bias" and "+1.96 SD" sat inside the data and
  covered points in both engines. They now sit in the right margin, outside
  the panel and level with their lines. The figure is 1 inch wider to give
  them room. The R y axis showed only two labelled ticks on the example data,
  and it now has about six.
- The Python `roc_curve` legend showed only the line of each marker, and the
  tick label "1.0" ended at the right edge of the image. The legend keys now
  show the hollow Youden point on the line, as in R, and the figure keeps a
  margin right of the axis.
- `swimmer` listed the arms in the order of the input rows, for example Arm C,
  Arm B, Arm A. Both engines now sort the arms as the other recipes sort their
  groups, so the colours and the legend follow that order. The R legend text
  was 10.4 pt and the title 13 pt. They are now 9 pt and 10 pt, as in Python.
- `forest_meta` in Python now truncates the DerSimonian-Laird tau2 at zero, as
  metafor does. Before, a homogeneous set of studies gave a negative tau2, a
  negative I2, negative study weights and a wrong random-effects estimate. On
  the 25-study stress case, Python now gives tau2 = 0, I2 = 0% and
  0.79 [0.71, 0.87], the same as R.
- `survival_km` in Python counts the number at risk as the subjects whose time
  is at or after each time point, as R and the definition do. Before, the
  lifelines helper counted the subjects after the time point. On
  `example/survival_trial.csv`, the treated arm now has 21 at risk at month 10
  in both engines.
- `survival_km` in R now shows the x axis up to the largest time, so the last
  event is on the figure.
- `paired_compare` in R now labels P < 1e-4 with four stars, as the other
  recipes and the Python engine do. The new `p_stars()` helper of `util.R`
  gives the same scale as `util.p_stars` of Python.
- `dose_response` in R no longer stops when the fit of one curve fails. Python
  no longer reports a meaningless EC50 for a flat curve.
- `two_group_compare --geom raincloud` in R puts the points to the left of the
  box, as Python does. The points no longer hide the box.
- `two_group_compare --geom violin` in R puts the bracket above the drawn
  violin. The bracket no longer cuts through the tail of the violin.
- `multi_group_compare` spaces the brackets by the span of the drawn axis. A
  star label no longer touches the next bracket line, also with `--geom bar`.
- `proportions` in Python stacks the bars in the order of the legend, as R
  does. The subtitle no longer touches the 100% tick label.
- Long category labels no longer overlap on the x axis of `two_group_compare`,
  `multi_group_compare`, `paired_compare` and `proportions`.
- `forest_meta` gives both engines the same 1-2-5 ticks on the log axis, with
  at least one tick on each side of the null line.

## [1.0.0] - 2026-10-07

This is the first major release. The rename from pubplot to pubready changes
the `figkit` surface, and the fixes below change figure output. The rule at the
top of this file calls for a major version for each of these changes.

### Added
- Seven recipes, in both engines, which take the catalog from 15 to 22:
  - `dose_response` fits a four-parameter logistic curve to each curve and
    reports the EC50, or the IC50 of a falling curve, with a 95% CI. With two or
    more curves, an extra-sum-of-squares F test asks whether one shared midpoint
    fits as well as one per curve. The shared fit is a profile over the midpoint
    with the plateaus solved exactly, so both engines reach the same optimum when
    the curves sit far apart.
  - `roc_curve` draws the ROC curves of one or more markers, with each AUC, its
    DeLong 95% CI, the Youden threshold and DeLong's test of each marker against
    the first. R uses pROC, and Python restates the rules of pROC. The new
    `--positive` option names the case level.
  - `bland_altman` gives the bias and the 95% limits of agreement of two methods,
    each with a CI (Bland and Altman 1999), and a regression test of proportional
    bias. The differences go to the QC panel.
  - `forest_meta` pools study estimates by inverse variance, as a fixed-effect
    and a DerSimonian-Laird random-effects model, with Cochran's Q, I2 and tau2.
    R uses metafor and Python statsmodels. The new `--ci_low`, `--ci_high`,
    `--se` and `--measure` options describe the input.
  - `timecourse` runs a mixed ANOVA (group between, time within) or a
    repeated-measures ANOVA, with the Greenhouse-Geisser correction when
    Mauchly's test rejects sphericity, and draws the mean +/- SEM per group over
    time. `--test gg | uncorrected` forces the choice.
  - `waterfall` draws the best change from baseline of each patient as sorted
    bars with the RECIST lines, from one row per visit or one row per patient.
  - `swimmer` draws the time on treatment of each patient as one sorted bar,
    with an arrow for an ongoing course. The new `--markers` option names the
    columns of event times, such as a response or a progression, that the bars
    carry as symbols.
- Example datasets for the new recipes: `dose_response.csv`, `biomarker_roc.csv`,
  `method_agreement.csv`, `meta_trials.csv`, `tumor_growth.csv` and
  `swimmer_patients.csv`, each from a seeded generator in `tools/make_examples/`.
- Public benchmarks. `benchmarks/` holds 29 dataset files and 122 reference
  values from the canonical R functions. `figkit test` now mounts the folder, and
  both unit suites compare every recipe with these values. `tests/benchmark_test.sh`
  runs the 31 cases end to end through `figkit plot` with both engines.
- Nine of these cases use real research data. The airway RNA-seq study of Himes
  et al. (2014) gives a DESeq2 table for `volcano`, a gene heatmap, a gene PCA, a
  GO over-representation table for `enrichment_dot` and a GO gene-set matrix for
  `upset`. Before these cases, `volcano`, `enrichment_dot` and `upset` had no
  benchmark. `MASS::anorexia` gives a paired design, `survival::veteran` gives a
  Kaplan-Meier and a Cox case with a four-level factor, and `state.x77` gives a
  correlation matrix with a space in two column names.
- The image adds pROC and metafor in a late layer, and GO.db in a later layer.
  Only the benchmark exporter uses GO.db.
- `spider_response`, the fifteenth recipe. It expresses each patient's
  measurement as a percentage change from that patient's own baseline, which
  is the first measurement in time, and draws one line per patient over the
  visits. Two dashed reference lines carry the RECIST 1.1 thresholds, +20%
  and -30%, and `--pd_threshold` and `--pr_threshold` move them. The stats
  table holds one row per patient with the baseline, the visit count, the
  best, worst and last change, and the category those put the patient in.
  The recipe runs no test, in the way `upset` runs none.
- The category reads the change from baseline alone, and the column is named
  `category_from_change` for that reason. A full RECIST 1.1 assessment also
  reads the nadir, the non-target lesions and any new lesion, and the methods
  paragraph says so rather than letting a reader assume otherwise.
- `example/spider_response.csv`, 20 patients in two arms measured at weeks 0,
  8, 16, 24 and 32, with three patients who stop after week 16. Generated by
  `tools/make_examples/make_spider_data.R` (seeded).
- Every bundle writes `packages_<ts>.txt` beside the figure and script: the packages the
  generated script loads directly first, then the complete installed set, each with its
  version and headed with the R or Python version and the container OS. R's
  `session_info.txt` listed only loaded packages, so this is the full environment record too.
- `tools/make_examples/check_determinism.sh`: regenerates the reference bundles in
  a temp tree and confirms every stats table matches the committed copy, so "same
  input, same numbers" is checkable rather than assumed.
- `reference/decision_tree.md` and `reference/assumptions.md`: the test-selection
  tree and assumption-check notes that the adapters and PRD point to. Written to
  match the recipe code, closing the gap where `reference/` held only a stub.
- Normality QC. `figkit diagnose` recommends a test and writes a QQ plot (per group,
  per variable, paired differences, or ANOVA residuals) labelled with the Shapiro-Wilk
  p, skewness, and kurtosis, without drawing the figure; every `figkit plot` bundle now
  carries the same `qc_normality_*` panel. The adapters walk `reference/decision_tree.md`
  and run `diagnose` when no test is named, so an unsure scientist is guided to the
  choice and shown the evidence. Both engines. Additive: no existing output changes.

### Changed
- The project is renamed from pubplot to pubready, because the GitHub repository
  is now `wguesdon/pubready`. The change is breaking. The default output
  directory is `pubready_output/`, the house theme value is `pubready_house`,
  the default image tag is `localhost/pubready:<version>`, and each `PUBPLOT_*`
  environment variable is now `PUBREADY_*`. The R package, the Python package
  and the container paths under `/opt/pubready` carry the same name. The command
  surface keeps the name `figkit`. A caller that sets `--theme pubplot_house`,
  that reads `pubplot_output/`, or that exports a `PUBPLOT_*` variable must
  change to the new name.
- The Claude Code skill is `pubready` and the opencode command is `/pubready`.
  The two earlier names do not resolve.
- Moved the example data generators and `generate_expected.sh` out of `example/`
  into `tools/make_examples/`. `example/` now holds only the input datasets, the
  `expected/` reference bundles, and its README, so the folder a reader opens
  first shows examples, not maintainer tooling.
- The committed volcano example uses a small synthetic DE table (1200 genes) in
  place of the 928 KB airway DESeq2 result. The airway benchmark stays available
  as the optional `tools/make_examples/make_volcano_airway.R`. Committed footprint
  drops by roughly 5 MB and the volcano figure is unchanged in kind. Output is not
  affected for any other recipe.
- The two engines report the same numbers on the two-group and paired recipes.
  Three differences are removed. `two_group_compare` in R checks equal variance
  with Levene's test rather than the F test, which is what the multi-group recipe,
  both Python recipes and the documentation already used. Both engines apply the
  p-value rule of `stats::wilcox.test` to the rank tests, so an untied sample
  below n = 50 takes the exact route in Python as well. The R stats table keeps
  the full-precision p value rather than the value that `rstatix` rounds to three
  significant digits. A caller who compares a new number with an older bundle
  will see the difference, which is largest for the Levene column.
- `paired_compare` in R no longer forces the normal approximation, so it agrees
  with the Python engine and with `two_group_compare --paired`.
- The methods paragraph of `two_group_compare` gives the effect size value and
  not only its name, as `paired_compare` already did.
- `tests/smoke_test.sh` asserts the parity. It runs the Mann-Whitney and the
  signed-rank cases on both engines and compares the p value and the effect size
  to six significant digits.
- The default image tag is `localhost/pubready:1.0.0`.

### Fixed
- The public benchmarks found five defects in the existing recipes. Each one
  changes a number in a stats table, so a caller who compares with an older
  bundle will see the difference.
- `factorial_anova` in R wrote the F statistic and the P value that
  `rstatix::anova_test` rounds to three decimals and to three significant digits.
  On ToothGrowth the table held F = 92 where `car::Anova` gives 91.99996. The
  parametric path now reads the Type II table from `car::Anova` at full precision.
- The two engines reported different effect sizes in `factorial_anova`. R reported
  the generalized eta squared, which equals the partial eta squared in a
  between-subjects design, and Python divided by the total sum of squares. Both
  engines now report the partial eta squared, name it in a new `effect_size_name`
  column and state it in the methods paragraph.
- `multi_group_compare` in R wrote the Tukey and the Games-Howell P values that
  `rstatix` rounds to three decimals. The R engine now takes Tukey from
  `stats::TukeyHSD`, computes Games-Howell with the formula of `rstatix` and
  `pingouin`, and takes the omnibus tests from `stats::anova`, `stats::oneway.test`
  and `stats::kruskal.test`. The Python table gains the `n1`, `n2` and `estimate`
  columns of the R table, so the two tables now have the same columns.
- The Python `two_group_compare` reported Cohen's d without its sign, because
  `pingouin` reports the magnitude. On the `sleep` data R gave -0.832 and Python
  0.832. The Python engine now computes the signed d with the formula of
  `rstatix::cohens_d`, negative when the first group is lower.
- The Python `cox_forest` stopped about 2e-6 short of the estimate of
  `survival::coxph`, because `lifelines` uses a precision of 1e-7 by default. The
  fit now uses 1e-11, and the hazard ratios on the NCCTG lung data match R to
  about 1e-10. The Python terms now carry the R names, for example `sexmale` in
  place of `sex[male]`.
- The Python engine gave the points of the `bar` geom the colors of the wrong
  groups. `ax.bar` colors by the plotted order, which is alphabetical, while the
  seaborn call mapped the same palette to the groups in the order they appear in
  the file. The blue bar carried red points whenever the two orders differed.
  Every seaborn call now passes `hue_order`, so a group keeps one color in the
  figure and across the two engines.
- The R engine wrote `NA` in the effect size column of every Mann-Whitney run,
  although the methods paragraph stated that the effect size was reported.
  `rstatix::wilcox_effsize` needs the `coin` package, which the image does not
  hold, so the call failed and the guard wrote `NA`. The rank-biserial
  correlation is now computed from the U statistic in both engines, with the
  sign that says which group ranks lower.
- A review of the figures of the nine real-data cases found the defects below.
  Each one changes a figure, and some also change a stats table or a methods
  paragraph.
- The Python `heatmap` put white at the middle of the data range and not at zero.
  On the airway genes, white fell at z = 0.19. The figure also showed only the
  first annotation column, drew no dendrogram, hid part of the gene names and
  clustered with average linkage, while R uses complete linkage. The Python
  figure now uses the R colors and the R limits of -2 and 2 on the z-score, every
  annotation column with a legend, both dendrograms, complete linkage and every
  row label. When PyComplexHeatmap fails and seaborn draws the figure, the
  manifest records the reason in `fallback_reason`.
- The Python `correlation_heatmap` clustered with average linkage, which is the
  seaborn default, so the two engines could order the variables differently.
  Python now uses complete linkage and the R colors.
- `enrichment_dot` colored the dots by the adjusted p on a linear scale. The
  values span orders of magnitude, so almost every term had the same color. Both
  engines now color by -log10 of the adjusted p, and the methods paragraph says
  so. The Python stats table listed the terms in the plot order; it now lists
  them in the order of the adjusted p, as R does. Ties keep the row order in both
  engines.
- `volcano` labelled the legend and the y axis with "p" when the y column held an
  adjusted p, such as `padj`. The legend and the axis now name the column, for
  example "FC and padj" and "-log10 padj".
- The two engines drew different `pca` figures. R drew the multivariate t ellipse
  of `stat_ellipse`, and Python drew a chi-squared ellipse. Both engines now draw
  the normal-theory ellipse of `stat_ellipse(type = "norm")` for each group of
  four samples or more. R sorted the groups by the locale and Python kept the
  order of the file. Both engines now sort the groups by code point, so a group
  keeps one color.
- The `pca` PERMANOVA drew 999 random labelings even when the design has fewer
  distinct labelings, for example 70 for two groups of four. The test is now
  exact when there are 5000 labelings or fewer, and the methods paragraph names
  the form that ran. On the airway samples both engines give p = 2/70 = 0.029.
- The log-rank P of `survival_km` sat inside the plot panel, where a curve could
  cross it. Both engines now put it in a subtitle. The Python number-at-risk table
  also showed the censored and the event counts, which is three rows for each
  group where R shows one. It now shows the number at risk only.
- The Python `cox_forest` labelled its log axis in powers of ten, named each term
  by its model name, such as `celltypelarge`, and gave one star to every P below
  0.05. It now labels the axis with plain numbers such as 0.5 and 2, writes
  `celltype: large`, and grades the stars as survminer's `ggforest` does in R,
  up to three stars.
- The R `upset` sorted the intersections by degree and printed no count on the
  bars. It now sorts them by size and prints each count, as the Python engine
  does. The Python figure cut the set names on the left, because the bundle
  writer reset the figure to its starting size after upsetplot had widened it.
- `figkit inspect` printed one line for every column of a wide table, which was
  508 lines for the airway PCA matrix. The text report now shows the first 25
  columns and counts the rest, and each candidate list stops at 20 names.
  `--format json` still gives every column.
- `figkit build` reused the cached image when only `VERSION` changed. Podman keys
  the LABEL step on its text and not on the value of the build argument, so the
  new tag pointed at the old image with the old version label. The Containerfile
  now writes the version to `/opt/pubready/VERSION` in a RUN step, which keys on
  the value.

## [0.4.1] - 2026-07-19

### Changed
- The volcano example now uses the airway RNA-seq benchmark (Himes et al. 2014,
  GSE52778): real DESeq2 results for dexamethasone-treated vs untreated airway
  smooth muscle cells, with the donor (cell line) blocked in the model and genes
  labeled by HGNC symbol. `example/de_results.csv` and the volcano reference
  bundles were regenerated, thresholded on the adjusted p.
- Container adds DESeq2, airway, and org.Hs.eg.db (Bioconductor), used only to
  build that example (not a recipe dependency). Default image tag is
  `localhost/pubplot:0.4.1`.

## [0.4.0] - 2026-07-19

### Added
- `volcano` recipe (both engines): an EnhancedVolcano-style scatter of a
  differential-expression table. Auto-detects DESeq2 / limma / edgeR columns
  (override with `--x`/`--y`/`--label`), a four-color scheme on the `--fc_cutoff`
  and `--p_cutoff` thresholds, and top hits labeled with ggrepel / adjustText.
- `enrichment_dot` recipe (both engines): a clusterProfiler-style GSEA / ORA dot
  plot. An NES column gives a GSEA dot plot; otherwise ORA. Size is the gene
  count, color the adjusted p, `--top_n` terms shown.
- `paired_compare` recipe (both engines): a paired before/after plot connecting
  each subject's two points (`--id`), with a paired t-test or Wilcoxon signed-rank
  chosen from the differences' normality.
- `proportions` recipe (both engines): a chi-square test of independence, or
  Fisher's exact when a cell is sparse, with an odds ratio for a 2x2 table, drawn
  as a 100% stacked proportion bar.
- `pca` recipe (both engines): a PC1/PC2 scatter colored by `--group`, with 95%
  confidence ellipses, variance explained on the axes, and a seeded PERMANOVA p.
- `--geom raincloud` and `--geom bar` (mean + SEM + individual points) on the
  two-group and multi-group comparison recipes.
- New CLI flags `--id`, `--group`, `--label`, `--fc_cutoff`, `--p_cutoff`,
  `--top_n`; seeded example generators and data, reference bundles, smoke cases,
  and unit tests for the five new recipes and both geoms.

### Changed
- Container adds R `ggrepel` + `ggdist` and Python `adjustText`. Default image tag
  is `localhost/pubplot:0.4.0`.

## [0.3.0] - 2026-07-19

### Added
- `correlation` recipe (both engines): a two-variable scatter with a linear fit
  and 95% CI band. Pearson, Spearman, or Kendall chosen from the normality of
  each variable, or forced with `--test`; the coefficient and p are annotated on
  the plot.
- `correlation_heatmap` recipe (both engines): pairwise correlation among the
  numeric columns of a table, drawn as a clustered heatmap with dendrograms and
  per-cell coefficients with BH-adjusted significance stars. ComplexHeatmap in R,
  seaborn `clustermap` in Python.
- `upset` recipe (both engines): an UpSet plot of set intersections from a binary
  membership matrix, with a descriptive set-size and intersection-size table (no
  hypothesis test). ComplexHeatmap `UpSet` in R, the `upsetplot` package in Python.
- `--theme prism`: the GraphPad Prism look on the comparison recipes and the
  correlation scatter (ggprism `theme_prism` in R, a matched matplotlib style in
  Python). The previously inert `appearance.theme` spec key is now honored.
- Example datasets and seeded generators (`correlation_xy.csv`,
  `correlation_vars.csv`, `set_membership.csv`), reference bundles, smoke cases,
  and unit tests for the three new recipes.

### Changed
- Python engine adds the `upsetplot` dependency; `ggprism` is now declared in the
  R `DESCRIPTION` Imports. Default image tag is `localhost/pubplot:0.3.0`.

## [0.2.0] - 2026-07-18

### Added
- The R engine (`core/r`) is now a source R package: `DESCRIPTION`, `NAMESPACE`,
  code under `R/`, and a `testthat` suite under `tests/testthat`.
- The Python engine (`core/python`) is now a declared package: `pyproject.toml`
  with pinned dependencies, a single-sourced version, `uv.lock`, and a `pytest`
  suite under `tests`.
- `figkit test` runs both engine suites inside the pinned image against the
  mounted source, so the tested code is exactly what figkit runs.
- `figkit version` prints the pubplot version, the image tag, id, and digest,
  the Podman version, and the git commit.
- Root `VERSION` file as the single source of truth for the version.
- `CHANGELOG.md` (this file).
- Container image now carries an `org.opencontainers.image.version` label, set
  from the `VERSION` file at build time.

### Changed
- `core/r/lib` and `core/r/recipes` moved to `core/r/R` (package layout).
  `bootstrap.R` now sources `R/`; the figkit command surface is unchanged.
- Default image tag is `localhost/pubplot:0.2.0`.

### Notes
- Engine code is still loaded from the runtime mount, not installed into the
  image, so editing a recipe reruns with no rebuild. Packaging adds testing,
  versioning, and dependency declaration without coupling code to the image.
- Two Python lockfiles coexist by design: the image freezes runtime dependencies
  to `/opt/pubplot/requirements.lock` (the reproducibility pin, with the image
  digest), while `core/python/uv.lock` reproduces a development install of the
  package. `pytest` lives in the image but is not in the runtime lock.

## [0.1.0] - 2026-07-18

### Added
- First working version. All six recipe families in both the R and Python
  engines: two-group comparison, multi-group ANOVA, factorial ANOVA,
  Kaplan-Meier, Cox hazard-ratio forest, and clustered heatmap.
- `figkit` CLI (`inspect`, `plot`, `render`, `build`, `shell`) running inside a
  pinned Podman image (rocker/r-ver:4.4.1 + uv).
- Self-contained reproducibility bundle per figure: figure, runnable script,
  stats table, input copy, manifest, methods paragraph, and `REPRODUCE.md`.
- Host adapters for Claude Code, Codex, and opencode on one shared core.

[Unreleased]: https://github.com/wguesdon/pubready/compare/v3.0.0...HEAD
[3.0.0]: https://github.com/wguesdon/pubready/compare/v2.0.0...v3.0.0
[2.0.0]: https://github.com/wguesdon/pubready/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/wguesdon/pubready/compare/v0.4.1...v1.0.0
[0.4.1]: https://github.com/wguesdon/pubready/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/wguesdon/pubready/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/wguesdon/pubready/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/wguesdon/pubready/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/wguesdon/pubready/releases/tag/v0.1.0
