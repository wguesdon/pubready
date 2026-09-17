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
| `upset` | binary membership matrix | `--data` | descriptive intersection sizes (no test) |
| `volcano` | DE table (log2FC + p) | `--data [--x --y --label]` | EnhancedVolcano-style scatter, top hits labeled |
| `enrichment_dot` | GSEA / ORA result table | `--data [--label]` | dot plot, size = count, color = adjusted p |
| `paired_compare` | condition, value, id | `--x --y --id` | paired t-test or Wilcoxon, connecting lines |
| `proportions` | two categoricals | `--x --y` | chi-square / Fisher, 100% stacked bar |
| `pca` | samples × features + group | `--data --group` | PC1/PC2 scatter, 95% ellipses, PERMANOVA |
| `spider_response` | patient, time, value | `--x --y --id [--group]` | descriptive change from baseline, RECIST lines (no test) |

Common options: `--engine r|python` (default r), `--geom box|violin|bar|raincloud`,
`--paired`, `--test auto|welch_t|wilcoxon|anova|kruskal|art` (correlation:
`pearson|spearman|kendall`; proportions: `chisq|fisher`), `--scale`/`--cluster`
(heatmaps), `--fc_cutoff`/`--p_cutoff`/`--top_n` (volcano), `--id --group --label`,
`--theme pubready_house|prism`, `--xlab --ylab --title --palette`, `--out DIR`.

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
test, and the smoke test asserts it. Three current gaps. The aligned rank transform
for non-parametric factorial ANOVA is R only, so `--test art --engine python` errors
and routes you to R. The Python heatmap keeps its per-feature significance stars in
the stats table rather than on the figure. In `multi_group_compare` the two engines
report different columns, and an extreme Tukey p underflows to 0 in R where Python
still prints the exponent.

## Read more

- `docs/how_it_works.md` — the full option table, an example for every recipe, and the bundle contents.
- `example/README.md` — worked commands and the committed reference figures.
- `reference/` — the statistical decision tree and assumption-check notes.

## Install

Codex reads the nearest `AGENTS.md` from the working directory up. To activate
this adapter for the repo, link it at the repo root:

```bash
ln -s skills/codex/AGENTS.md AGENTS.md
```

If the repo already has a root `AGENTS.md`, append the contents of this file to it
instead of overwriting.
