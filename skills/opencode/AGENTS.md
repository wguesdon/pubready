# pubplot (opencode adapter)

When the user has tidy data and wants a publication figure, use pubplot. It makes
a publication-ready figure with the statistical test drawn on top, the way a
scientist uses publication figure tools publication figures. A small library of tested, parameterized recipes
does the plotting and the statistics. Your job is to choose the right test with
the user and call the recipe. Do not write ggplot or matplotlib by hand. Write
original plotting or stats code only when no recipe fits, and say so when you do.

## The one command surface

Everything goes through `./cli/figkit`, which runs the pinned Podman image. Never
call R or Python directly.

```
figkit build                                  # once, builds the image
figkit inspect --data FILE [--format json]    # columns, types, group sizes, candidate x/y
figkit plot --recipe NAME --data FILE ...     # draw + test, writes a bundle
figkit render --config plot_config.yaml       # redraw from an edited spec, recompute stats
```

If a command reports the image is missing, run `figkit build` first.

## Workflow

1. Run `figkit inspect` on the data. Read the real column types and group sizes.
   Do not guess them.
2. Talk through the design with the user: how many groups, paired or not, one
   factor or several, time to event, or a matrix. Confirm before drawing.
3. Call `figkit plot` with the recipe and columns. Leave `--test auto` unless the
   user forces a test.
4. Point them at the bundle folder under `pubplot_output/`: the figure, the
   standalone script, the stats table, the input copy, and a methods paragraph.
5. To restyle, edit `plot_config.yaml` in the bundle and run `figkit render`. It
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

Common options: `--engine r|python` (default r), `--geom box|violin|bar`,
`--paired`, `--test auto|welch_t|wilcoxon|anova|kruskal|art`, `--scale`/`--cluster`
(heatmap), `--xlab --ylab --title --palette`, `--out DIR`.

## Choosing the test

`--test auto` measures before it picks. For the two-group and multi-group recipes
it checks normality per group with Shapiro-Wilk and equal variance with Levene,
then goes parametric when both hold, non-parametric otherwise, with the Welch
correction when variances differ. Walk this with the user and override any step
with `--test`. Event columns code `1 = event, 0 = censored`. `survival_km` needs
`--x` for the grouping arm.

## Engine limits

Pick the engine per figure with `--engine`. Both write the same bundle and, where
a matching method exists, the same numbers. Two current gaps. The aligned rank
transform for non-parametric factorial ANOVA is R only, so `--test art --engine
python` errors and routes you to R. The Python heatmap keeps its per-feature
significance stars in the stats table rather than on the figure.

## Read more

- `docs/how_it_works.md` — the full option table, an example for every recipe, and the bundle contents.
- `example/README.md` — worked commands and the committed reference figures.
- `reference/` — the statistical decision tree and assumption-check notes.
