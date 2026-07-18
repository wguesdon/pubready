# How pubplot works

pubplot turns tidy data into a publication-ready figure the way a scientist uses
publication figure tools publication figures: drop the data in, decide on the right statistical test, and get a
figure with the significance annotation drawn on top. The difference is that
pubplot runs as an agent skill and ships a full reproducible record with every
figure.

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
- **Recipes** (`core/r/recipes`, `core/python/pubplot/recipes`) are the fixed
  functions. Each one owns its test choice, its figure, its stats table, its
  methods paragraph, and the standalone script it emits.
- **The container** (`container/Containerfile`) is one Podman image with R and
  Python and pinned versions. The environment is pinned by the image digest; the
  code is pinned by the pubplot git commit. Both are recorded in every bundle.
- **The bundle writer** takes what a recipe returns and writes the artifact
  folder. It is engine agnostic.

The agent stays thin because the recipes are thick. When a recipe covers the
question, the agent fills in arguments. When no recipe fits, the agent may write
original code, and that code is flagged so the scientist knows it left the tested
path.

## Using figkit

```
figkit inspect --data FILE [--format text|json]
figkit plot --recipe NAME --data FILE [options]
figkit render --config plot_config.yaml [--data FILE]
figkit build          # build or rebuild the container image
figkit shell          # open a shell in the image
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
| `--geom box\|violin\|bar` | figure geometry |
| `--scale`, `--cluster` | heatmap scaling and clustering |
| `--xlab`, `--ylab`, `--title`, `--palette` | labels and colors |
| `--out DIR` | output root (default `pubplot_output`) |

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

### Test selection

`--test auto` walks a short decision tree rather than picking a test from a single
prompt. For the two-group and multi-group recipes it checks normality per group
with the Shapiro-Wilk test and equal variance with Levene's test, then chooses:
parametric when both hold, non-parametric otherwise, with the Welch correction
when variances differ. The scientist can override any step by passing `--test`
explicitly.

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
```

## The artifact bundle

Every `figkit plot` run writes one folder. That folder is the reproducible unit.

- `figure_*.pdf`, `.png`, `.svg` — the figure in each format.
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
- `manifest_*.json` — structured provenance: recipe and arguments, the container
  pinned by digest, the pubplot git commit, and every package version.
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

Where Python has no faithful equivalent, pubplot routes you to R rather than
running something different. The one current case is the aligned rank transform
for non-parametric factorial ANOVA: `--test art --engine python` errors and points
to `--engine r`. The Python `heatmap` draws the clustered matrix and the group
bar but keeps the per-feature significance stars in the stats table rather than on
the figure.

## Setup

pubplot needs Podman and one built image.

```bash
figkit build          # builds localhost/pubplot:0.1.0 from container/Containerfile
```

The image carries R (ggpubr, rstatix, ggprism, survival, survminer, ARTool,
ComplexHeatmap) and Python (matplotlib, seaborn, statannotations, pingouin,
statsmodels, scikit-posthocs, lifelines, PyComplexHeatmap), all pinned. The
`core/` scripts are mounted at runtime, so editing a recipe does not need a
rebuild; only changing the dependencies does.

## Adding a recipe

A recipe is one function, `recipe_<name>`, in `core/r/recipes/<name>.R` and
`core/python/pubplot/recipes/<name>.py`. It receives the raw data frame and the
spec, and returns:

- the figure (a ggplot in R, a matplotlib figure in Python),
- a tidy stats table,
- a test-metadata dict for the manifest,
- a methods paragraph,
- a builder for the standalone script,
- the cleaned data it used and the cleaning steps.

The bundle writer does the rest. Because recipes own their own methods text and
script, new ones plug in without touching the writer. Register the recipe by
adding an example dataset, a smoke case in `tests/smoke_test.sh`, and a reference
bundle in `example/generate_expected.sh`.

## Design principles

- Thin agent, thick recipes. The model reasons and fills arguments; fixed code
  does the analysis.
- Measure the data before choosing a test.
- Put the statistic on the figure, not just in a table.
- Ship the record with the result. Every figure is reproducible or it does not
  ship.
- When an engine cannot do something faithfully, say so and route to the one that
  can.
