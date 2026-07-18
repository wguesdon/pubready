# pubplot — Product Requirements Document

Status: draft v0.1
Owner: wguesdon
Last updated: 2026-07-18

## Problem

Bench scientists make figures for papers in tools like publication figure tools publication figures. They
paste in tidy data, click through a statistical test, and get a plot with the
significance annotation drawn on top. The workflow is fast, but it is a closed
tool: no reproducible script, no version control, and choices about which test
to run are made through menus rather than reasoning.

pubplot gives the same experience through a chat agent. The scientist drops in
tidy data, talks through the right test, and gets a publication-ready figure
plus the exact script and numbers behind it.

## Users

Primary: a wet-lab or clinical scientist comfortable with Excel and tidy data,
not necessarily fluent in R or Python. They know their experimental design.
They want the correct test and a clean figure, not a plotting tutorial.

Secondary: a computational collaborator who wants the emitted script to be
readable, correct, and reproducible so it can go into a methods section.

## Core principle: thin LLM

The language model does as little code authoring as possible. A curated library
of parameterized functions does the plotting and the statistics. The model's job
is to converse to choose the correct test, then write a short script that calls
those functions with the right arguments.

The model writes original analysis code only when no recipe covers the question.
When it does, that code is clearly flagged in the output so the scientist knows
it left the validated path and should review it.

## Dual engine

Every recipe exists in two backends. The scientist picks per figure.

- R: `ggpubr` + `rstatix` + `ggprism`. The strongest path for significance
  brackets and p-values positioned above groups, and the closest match to the
  publication figures look.
- Python: `matplotlib` / `seaborn` + `pingouin` + `statannotations`. For users
  who prefer a Python toolchain.

Default suggestion is R for the figure quality, but the choice is always the
user's.

## Workflow

1. The scientist drops a tidy CSV or XLSX into a project folder.
2. A fixed `inspect` command reports columns, types, group sizes, and quick
   normality and variance checks. The data shape is measured, not guessed.
3. The agent chats through the experimental design: paired or unpaired, number
   of groups, parametric assumptions, one or two factors, repeated measures. A
   decision-tree reference doc guides the questions.
4. The agent recommends a test, explains why in one or two lines, and the user
   confirms or overrides.
5. The matching recipe renders inside a Podman container.
6. Output: one self-contained artifact bundle folder (see below) holding the
   figure, the script, the stats, a copy of the input, the full environment, a
   manifest, and a manuscript-ready methods paragraph.

## Statistics on top

This is the feature that makes it feel like publication figures.

- R engine: `rstatix` computes the test and adjusted p-values, `ggpubr`
  (`stat_pvalue_manual`) and `ggprism` place the brackets and labels above the
  relevant groups.
- Python engine: `statannotations` draws the same brackets on matplotlib or
  seaborn axes, with statistics from `pingouin` or `scipy`.

## Repository layout

- `core/` — the recipe library, mirrored across R and Python. Recipes such as
  `two_group_compare`, `multi_group_anova`, `paired_compare`, `correlation`,
  and grouped box / violin / bar with brackets.
- `container/` — one Podman image with R and Python and pinned versions, so a
  figure rendered today renders the same next year.
- `cli/` — the stable command surface the agent calls. Example:
  `figkit inspect data.csv`,
  `figkit plot --engine r --recipe two_group_compare --x group --y value --test auto`.
- `skills/` — thin adapters that point at the same CLI and reference docs, one
  per host tool: Claude Code `SKILL.md`, Codex `AGENTS.md`, opencode.
- `reference/` — the statistical test decision tree and assumption-check notes,
  shared by every host tool.

## Reproducibility: the artifact bundle

Reproducibility is the trust anchor, and it is enforced by structure, not by
convention. Every run of `figkit plot` writes one self-contained output folder.
That folder is the artifact. It holds everything needed to regenerate the exact
figure on another machine years later.

Bundle folder name is snake_case with a timestamp, for example
`two_group_compare_tumor_volume_2026_07_18_140355/`.

Contents of every bundle:

- `figure_<name>_<timestamp>.pdf`, `.png`, `.svg` — the figure in each format.
- `script_<timestamp>.R` or `.py` — the exact, runnable script that produced the
  figure. Standalone: reads the input, runs the test, draws the plot.
- `stats_<timestamp>.csv` — the statistics table (test, statistic, p-value,
  adjusted p-value, effect size, group n).
- `input_<name>.csv` — a copy of the exact input data used, so the bundle does
  not depend on an external file that may change.
- `manifest_<timestamp>.json` — structured provenance (schema below).
- `session_info.txt` — full environment dump: R `sessionInfo()` or Python
  `pip freeze`. The complete package set, not just the headline ones.
- `methods_<timestamp>.md` — a manuscript-ready methods paragraph (see below).
- `REPRODUCE.md` — the one command to rerun this bundle in the pinned container.

### Methods paragraph

Every bundle includes a methods file written in the register of a journal
methods section, ready to paste into a manuscript. It is templated from the same
manifest that drives the figure, so the reported test, software, versions,
correction method, and sample sizes always match what was actually run. No
hand-copying numbers between the figure and the paper.

Example output:

> Tumor volume was compared between the two groups with a Welch two-sample
> t-test. Normality was assessed by the Shapiro-Wilk test and equal variance was
> not assumed. Effect size is reported as Cohen's d. Significance was set at
> P < 0.05. Analyses were performed in R 4.4.1 with rstatix 0.7.2, and figures
> were produced with ggpubr 0.6.0 and ggprism 1.0.5. Group sizes were n = 21 and
> n = 21.

The paragraph is generated deterministically by the recipe from the manifest.
The agent may refine the wording on request, but the reported facts come from the
run, not from the model.

### Manifest schema

```json
{
  "pubplot_version": "0.1.0",
  "pubplot_git_commit": "<sha>",
  "created_utc": "2026-07-18T14:03:55Z",
  "engine": "r",
  "recipe": "two_group_compare",
  "recipe_version": "1",
  "arguments": { "x": "group", "y": "volume", "test": "auto" },
  "input": { "file": "input_tumor_volume.csv", "sha256": "<hash>", "n_rows": 42 },
  "statistical_test": {
    "name": "Welch two-sample t-test",
    "statistic": 3.14,
    "p_value": 0.004,
    "adjustment": "none",
    "effect_size": { "name": "Cohen's d", "value": 0.98 }
  },
  "container": {
    "image": "localhost/pubplot:0.1.0",
    "digest": "sha256:<digest>",
    "podman_version": "4.9.3"
  },
  "environment": {
    "language": "R",
    "language_version": "4.4.1",
    "packages": { "ggpubr": "0.6.0", "rstatix": "0.7.2", "ggprism": "1.0.5" }
  },
  "outputs": ["figure_...pdf", "figure_...png", "figure_...svg", "stats_...csv"]
}
```

The container is pinned by digest, not just tag, so `REPRODUCE.md` pulls the
identical image. The input is copied in and checksummed. The full package set is
captured. That combination is what makes the bundle reproducible rather than
merely re-runnable.

## Recipe catalog (initial target)

Each recipe is a named function with a fixed argument surface, present in both
engines.

- `two_group_compare` — t-test or Mann-Whitney, unpaired. Box or violin with
  bracket.
- `paired_compare` — paired t-test or Wilcoxon signed-rank. Before/after plot.
- `multi_group_anova` — one-way ANOVA or Kruskal-Wallis plus post-hoc, brackets
  for chosen comparisons.
- `two_factor` — two-way ANOVA, grouped bar or box.
- `correlation` — Pearson or Spearman, scatter with fit and CI.
- `proportions` — chi-square or Fisher, bar of proportions.

Survival, dose-response, and mixed-effects models are candidates for later
versions and are out of scope for v1.

## Test selection logic

The agent does not pick a test from a single prompt. It walks a decision tree in
`reference/`:

1. Outcome type: continuous, count, proportion, time-to-event.
2. Number of groups and whether they are paired or repeated.
3. Parametric assumptions, checked against the `inspect` output rather than
   assumed. Normality and equal-variance results steer parametric versus
   non-parametric.
4. Multiple-comparison correction when there is more than one comparison.

The recommendation names the test, the assumption checks that led to it, and the
correction method. The user can override any step.

## Host tool support

All three from the start, built on one shared core.

- Claude Code: `skills/claude-code/SKILL.md`.
- Codex: `skills/codex/AGENTS.md`.
- opencode: `skills/opencode/`.

Each adapter is thin. It tells the host tool how to call `figkit` and where the
reference docs live. The reasoning and the code live in the shared core so the
three stay in sync.

## Out of scope for v1

- A graphical user interface.
- Cloud hosting or a hosted service.
- Very large datasets (streaming, out-of-core).
- Exotic models beyond the documented escape hatch: mixed-effects, Bayesian,
  multivariate.

## Resolved decisions

- R core: plain `.R` scripts sourced by the CLI, not a formal R package. Easy for
  the agent to read and extend. Revisit packaging once several recipes exist.
- Figure theme: one clean house style for v1. Journal presets (Nature, Cell,
  PLOS) come later.

## Open questions

- Excel ingestion: how much cleaning to do for messy sheets before declaring the
  data not tidy and asking the user to fix it.
- How the CLI is distributed to users: a shell entrypoint into the container, or
  an installed command that shells out to Podman.

## Milestones

1. Container with R and Python engines and pinned deps.
2. `figkit inspect`, one end-to-end recipe (`two_group_compare`) in R, and the
   artifact bundle writer (figure, script, stats, input copy, manifest, session
   info, methods paragraph, REPRODUCE.md).
3. The same recipe in Python.
4. Decision-tree reference and the Claude Code skill wired to the CLI.
5. Codex and opencode adapters.
6. Remaining recipes.
