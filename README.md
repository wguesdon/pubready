# pubplot

An agentic skill for making publication-ready figures from tidy data, with the
right statistical test chosen and drawn on top.

Drop in a tidy CSV or Excel file, chat through the right statistical test, and
get a figure with the significance annotation drawn on top, plus the exact
script and numbers behind it.

## How it works

- **Thin LLM.** A curated library of parameterized functions does the plotting
  and the statistics. The agent chats to choose the correct test, then writes a
  short script that calls those functions. It writes original code only when no
  recipe covers the question, and flags it when it does.
- **Two engines, your choice.** R (`ggpubr` + `rstatix` + `ggprism`) or Python
  (`matplotlib`/`seaborn` + `pingouin` + `statannotations`), picked per figure.
- **Reproducible.** Everything runs in a Podman container with pinned versions.
  Every figure ships with its script and its stats table.
- **Runs where you work.** Claude Code, Codex, and opencode, on one shared core.

## Setup

Two steps, once.

1. Build the container image. Needs Podman.

```bash
./cli/figkit build
```

2. Install the skill in your agent.

| Host | Install (run from the repo root) |
|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/pubplot` |
| Codex | `ln -s skills/codex/AGENTS.md AGENTS.md` (append if a root `AGENTS.md` exists) |
| opencode | add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`; for the `/pubplot` command also `mkdir -p .opencode/command && ln -s ../../skills/opencode/command/pubplot.md .opencode/command/pubplot.md` |

Now open the agent in this folder and it can make figures.

## Talk to it

You do not name a recipe or a test, and you do not need the skill name. Drop a
tidy CSV or Excel file in the folder and describe what you want. The skill reads
the columns, confirms the design with you, runs the test, and hands back the
figure. These prompts use the files in `example/`, so you can try them as is.

- "I have `example/tumor_volume.csv` with a control and a treated group. Is the treatment significant? Give me a box plot with the p-value on it."
- "`example/gene_expression.csv` has four genotypes. Which differ in expression? I want the significant pairs bracketed."
- "`example/twoway_response.csv` is genotype crossed with treatment. Run the ANOVA and show the interaction."
- "`example/survival_trial.csv` has time, event, and arm. Draw a Kaplan-Meier with the log-rank p."
- "From `example/survival_trial.csv`, what is the hazard ratio for arm adjusting for age, sex, and stage? A forest plot please."
- "Cluster `example/expression_matrix.csv` and annotate the samples with `example/expression_annotation.csv`."

The shortest useful prompt is the filename plus what to compare: "plot
`example/gene_expression.csv`, expression by genotype, with stats." Mention the
data and a figure or a comparison and the skill takes over. In Claude Code you can
force it with `/pubplot`, but you rarely need to.

Each run writes a self-contained bundle: the figure, the exact script, the stats,
a copy of the input, a manifest, and a manuscript-ready methods paragraph. The
figure lands at `pubplot_output/<run>/figure_*.png`, plus PDF and SVG.

## Drive figkit directly

Prefer the command line, or scripting a batch? Call `figkit` yourself.

```bash
figkit inspect --data data.csv                  # see the columns and candidate x / y
figkit plot --recipe two_group_compare --data data.csv --x group --y value
```

## Documentation

- [docs/how_it_works.md](docs/how_it_works.md) — the workflow, the recipes, the
  bundle, the two engines, and how to add a recipe.
- [example/README.md](example/README.md) — worked examples with R and Python
  figures for every recipe.
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release.
- [PRD.md](PRD.md) — the product plan.

## Status

All six recipes are built in both the R and Python engines: two-group comparison,
multi-group ANOVA, factorial ANOVA, Kaplan-Meier, Cox hazard ratios, and
clustered heatmaps. Host adapters ship for Claude Code, Codex, and opencode. Both
engines are packaged with test suites (`figkit test`) and versioned in lockstep
with the image, tracked in [CHANGELOG.md](CHANGELOG.md). See
[SESSION_SUMMARY.md](SESSION_SUMMARY.md) for progress and what comes next.

## Layout

- `core/` — the recipe library, mirrored in R and Python. `core/r` is a source
  R package and `core/python` a `pyproject.toml` package; both carry test suites,
  run with `figkit test`.
- `container/` — the Podman image with both engines and pinned deps.
- `cli/` — the stable command surface the agent calls (`figkit`).
- `skills/` — thin adapters for Claude Code, Codex, and opencode.
- `reference/` — statistical test decision tree and assumption-check notes.
