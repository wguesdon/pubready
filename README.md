# pubplot

An agentic skill for making publication-ready figures from tidy data, the way a
scientist uses publication figure tools publication figures.

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

## Quick start

```bash
figkit build                                    # build the container image once
figkit inspect --data data.csv                  # see the columns and candidate x / y
figkit plot --recipe two_group_compare --data data.csv --x group --y value
```

Each run writes a self-contained bundle: the figure, the exact script, the stats,
a copy of the input, a manifest, and a manuscript-ready methods paragraph.

## Documentation

- [docs/how_it_works.md](docs/how_it_works.md) — the workflow, the recipes, the
  bundle, the two engines, and how to add a recipe.
- [example/README.md](example/README.md) — worked examples with R and Python
  figures for every recipe.
- [PRD.md](PRD.md) — the product plan.

## Status

All six recipes are built in both the R and Python engines: two-group comparison,
multi-group ANOVA, factorial ANOVA, Kaplan-Meier, Cox hazard ratios, and
clustered heatmaps. See [SESSION_SUMMARY.md](SESSION_SUMMARY.md) for progress and
what comes next (the agent skill wrappers).

## Layout

- `core/` — the recipe library, mirrored in R and Python.
- `container/` — the Podman image with both engines and pinned deps.
- `cli/` — the stable command surface the agent calls (`figkit`).
- `skills/` — thin adapters for Claude Code, Codex, and opencode.
- `reference/` — statistical test decision tree and assumption-check notes.
