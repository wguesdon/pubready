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

## Status

Early development. See [PRD.md](PRD.md) for the plan and
[SESSION_SUMMARY.md](SESSION_SUMMARY.md) for current progress.

## Layout

- `core/` — the recipe library, mirrored in R and Python.
- `container/` — the Podman image with both engines and pinned deps.
- `cli/` — the stable command surface the agent calls (`figkit`).
- `skills/` — thin adapters for Claude Code, Codex, and opencode.
- `reference/` — statistical test decision tree and assumption-check notes.
