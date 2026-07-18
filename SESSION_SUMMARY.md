# Session summary

Resume point for work on pubplot. Update as work progresses.

## 2026-07-18

Kickoff session. Agreed the concept and drafted the PRD.

Decisions made:
- Repo name: pubplot. Private, on GitHub under wguesdon.
- Goal: an agentic skill that turns tidy CSV/XLSX into publication-ready figures
  with statistics drawn on top, used the way a scientist uses publication figure tools publication figures.
- Thin LLM. A curated library of parameterized functions does plotting and
  stats. The model chats to pick the test, then writes a short script calling
  those functions. It writes original code only when no recipe fits, and flags
  it when it does.
- Dual engine, user picks per figure: R (ggpubr + rstatix + ggprism) and Python
  (matplotlib/seaborn + pingouin + statannotations).
- Runs in a Podman container with pinned deps for reproducibility.
- Host tools: Claude Code, Codex, opencode, all from the start, on one shared
  core.
- Every figure ships with its script and its stats table for reproducibility.

Done this session:
- Created private repo, cloned locally.
- Wrote PRD.md.
- Scaffolded repo directories: core/, container/, cli/, skills/, reference/.

Next steps (not yet started, awaiting PRD sign-off before building):
- Resolve the open questions in PRD.md (figure theme presets, Excel cleaning,
  R packaging, CLI distribution).
- Milestone 1: container with both engines and pinned deps.
- Milestone 2: figkit inspect + two_group_compare recipe in R, end to end.

Open questions for the user:
- Figure theme: single house style or journal presets from the start?
- R core packaged as a real R package or scripts sourced by the CLI?
