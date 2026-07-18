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
- Built the container (rocker/r-ver:4.4.1 + uv), image localhost/pubplot:0.1.0.
- Implemented the R engine: figkit CLI (inspect/plot/render/build/shell), the
  two_group_compare recipe, the assumption-driven test resolver, the house
  theme, and the full artifact bundle writer.
- Verified end to end on example data: t-test path, Mann-Whitney path, and a
  render round-trip from an edited config. Confirmed the emitted standalone
  script reproduces the figure in the pinned container.
- Added example test cases (tumor_volume.csv, cytokine_pg_ml.csv,
  edited_config.yaml) and tests/smoke_test.sh; all three cases pass.

Milestone status:
- M1 (container) DONE. M2 (inspect + two_group_compare in R + bundle) DONE.
- Next: M3 Python engine (same recipe), then M4 reference decision tree + Claude
  Code skill, M5 Codex/opencode adapters, M6 remaining recipes.

How to run (from repo root, image already built):
- ./cli/figkit inspect --data example/tumor_volume.csv
- ./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv --x group --y volume
- ./tests/smoke_test.sh
- Rebuild image: ./cli/figkit build

Resolved this session:
- R core: plain .R scripts sourced by the CLI, not a formal package.
- Figure theme: one house style for v1, journal presets later.
- Reproducibility artifact bundle: every figkit plot run writes one
  self-contained folder with the figure, runnable script, stats CSV, a copy of
  the input, manifest.json (provenance + container digest + package versions),
  session_info.txt, a manuscript-ready methods_<ts>.md, and REPRODUCE.md.
  See PRD.md "Reproducibility: the artifact bundle" for the manifest schema.
- Bundle also includes data_log.md (original filename + md5 + sha256 + column
  types + cleaning steps with row counts) and plot_config.yaml (editable spec:
  engine choice, palette, labels, axis limits, theme, bracket options). The
  agent co-edits plot_config.yaml with the user; figkit render --config
  regenerates the figure without redoing the chat, recomputing stats each time.

Environment checked:
- Podman 4.9.3, rootless, overlay driver. podman-compose present.

Still open:
- Excel ingestion cleaning threshold.
- CLI distribution (container entrypoint vs installed command shelling to Podman).
