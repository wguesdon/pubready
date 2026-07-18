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
- Recipe library expansion (from user's test list): ANOVA family first, then
  survival (KM + Cox + hazard-ratio forest), then heatmap (ComplexHeatmap /
  PyComplexHeatmap, matrix input). Post-hoc default: significant pairs only.
- multi_group_compare BUILT: one-way ANOVA / Welch / Kruskal-Wallis from the
  assumptions, Tukey / Games-Howell / Dunn post hoc, brackets for significant
  pairs only, omnibus test shown as a subtitle. Example gene_expression.csv +
  reference bundle + smoke case added.
- Refactor: recipes now return methods + build_script; write_bundle no longer
  hard-codes two-group logic, so new recipes plug in cleanly.
- factorial_anova BUILT: two-way and three-way ANOVA (Type II), aligned rank
  transform (ARTool) for the non-parametric factorial path, grouped box/bar,
  faceted by the third factor, effects shown as a wrapped subtitle. Added ARTool
  + emmeans to the image (image rebuilt). Added --fill / --facet CLI options and
  per-recipe figure sizing (recipes return width/height). Examples
  twoway_response.csv + threeway_response.csv (from make_factorial_data.R,
  seeded), reference bundles, and smoke cases [5]/[6].
- ANOVA family complete (one-way + factorial).
- Survival family BUILT: survival_km (Kaplan-Meier + censoring ticks + log-rank
  p on plot + number-at-risk table via survminer::ggsurvplot, combined with
  ggpubr::ggarrange) and cox_forest (Cox PH + survminer::ggforest hazard-ratio
  forest + cox.zph check + broom HR table). Added survival + survminer to the
  image (rebuilt). New CLI options --time/--event/--covariates. Example
  survival_trial.csv (make_survival_data.R, seeded), reference bundles, smoke
  cases [7]/[8].
- Refactor: cleaning moved into recipes (each returns df_used + clean_steps);
  recipes can supply a bundle label (survival has no x/y) and figure width/height.
  Use survminer::surv_fit (not survfit) to avoid the "symbol not subsettable"
  ggsurvplot error with a variable formula.
- Next per user's list: heatmap (ComplexHeatmap / PyComplexHeatmap, matrix input,
  Bioconductor rebuild). After recipes: Python engine, reference decision tree +
  Claude Code skill, Codex/opencode adapters.
- Deferred: factorial pairwise post-hoc brackets (emmeans/ART contrasts, pkgs
  installed); survival_km currently requires a grouping column (--x).

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
