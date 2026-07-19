# Session summary

Resume point for work on pubplot. Update as work progresses.

## 2026-07-19 — new recipes + Prism theme (v0.3.0)

Added three recipes and a style option, mirrored across both engines, taking the
catalog from six to nine.

- `correlation` (scatter + linear fit + CI, r and p annotated). Pearson /
  Spearman / Kendall chosen from per-variable normality, forced with `--test`.
  R: ggplot + `ggpubr::stat_cor`. Python: `pg.corr` + `sns.regplot`. Numbers match
  (r = 0.71, p = 2.09e-10 on the example).
- `correlation_heatmap` (clustered correlation matrix, coefficients + BH stars in
  each cell, dendrograms). Input is a table of numeric variables; correlation is
  among columns. R: ComplexHeatmap draw path. Python: `sns.clustermap`. No new dep.
- `upset` (UpSet plot from a binary membership matrix, descriptive intersection
  sizes, no test). R: `ComplexHeatmap::make_comb_mat` + `UpSet` (no new dep).
  Python: new `upsetplot` dependency (only reason for a rebuild).
- `--theme prism`: wired the previously inert `appearance.theme` spec key through
  both spec builders and both entry scripts. R adds `pub_theme()` in `theme.R`
  branching to `ggprism::theme_prism`; Python adds `apply_prism_style` in
  `theme.py`. Applied to the four ggplot comparison recipes + correlation scatter.

Plumbing: `--theme` option added to `entry/plot.{R,py}`; `ggprism` added to
DESCRIPTION Imports; `upsetplot` added to `requirements-extra.txt` + pyproject
(+ image rebuilt as `localhost/pubplot:0.3.0`). Bumped VERSION / DESCRIPTION /
version.{R,py} and the runtime image-tag fallbacks to 0.3.0.

Examples/tests/docs: seeded generators `make_correlation_data.R` +
`make_upset_data.R` → `correlation_xy.csv`, `correlation_vars.csv`,
`set_membership.csv`; reference bundles regenerated under `example/expected/`
(R + Python, plus a Prism showcase); `example/README.md` gallery updated; smoke
cases [10]-[13] + Python `pyc` lines; unit tests for all three recipes and a
theme spec test in both engines. Updated `cli/figkit`, `docs/how_it_works.md`,
the three skill adapters + the opencode command, PRD recipe catalog, CHANGELOG.

Verified: `figkit test` green (R + 20 Python), full `tests/smoke_test.sh` green
(13 R cases + all Python), every new recipe eyeballed on both engines.

Deferred: `proportions` recipe (chi-square / Fisher) remains the next candidate.
Python UpSet emits a harmless pandas FutureWarning from inside the upsetplot
library. cox_forest is not Prism-themed (survminer forest, low themeability).

## 2026-07-18 — packaging and versioning (v0.2.0)

Turned both engines into proper packages and added release plumbing. Key
decision: keep loading engine code from the runtime mount (fast iteration, code
pinned by git commit) rather than installing it into the image. Packaging is for
tests, versioning, and dependency declaration, not for baking code in.

- R engine (`core/r`) is now a source R package: `DESCRIPTION`, `NAMESPACE`,
  code moved from `lib/`+`recipes/` to `R/`, testthat suite in `tests/testthat`.
  `bootstrap.R` sources `R/`; figkit command surface unchanged.
- Python engine (`core/python`) has `pyproject.toml` (deps mirrored, version
  single-sourced from `version.py`), `uv.lock` (now committed; removed from
  `.gitignore`), and a pytest suite in `tests`.
- `figkit test` runs both suites in-container against the mounted source;
  `figkit version` prints version + image tag/digest + podman + commit.
- Root `VERSION` file is the single source of truth; bumped to 0.2.0 across R
  DESCRIPTION, python version.py, image tag, and runtime fallbacks. Added
  `CHANGELOG.md` (Keep a Changelog). Container gained a test-tooling layer
  (testthat + pytest, late so heavy layers stay cached) and a version label.
- Validated: `figkit test` green (R 4 files, Python 14 passed), full
  `tests/smoke_test.sh` green (all recipes, both engines). Image rebuilt as
  `localhost/pubplot:0.2.0` in ~28s (cache reuse).
- Open decision deferred to the user: license choice (leaning MIT). DESCRIPTION
  uses `License: file LICENSE` as a placeholder; no LICENSE file added yet.

## 2026-07-18

Kickoff session. Agreed the concept and drafted the PRD.

Decisions made:
- Repo name: pubplot. Private, on GitHub under wguesdon.
- Goal: an agentic skill that turns tidy CSV/XLSX into publication-ready figures
  with statistics drawn on top.
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
- heatmap BUILT: ComplexHeatmap clustered heatmap from a matrix CSV (first col =
  feature id) + optional annotation CSV (--annotation, --scale, --cluster). Row
  z-score + cluster both by default; per-feature Welch t-test with significance
  stars when the annotation has 2 groups. Added ComplexHeatmap + circlize to the
  image (rebuilt). Refactor: save_figure gained a draw-function path (base-
  graphics devices, not just ggsave) and write_bundle copies extra_inputs (the
  annotation) into the bundle. Example expression_matrix.csv +
  expression_annotation.csv (make_heatmap_data.R, seeded), reference, smoke [9].
- ALL SIX recipe families built in R AND Python: two_group, multi_group,
  factorial_anova, survival_km, cox_forest, heatmap.
- Python engine BUILT under core/python/pubplot (package) + core/python/entry.
  figkit dispatches on --engine r|python (default r). Python stack: pingouin +
  statannotations (two/multi-group), statsmodels (factorial, parametric only —
  no ART), lifelines (KM + Cox), PyComplexHeatmap w/ seaborn clustermap fallback.
  Numbers match R (KM p=0.007, Cox HR 0.545, ANOVA F identical). Added
  statsmodels/scikit-posthocs/lifelines/PyComplexHeatmap to the image (new layer
  after Bioconductor so R layers stay cached).
- Python bundles are byte-identical in structure to R (engine="python" in
  manifest, script_*.py, session_info = pip list). Reference bundles generated
  with _py stamps alongside R in example/expected/. Smoke test now 9 R + 6 Python.
- Gotchas fixed: entry/inspect.py shadowed stdlib inspect (renamed inspect_data.py);
  pingouin 0.6.1 uses underscore columns (p_val, cohen_d, U_val).
- Done: host adapters built — skills/claude-code/SKILL.md, skills/codex/AGENTS.md,
  skills/opencode/ (AGENTS.md + /pubplot command). Each carries the same figkit
  operating guide; detail stays in docs/how_it_works.md so they do not drift.
- Next: reference decision tree (reference/decision_tree.md, assumptions.md).
  Deferred: factorial pairwise post-hoc brackets; survival_km needs --x;
  Python heatmap has no per-feature stars on the figure (in stats CSV only).

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
