# Changelog

All notable changes to pubplot are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because reproducibility is the point of this tool, "breaking" here means either a
change to the `figkit` command surface or a change that alters figure output
(a default test, a post-hoc rule, or a plotting detail). Those bump the major
version; new recipes bump the minor; fixes that leave outputs unchanged bump the
patch. The repo version, both engine package versions, and the container image
tag move together.

## [Unreleased]

### Changed
- Moved the example data generators and `generate_expected.sh` out of `example/`
  into `tools/make_examples/`. `example/` now holds only the input datasets, the
  `expected/` reference bundles, and its README, so the folder a reader opens
  first shows examples, not maintainer tooling.
- The committed volcano example uses a small synthetic DE table (1200 genes) in
  place of the 928 KB airway DESeq2 result. The airway benchmark stays available
  as the optional `tools/make_examples/make_volcano_airway.R`. Committed footprint
  drops by roughly 5 MB and the volcano figure is unchanged in kind. Output is not
  affected for any other recipe.

### Added
- `tools/make_examples/check_determinism.sh`: regenerates the reference bundles in
  a temp tree and confirms every stats table matches the committed copy, so "same
  input, same numbers" is checkable rather than assumed.

## [0.4.1] - 2026-07-19

### Changed
- The volcano example now uses the airway RNA-seq benchmark (Himes et al. 2014,
  GSE52778): real DESeq2 results for dexamethasone-treated vs untreated airway
  smooth muscle cells, with the donor (cell line) blocked in the model and genes
  labeled by HGNC symbol. `example/de_results.csv` and the volcano reference
  bundles were regenerated, thresholded on the adjusted p.
- Container adds DESeq2, airway, and org.Hs.eg.db (Bioconductor), used only to
  build that example (not a recipe dependency). Default image tag is
  `localhost/pubplot:0.4.1`.

## [0.4.0] - 2026-07-19

### Added
- `volcano` recipe (both engines): an EnhancedVolcano-style scatter of a
  differential-expression table. Auto-detects DESeq2 / limma / edgeR columns
  (override with `--x`/`--y`/`--label`), a four-color scheme on the `--fc_cutoff`
  and `--p_cutoff` thresholds, and top hits labeled with ggrepel / adjustText.
- `enrichment_dot` recipe (both engines): a clusterProfiler-style GSEA / ORA dot
  plot. An NES column gives a GSEA dot plot; otherwise ORA. Size is the gene
  count, color the adjusted p, `--top_n` terms shown.
- `paired_compare` recipe (both engines): a paired before/after plot connecting
  each subject's two points (`--id`), with a paired t-test or Wilcoxon signed-rank
  chosen from the differences' normality.
- `proportions` recipe (both engines): a chi-square test of independence, or
  Fisher's exact when a cell is sparse, with an odds ratio for a 2x2 table, drawn
  as a 100% stacked proportion bar.
- `pca` recipe (both engines): a PC1/PC2 scatter colored by `--group`, with 95%
  confidence ellipses, variance explained on the axes, and a seeded PERMANOVA p.
- `--geom raincloud` and `--geom bar` (mean + SEM + individual points) on the
  two-group and multi-group comparison recipes.
- New CLI flags `--id`, `--group`, `--label`, `--fc_cutoff`, `--p_cutoff`,
  `--top_n`; seeded example generators and data, reference bundles, smoke cases,
  and unit tests for the five new recipes and both geoms.

### Changed
- Container adds R `ggrepel` + `ggdist` and Python `adjustText`. Default image tag
  is `localhost/pubplot:0.4.0`.

## [0.3.0] - 2026-07-19

### Added
- `correlation` recipe (both engines): a two-variable scatter with a linear fit
  and 95% CI band. Pearson, Spearman, or Kendall chosen from the normality of
  each variable, or forced with `--test`; the coefficient and p are annotated on
  the plot.
- `correlation_heatmap` recipe (both engines): pairwise correlation among the
  numeric columns of a table, drawn as a clustered heatmap with dendrograms and
  per-cell coefficients with BH-adjusted significance stars. ComplexHeatmap in R,
  seaborn `clustermap` in Python.
- `upset` recipe (both engines): an UpSet plot of set intersections from a binary
  membership matrix, with a descriptive set-size and intersection-size table (no
  hypothesis test). ComplexHeatmap `UpSet` in R, the `upsetplot` package in Python.
- `--theme prism`: the GraphPad Prism look on the comparison recipes and the
  correlation scatter (ggprism `theme_prism` in R, a matched matplotlib style in
  Python). The previously inert `appearance.theme` spec key is now honored.
- Example datasets and seeded generators (`correlation_xy.csv`,
  `correlation_vars.csv`, `set_membership.csv`), reference bundles, smoke cases,
  and unit tests for the three new recipes.

### Changed
- Python engine adds the `upsetplot` dependency; `ggprism` is now declared in the
  R `DESCRIPTION` Imports. Default image tag is `localhost/pubplot:0.3.0`.

## [0.2.0] - 2026-07-18

### Added
- The R engine (`core/r`) is now a source R package: `DESCRIPTION`, `NAMESPACE`,
  code under `R/`, and a `testthat` suite under `tests/testthat`.
- The Python engine (`core/python`) is now a declared package: `pyproject.toml`
  with pinned dependencies, a single-sourced version, `uv.lock`, and a `pytest`
  suite under `tests`.
- `figkit test` runs both engine suites inside the pinned image against the
  mounted source, so the tested code is exactly what figkit runs.
- `figkit version` prints the pubplot version, the image tag, id, and digest,
  the Podman version, and the git commit.
- Root `VERSION` file as the single source of truth for the version.
- `CHANGELOG.md` (this file).
- Container image now carries an `org.opencontainers.image.version` label, set
  from the `VERSION` file at build time.

### Changed
- `core/r/lib` and `core/r/recipes` moved to `core/r/R` (package layout).
  `bootstrap.R` now sources `R/`; the figkit command surface is unchanged.
- Default image tag is `localhost/pubplot:0.2.0`.

### Notes
- Engine code is still loaded from the runtime mount, not installed into the
  image, so editing a recipe reruns with no rebuild. Packaging adds testing,
  versioning, and dependency declaration without coupling code to the image.
- Two Python lockfiles coexist by design: the image freezes runtime dependencies
  to `/opt/pubplot/requirements.lock` (the reproducibility pin, with the image
  digest), while `core/python/uv.lock` reproduces a development install of the
  package. `pytest` lives in the image but is not in the runtime lock.

## [0.1.0] - 2026-07-18

### Added
- First working version. All six recipe families in both the R and Python
  engines: two-group comparison, multi-group ANOVA, factorial ANOVA,
  Kaplan-Meier, Cox hazard-ratio forest, and clustered heatmap.
- `figkit` CLI (`inspect`, `plot`, `render`, `build`, `shell`) running inside a
  pinned Podman image (rocker/r-ver:4.4.1 + uv).
- Self-contained reproducibility bundle per figure: figure, runnable script,
  stats table, input copy, manifest, methods paragraph, and `REPRODUCE.md`.
- Host adapters for Claude Code, Codex, and opencode on one shared core.

[Unreleased]: https://github.com/wguesdon/pubplot/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/wguesdon/pubplot/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/wguesdon/pubplot/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/wguesdon/pubplot/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/wguesdon/pubplot/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/wguesdon/pubplot/releases/tag/v0.1.0
