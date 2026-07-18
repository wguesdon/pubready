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

[Unreleased]: https://github.com/wguesdon/pubplot/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/wguesdon/pubplot/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/wguesdon/pubplot/releases/tag/v0.1.0
