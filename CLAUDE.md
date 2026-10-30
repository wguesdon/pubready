# Repository instructions

For a publication figure request, read and follow `skills/claude-code/SKILL.md`
before you work. It defines the pubready workflow and it requires the
`./cli/figkit` command surface. The same guide is available as a project skill
if `.claude/skills/pubready` is linked, but read the file directly when it is not.

For work on the repository itself, follow the rules below.

## Repository layout

| Path | Contents |
|---|---|
| `cli/figkit` | The one command surface. It runs the pinned Podman image. |
| `pyproject.toml`, `src/pubready_cli/` | The uv package. Its `figkit` launcher runs `cli/figkit` from the files that the wheel carries. |
| `tests/install_test.sh` | The end-to-end test of a uv tool install of HEAD. |
| `tests/doctor_test.sh` | The test of `figkit doctor` and of the Podman errors, with stub `podman` and `uname` scripts. |
| `tools/render_audit/` | The stress datasets and the runner that renders every recipe with both engines for a visual review. |
| `core/r/` | R recipes, built on ggpubr, rstatix and ggprism. |
| `core/python/` | Python recipes, built on matplotlib, seaborn and pingouin. |
| `core/recipe_options.txt` | The options that each recipe reads. Both engines stop on any other option. |
| `container/` | The `Containerfile` and the pinned R and Python dependencies. |
| `skills/claude-code/SKILL.md` | The Claude Code adapter. |
| `skills/codex/AGENTS.md` | The Codex adapter, linked from the root `AGENTS.md`. |
| `skills/opencode/` | The opencode adapter and the `/pubready` command. |
| `reference/` | The statistical decision tree and the assumption-check notes. |
| `reference/decision_tree.yaml` | The versioned tree that `figkit advise` walks. |
| `core/python/pubready/advise.py`, `core/python/entry/advise.py` | The shared walker and its command entry. |
| `tests/advise_test.sh` | The end-to-end advice and emitted-command checks. |
| `docs/how_it_works.md` | The full option table, one example per recipe, and the bundle contents. |
| `example/` | The committed input data and the reference figures. |
| `tests/smoke_test.sh` | The end-to-end test over the example cases. |
| `benchmarks/` | Public benchmark datasets and the reference values that each recipe must reproduce. |
| `tests/benchmark_test.sh` | The end-to-end run of every benchmark case with both engines. |
| `tools/make_benchmarks/` | The exporter that writes the benchmark datasets and computes the reference values. |

## Rules for work in this repository

1. Never call R or Python directly, and never call `podman` directly. Every run
   goes through `./cli/figkit`. If a command reports that the image is missing,
   run `./cli/figkit build` first. If it reports a Podman problem, run
   `./cli/figkit doctor`.
2. Keep the two engines at parity. When you change a recipe in `core/r/`, make
   the matching change in `core/python/`, and keep the argument names and the
   reported numbers the same. Record a gap in the "Engine limits" section of each
   adapter when parity is not possible.
3. Update the three adapters together. `skills/claude-code/SKILL.md`,
   `skills/codex/AGENTS.md` and `skills/opencode/AGENTS.md` carry the same
   operating guide, so a change to the `figkit` surface goes into all three in
   one commit. Put the detail in `docs/how_it_works.md`, not in the adapters.
4. Run `./tests/smoke_test.sh` after a change to `cli/figkit`, to a recipe or to
   the container. It asserts that each example writes a complete bundle and that
   the expected test is chosen. Run `./cli/figkit test` too. It compares every
   recipe, in both engines, with the reference values in `benchmarks/`. After a
   change to `cli/figkit`, run `./tests/doctor_test.sh`. It needs no image. After a
   change to `cli/figkit`, `pyproject.toml` or `src/pubready_cli/`, commit and
   run `./tests/install_test.sh`. It installs HEAD with uv and draws a figure
   from the installed `figkit`. A new top-level folder that `cli/figkit` mounts
   must also go into the `force-include` table of `pyproject.toml`.
   Run `./tests/advise_test.sh` after a change to either decision tree, the
   advise code, accepted test values or `cli/figkit`.
5. Add a new recipe to the recipe table in all three adapters, to
   `docs/how_it_works.md`, to `example/README.md` and to the smoke test. Add a
   line with its options to `core/recipe_options.txt`. Add a public benchmark
   case for it as `benchmarks/README.md` describes.
6. Write output only under `pubready_output/` or a directory that `--out` names.
   These directories are ignored by git.
7. Record a user-visible change in `CHANGELOG.md` and bump `VERSION` when the
   release goes out.

## Style

Use Simplified Technical English in the Markdown documentation, the adapters and
the reference notes. Keep the normal style of the language in code comments, in
roxygen blocks and in Python docstrings. Python docstrings use the Google style.
