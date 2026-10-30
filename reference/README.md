# reference

Shared knowledge for choosing a statistical test, used by every host tool. The
agent walks these during the chat rather than picking a test from a single prompt.

- [decision_tree.md](decision_tree.md) — by outcome type, group count, and
  pairing: the test each recipe reaches under `--test auto`, the flag to override
  it, and how multiple comparisons are corrected.
- [assumptions.md](assumptions.md) — the normality (Shapiro-Wilk) and
  equal-variance (Levene) checks the recipes run, their thresholds, and
  how each result steers the parametric or nonparametric choice.

Every rule in these files matches the recipe code, so the reference and the
analysis stay in step.

`decision_tree.yaml` is the machine-readable tree for `figkit advise`.
The prose tree and YAML are edited together. Tests check the recipes, headings,
canonical test values and recipe option lists.
