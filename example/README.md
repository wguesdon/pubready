# Example data and test cases

These files double as documentation and as the smoke-test inputs run by
`tests/smoke_test.sh`. Each one exercises a different path through the engine.

## Datasets

### `tumor_volume.csv` — parametric path
Two groups (control, treated), one continuous outcome, both groups roughly
normal with similar variance. The `auto` test resolver picks a **Student's
two-sample t-test**.

```bash
figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)"
```

### `cytokine_pg_ml.csv` — non-parametric path
Two groups with right-skewed values (each has a high outlier). The Shapiro-Wilk
check fails, so the `auto` resolver falls back to a **Mann-Whitney U test**.

```bash
figkit plot --recipe two_group_compare --data example/cytokine_pg_ml.csv \
  --x group --y concentration --ylab "IL-6 (pg/mL)"
```

### `gene_expression.csv` — multi-group (ANOVA)
Four groups (`wt`, `het`, `ko`, `rescue`), one continuous outcome, roughly normal
with similar variance. The `auto` resolver picks a **one-way ANOVA** with
**Tukey's HSD** post hoc. Only pairwise comparisons that reach significance are
bracketed; here `wt` vs `rescue` is not significant and is left off.

```bash
figkit plot --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)"
```

## Config

### `edited_config.yaml` — render path
A hand-edited `plot_config.yaml` that switches the geometry to a violin, sets a
custom palette, and uses a numeric p-value label. Regenerate the figure without
rerunning the whole conversation:

```bash
figkit render --config example/edited_config.yaml
```

This file deliberately uses an unquoted `y:` key to cover the YAML 1.1 quirk
where a bare `y` is read as a boolean. The reader normalizes it.

## Expected results

| Input                 | Test chosen                 | Significant            |
|-----------------------|-----------------------------|------------------------|
| `tumor_volume.csv`    | Student's two-sample t-test | yes                    |
| `cytokine_pg_ml.csv`  | Mann-Whitney U test         | yes                    |
| `gene_expression.csv` | One-way ANOVA + Tukey HSD   | 5 of 6 pairs           |
| `edited_config.yaml`  | (same as tumor_volume)      | yes                    |

## Generated figures

The complete reference bundles are committed under `expected/`. Each folder
holds the figure, the standalone script, the stats, the input copy, the data
log, the config, the manifest, and the methods paragraph. Regenerate them with
`./example/generate_expected.sh` (deterministic names, so they diff cleanly).

Parametric (t-test):

![t-test](expected/two_group_compare_volume_ttest/figure_volume_ttest.png)

Non-parametric (Mann-Whitney):

![Mann-Whitney](expected/two_group_compare_concentration_mwu/figure_concentration_mwu.png)

Render from edited config (violin, custom palette, numeric p):

![violin](expected/two_group_compare_volume_violin/figure_volume_violin.png)

Multi-group (one-way ANOVA, Tukey HSD, significant pairs only):

![ANOVA](expected/multi_group_compare_expression_anova/figure_expression_anova.png)
