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

| Input                 | Test chosen                 | Significant |
|-----------------------|-----------------------------|-------------|
| `tumor_volume.csv`    | Student's two-sample t-test | yes         |
| `cytokine_pg_ml.csv`  | Mann-Whitney U test         | yes         |
| `edited_config.yaml`  | (same as tumor_volume)      | yes         |
