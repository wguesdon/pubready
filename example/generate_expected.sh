#!/usr/bin/env bash
# Regenerate the committed reference bundles under example/expected/.
# Uses fixed --stamp names and a fixed created date so the outputs are
# deterministic and diff cleanly. Run from anywhere; requires the built image.
set -e
cd "$(dirname "$0")/.."

export PUBPLOT_CREATED="2026-07-18T00:00:00Z"
OUT="example/expected"
rm -rf "$OUT"

# Parametric path -> Student's t-test
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --xlab "Group" \
  --stamp ttest --out "$OUT"

# Non-parametric path -> Mann-Whitney U
./cli/figkit plot --recipe two_group_compare --data example/cytokine_pg_ml.csv \
  --x group --y concentration --ylab "IL-6 (pg/mL)" --xlab "Group" \
  --stamp mwu --out "$OUT"

# Render path from an edited config -> violin, custom palette, numeric label
./cli/figkit render --config example/edited_config.yaml \
  --stamp violin --out "$OUT"

# Multi-group path -> one-way ANOVA + Tukey, significant-only brackets
./cli/figkit plot --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)" --xlab "Genotype" \
  --stamp anova --out "$OUT"

# Two-way factorial ANOVA -> grouped box with effects on top
./cli/figkit plot --recipe factorial_anova --data example/twoway_response.csv \
  --y response --x genotype --fill treatment --ylab "Response (a.u.)" \
  --stamp twoway --out "$OUT"

# Three-way factorial ANOVA -> faceted grouped box with effects on top
./cli/figkit plot --recipe factorial_anova --data example/threeway_response.csv \
  --y response --x genotype --fill treatment --facet sex --ylab "Response (a.u.)" \
  --stamp threeway --out "$OUT"

# Kaplan-Meier survival curves with log-rank p and risk table
./cli/figkit plot --recipe survival_km --data example/survival_trial.csv \
  --time time --event event --x arm --xlab "Months" \
  --stamp km --out "$OUT"

# Cox proportional-hazards forest plot of hazard ratios
./cli/figkit plot --recipe cox_forest --data example/survival_trial.csv \
  --time time --event event --covariates "arm,age,sex,stage" \
  --stamp cox --out "$OUT"

# Clustered heatmap (row z-score, cluster both) with group annotation + sig stars
./cli/figkit plot --recipe heatmap --data example/expression_matrix.csv \
  --annotation example/expression_annotation.csv \
  --stamp demo --out "$OUT"

# ---------------------------------------------------------------------------
# Python engine references (engine=python, _py stamp) so each recipe has both.
# ---------------------------------------------------------------------------
PY="--engine python"
./cli/figkit plot $PY --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --xlab "Group" --stamp ttest_py --out "$OUT"
./cli/figkit plot $PY --recipe two_group_compare --data example/cytokine_pg_ml.csv \
  --x group --y concentration --ylab "IL-6 (pg/mL)" --xlab "Group" --stamp mwu_py --out "$OUT"
./cli/figkit plot $PY --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)" --xlab "Genotype" --stamp anova_py --out "$OUT"
./cli/figkit plot $PY --recipe factorial_anova --data example/twoway_response.csv \
  --y response --x genotype --fill treatment --ylab "Response (a.u.)" --stamp twoway_py --out "$OUT"
./cli/figkit plot $PY --recipe factorial_anova --data example/threeway_response.csv \
  --y response --x genotype --fill treatment --facet sex --ylab "Response (a.u.)" --stamp threeway_py --out "$OUT"
./cli/figkit plot $PY --recipe survival_km --data example/survival_trial.csv \
  --time time --event event --x arm --xlab "Months" --stamp km_py --out "$OUT"
./cli/figkit plot $PY --recipe cox_forest --data example/survival_trial.csv \
  --time time --event event --covariates "arm,age,sex,stage" --stamp cox_py --out "$OUT"
./cli/figkit plot $PY --recipe heatmap --data example/expression_matrix.csv \
  --annotation example/expression_annotation.csv --stamp demo_py --out "$OUT"

echo "Reference bundles written under $OUT/"
