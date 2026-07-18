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

echo "Reference bundles written under $OUT/"
