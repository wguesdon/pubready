#!/usr/bin/env bash
# Regenerate the committed reference bundles under example/expected/.
# Uses fixed --stamp names and a fixed created date so the outputs are
# deterministic and diff cleanly. Run from anywhere; requires the built image.
# Set PUBPLOT_EXPECTED_OUT to write elsewhere (the determinism check points it
# at a temp tree so it can diff against the committed bundles).
set -e
cd "$(dirname "$0")/../.."

export PUBPLOT_CREATED="2026-07-18T00:00:00Z"
OUT="${PUBPLOT_EXPECTED_OUT:-example/expected}"
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

# Correlation scatter (auto -> Pearson) with a linear fit and CI band
./cli/figkit plot --recipe correlation --data example/correlation_xy.csv \
  --x gene_a --y gene_b --xlab "Gene A (a.u.)" --ylab "Gene B (a.u.)" \
  --stamp corr --out "$OUT"

# Clustered correlation heatmap (Pearson) with coefficients and significance stars
./cli/figkit plot --recipe correlation_heatmap --data example/correlation_vars.csv \
  --stamp corrhm --out "$OUT"

# UpSet plot of set intersections from a binary membership matrix
./cli/figkit plot --recipe upset --data example/set_membership.csv \
  --stamp upset --out "$OUT"

# Prism theme (ggprism) on the two-group comparison
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --xlab "Group" --theme prism \
  --stamp prism --out "$OUT"

# Volcano plot (EnhancedVolcano style) from the synthetic DE table, FDR-thresholded.
# The ylab uses the --flag=value form because optparse reads a leading-dash value
# (-log10 ...) as another flag.
./cli/figkit plot --recipe volcano --data example/de_results.csv \
  --y padj "--ylab=-log10 adjusted p" --stamp volcano --out "$OUT"

# Enrichment dot plot (ORA) from a clusterProfiler-like table
./cli/figkit plot --recipe enrichment_dot --data example/enrichment_results.csv --stamp enrich --out "$OUT"

# Paired before/after comparison with connecting lines
./cli/figkit plot --recipe paired_compare --data example/paired_response.csv \
  --x condition --y value --id subject --ylab "Marker (a.u.)" --stamp paired --out "$OUT"

# Proportions (chi-square) as a 100% stacked bar
./cli/figkit plot --recipe proportions --data example/response_by_arm.csv \
  --x arm --y response --stamp prop --out "$OUT"

# PCA scatter with 95% ellipses and a PERMANOVA p
./cli/figkit plot --recipe pca --data example/pca_samples.csv --group group --stamp pca --out "$OUT"

# Spider plot: the change from baseline of each patient, with the RECIST lines
./cli/figkit plot --recipe spider_response --data example/spider_response.csv \
  --x week --y target_lesion_mm --id patient --group arm --xlab "Week" \
  --stamp spider --out "$OUT"

# Raincloud and bar geoms on the comparison recipes
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --geom raincloud --stamp raincloud --out "$OUT"
./cli/figkit plot --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)" --geom bar --stamp bar --out "$OUT"

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
./cli/figkit plot $PY --recipe correlation --data example/correlation_xy.csv \
  --x gene_a --y gene_b --xlab "Gene A (a.u.)" --ylab "Gene B (a.u.)" --stamp corr_py --out "$OUT"
./cli/figkit plot $PY --recipe correlation_heatmap --data example/correlation_vars.csv \
  --stamp corrhm_py --out "$OUT"
./cli/figkit plot $PY --recipe upset --data example/set_membership.csv \
  --stamp upset_py --out "$OUT"
./cli/figkit plot $PY --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --xlab "Group" --theme prism --stamp prism_py --out "$OUT"
./cli/figkit plot $PY --recipe volcano --data example/de_results.csv \
  --y padj "--ylab=-log10 adjusted p" --stamp volcano_py --out "$OUT"
./cli/figkit plot $PY --recipe enrichment_dot --data example/enrichment_results.csv --stamp enrich_py --out "$OUT"
./cli/figkit plot $PY --recipe paired_compare --data example/paired_response.csv \
  --x condition --y value --id subject --ylab "Marker (a.u.)" --stamp paired_py --out "$OUT"
./cli/figkit plot $PY --recipe proportions --data example/response_by_arm.csv \
  --x arm --y response --stamp prop_py --out "$OUT"
./cli/figkit plot $PY --recipe pca --data example/pca_samples.csv --group group --stamp pca_py --out "$OUT"
./cli/figkit plot $PY --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --geom raincloud --stamp raincloud_py --out "$OUT"
./cli/figkit plot $PY --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)" --geom bar --stamp bar_py --out "$OUT"
./cli/figkit plot $PY --recipe spider_response --data example/spider_response.csv \
  --x week --y target_lesion_mm --id patient --group arm --xlab "Week" --stamp spider_py --out "$OUT"

echo "Reference bundles written under $OUT/"
