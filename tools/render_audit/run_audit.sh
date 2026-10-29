#!/usr/bin/env bash
# Render every recipe on the stress datasets with both engines, then put the R
# and the Python figure of each case side by side for a visual review.
#
# Usage: tools/render_audit/run_audit.sh [CASE_FILTER]
#   CASE_FILTER  run only the cases whose name contains this text.
#
# The datasets come from make_stress_data.R. The bundles go to
# pubready_output/render_audit/bundles and the side-by-side images to
# pubready_output/render_audit/pairs. A run that fails is a defect too, so the
# script records it in failures.txt and goes on with the next case.
set -u
cd "$(dirname "$0")/../.."

ROOT="pubready_output/render_audit"
DATA="$ROOT/data"
BUNDLES="$ROOT/bundles"
PAIRS="$ROOT/pairs"
FILTER="${1:-}"
VERSION="$(cat VERSION)"
IMAGE="${PUBREADY_IMAGE:-localhost/pubready:$VERSION}"

mkdir -p "$DATA" "$BUNDLES" "$PAIRS"
if [ ! -f "$DATA/two_group_long_labels.csv" ]; then
  podman --log-level=error run --rm -v "$PWD":/work -w /work "$IMAGE" \
    Rscript tools/render_audit/make_stress_data.R "$DATA"
fi

# name | recipe | data file | extra figkit arguments
CASES=(
  "two_group_long_labels|two_group_compare|two_group_long_labels.csv|--x group --y value"
  "two_group_small_n|two_group_compare|two_group_small_n.csv|--x group --y value"
  "two_group_large_n_violin|two_group_compare|two_group_large_n.csv|--x group --y value --geom violin"
  "two_group_large_n_raincloud|two_group_compare|two_group_large_n.csv|--x group --y value --geom raincloud"
  "two_group_large_values|two_group_compare|two_group_large_values.csv|--x group --y value --geom bar"
  "multi_group_eight|multi_group_compare|multi_group_eight.csv|--x genotype --y expression"
  "multi_group_eight_bar|multi_group_compare|multi_group_eight.csv|--x genotype --y expression --geom bar"
  "factorial_four_by_three|factorial_anova|factorial_four_by_three.csv|--x strain --fill treatment --y response"
  "survival_four_arms|survival_km|survival_four_arms.csv|--time months --event event --x arm"
  "cox_many_levels|cox_forest|survival_four_arms.csv|--time months --event event --covariates arm,age,sex,ecog,histology,ldh_u_l"
  "heatmap_large|heatmap|heatmap_large_matrix.csv|--annotation DATA/heatmap_large_annotation.csv"
  "correlation_large_n|correlation|correlation_large_n.csv|--x dose_mg_per_kg --y plasma_concentration_ng_ml"
  "correlation_heatmap_twenty|correlation_heatmap|correlation_heatmap_twenty.csv|"
  "upset_seven_sets|upset|upset_seven_sets.csv|"
  "volcano_extreme|volcano|volcano_extreme.csv|--y padj"
  "enrichment_long_terms|enrichment_dot|enrichment_long_terms.csv|--top_n 25"
  "paired_many_subjects|paired_compare|paired_many_subjects.csv|--x condition --y value --id subject"
  "proportions_five_by_four|proportions|proportions_five_by_four.csv|--x arm --y response"
  "pca_six_groups|pca|pca_six_groups.csv|--group group"
  "spider_sixty|spider_response|spider_sixty_patients.csv|--x week --y lesion_mm --id patient --group arm"
  "waterfall_sixty|waterfall|spider_sixty_patients.csv|--x week --y lesion_mm --id patient --group arm"
  "dose_response_six|dose_response|dose_response_six.csv|--x concentration_um --y viability_pct --group compound"
  "roc_four_markers|roc_curve|roc_four_markers.csv|--y outcome --x troponin,nt_probnp,crp,inverse_marker --positive case"
  "bland_altman_large|bland_altman|bland_altman_large.csv|--x candidate --y reference"
  "forest_meta_twenty_five|forest_meta|forest_meta_twenty_five.csv|--y hr --ci_low lower --ci_high upper --label study"
  "swimmer_seventy|swimmer|swimmer_seventy.csv|--id patient --time months --event ended --group arm --markers response_month,progression_month"
  "timecourse_four_arms|timecourse|timecourse_four_arms.csv|--x day --y volume_mm3 --id mouse --group arm"
  "factorial_interaction_four_by_three|factorial_anova|factorial_four_by_three.csv|--x strain --fill treatment --y response --geom interaction"
  "factorial_bar_four_by_three|factorial_anova|factorial_four_by_three.csv|--x strain --fill treatment --y response --geom bar"
  "estimation_eight_groups|estimation_plot|multi_group_eight.csv|--x genotype --y expression --ref WT"
  "estimation_long_labels|estimation_plot|two_group_long_labels.csv|--x group --y value"
  "enzyme_five_conditions|enzyme_kinetics|enzyme_five_conditions.csv|--x substrate_um --y rate --group condition"
  "standard_curve_plate|standard_curve|standard_curve_plate.csv|--x conc_ng_ml --y od450 --id well"
  "cumulative_incidence_three_causes|cumulative_incidence|cumulative_incidence_three_causes.csv|--time months --event status --x arm"
  "calibration_three_models|calibration|calibration_three_models.csv|--y outcome --x model_calibrated,model_overfit,model_shifted --positive event"
  "sankey_five_visits|sankey|sankey_five_visits.csv|--id patient --x cycle --y best_response"
  "ma_extreme|ma_plot|ma_extreme.csv|--y padj"
  "manhattan_200k|manhattan|manhattan_200k.csv|--group chrom --x bp --y pval --label snp"
  "feature_panel_wide|feature_panel|feature_panel_wide.csv|--x group"
  "embedding_25_clusters|embedding|embedding_25_clusters.csv|--x UMAP_1 --y UMAP_2 --group seurat_clusters"
  "distribution_eight_groups|distribution|distribution_eight_groups.csv|--x cohort --y value"
  "composition_forty_samples|composition|composition_forty_samples.csv|--id sample --x cell_type --y n_cells --group condition"
)

: > "$ROOT/failures.txt"
for entry in "${CASES[@]}"; do
  IFS='|' read -r name recipe file extra <<< "$entry"
  if [ -n "$FILTER" ] && [[ "$name" != *"$FILTER"* ]]; then continue; fi
  extra="${extra//DATA/$DATA}"
  for engine in r python; do
    out="$BUNDLES/${name}_${engine}"
    rm -rf "$out"
    # shellcheck disable=SC2086
    if ! ./cli/figkit plot --engine "$engine" --recipe "$recipe" --data "$DATA/$file" \
        $extra --stamp "$name" --out "$out" > "$out.log" 2>&1; then
      echo "$name $engine: $(tail -1 "$out.log")" >> "$ROOT/failures.txt"
      echo "FAIL: $name ($engine)"
    else
      echo "ok:   $name ($engine)"
    fi
  done
  if command -v magick > /dev/null; then
    r_png="$(ls "$BUNDLES/${name}_r"/*/figure_*.png 2> /dev/null | head -1)"
    p_png="$(ls "$BUNDLES/${name}_python"/*/figure_*.png 2> /dev/null | head -1)"
    if [ -n "$r_png" ] && [ -n "$p_png" ]; then
      magick "$r_png" "$p_png" -resize x900 +append -resize '2000x>' "$PAIRS/$name.png"
    fi
  fi
done

echo
if [ -s "$ROOT/failures.txt" ]; then
  echo "Failed runs:"
  cat "$ROOT/failures.txt"
else
  echo "Every run wrote a bundle."
fi
