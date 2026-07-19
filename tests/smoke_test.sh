#!/usr/bin/env bash
# End-to-end smoke test. Runs the example cases through figkit and asserts each
# writes a complete artifact bundle and picks the expected test. Requires the
# image: build it first with `figkit build`.
cd "$(dirname "$0")/.."

OUT="test_output"
rm -rf "$OUT"
fail=0

expect_files=(
  "figure_*.pdf" "figure_*.png" "figure_*.svg"
  "stats_*.csv" "manifest_*.json" "methods_*.md"
  "data_log.md" "plot_config.yaml" "session_info.txt"
  "REPRODUCE.md" "script_*" "input_*"
)

check_bundle() {  # $1 = output root, $2 = label
  local d
  d="$(find "$1" -mindepth 1 -maxdepth 1 -type d | head -1)"
  if [ -z "$d" ]; then echo "FAIL: no bundle written for $2"; fail=1; return; fi
  local missing=0
  for pat in "${expect_files[@]}"; do
    if ! compgen -G "$d/$pat" > /dev/null; then
      echo "  MISSING: $pat"; missing=1; fail=1
    fi
  done
  [ "$missing" -eq 0 ] && echo "PASS: $2 wrote a complete bundle -> $d"
}

echo "[1] parametric path (expect t-test)"
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --ylab "Tumor volume (mm^3)" --out "$OUT/t_test"
check_bundle "$OUT/t_test" "tumor_volume"
if compgen -G "$OUT/t_test/*/stats_*.csv" > /dev/null && grep -q "t-test" "$OUT"/t_test/*/stats_*.csv; then
  echo "PASS: chose a t-test"
else
  echo "FAIL: expected a t-test"; fail=1
fi

echo "[2] non-parametric path (expect Mann-Whitney)"
./cli/figkit plot --recipe two_group_compare --data example/cytokine_pg_ml.csv \
  --x group --y concentration --ylab "IL-6 (pg/mL)" --out "$OUT/mwu"
check_bundle "$OUT/mwu" "cytokine_pg_ml"
if compgen -G "$OUT/mwu/*/stats_*.csv" > /dev/null && grep -q "Mann-Whitney" "$OUT"/mwu/*/stats_*.csv; then
  echo "PASS: chose Mann-Whitney"
else
  echo "FAIL: expected Mann-Whitney"; fail=1
fi

echo "[3] render from edited config"
./cli/figkit render --config example/edited_config.yaml --out "$OUT/render"
check_bundle "$OUT/render" "edited_config"

echo "[4] multi-group path (expect one-way ANOVA + Tukey)"
./cli/figkit plot --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --ylab "Expression (a.u.)" --out "$OUT/anova"
check_bundle "$OUT/anova" "gene_expression"
if compgen -G "$OUT/anova/*/methods_*.md" > /dev/null && grep -q "one-way ANOVA" "$OUT"/anova/*/methods_*.md; then
  echo "PASS: chose one-way ANOVA"
else
  echo "FAIL: expected one-way ANOVA"; fail=1
fi

echo "[5] two-way factorial ANOVA"
./cli/figkit plot --recipe factorial_anova --data example/twoway_response.csv \
  --y response --x genotype --fill treatment --ylab "Response (a.u.)" --out "$OUT/twoway"
check_bundle "$OUT/twoway" "twoway_response"
if compgen -G "$OUT/twoway/*/stats_*.csv" > /dev/null && grep -q "genotype:treatment" "$OUT"/twoway/*/stats_*.csv; then
  echo "PASS: two-way interaction term present"
else
  echo "FAIL: expected interaction term"; fail=1
fi

echo "[6] three-way factorial ANOVA"
./cli/figkit plot --recipe factorial_anova --data example/threeway_response.csv \
  --y response --x genotype --fill treatment --facet sex --ylab "Response (a.u.)" --out "$OUT/threeway"
check_bundle "$OUT/threeway" "threeway_response"
if compgen -G "$OUT/threeway/*/stats_*.csv" > /dev/null && grep -q "genotype:treatment:sex" "$OUT"/threeway/*/stats_*.csv; then
  echo "PASS: three-way term present"
else
  echo "FAIL: expected three-way term"; fail=1
fi

echo "[7] Kaplan-Meier + log-rank"
./cli/figkit plot --recipe survival_km --data example/survival_trial.csv \
  --time time --event event --x arm --xlab "Months" --out "$OUT/km"
check_bundle "$OUT/km" "survival_km"
if compgen -G "$OUT/km/*/stats_*.csv" > /dev/null && grep -q "log_rank_p" "$OUT"/km/*/stats_*.csv; then
  echo "PASS: KM with log-rank"
else
  echo "FAIL: expected KM stats"; fail=1
fi

echo "[8] Cox hazard-ratio forest"
./cli/figkit plot --recipe cox_forest --data example/survival_trial.csv \
  --time time --event event --covariates "arm,age,sex,stage" --out "$OUT/cox"
check_bundle "$OUT/cox" "cox_forest"
if compgen -G "$OUT/cox/*/stats_*.csv" > /dev/null && grep -q "armtreated" "$OUT"/cox/*/stats_*.csv; then
  echo "PASS: Cox HR table"
else
  echo "FAIL: expected Cox HR table"; fail=1
fi

echo "[9] clustered heatmap"
./cli/figkit plot --recipe heatmap --data example/expression_matrix.csv \
  --annotation example/expression_annotation.csv --out "$OUT/heatmap"
check_bundle "$OUT/heatmap" "heatmap"
if compgen -G "$OUT/heatmap/*/stats_*.csv" > /dev/null && grep -q "p_adj" "$OUT"/heatmap/*/stats_*.csv; then
  echo "PASS: heatmap per-feature stats"
else
  echo "FAIL: expected heatmap stats"; fail=1
fi
# the annotation must be copied into the bundle for standalone reproduction
if compgen -G "$OUT/heatmap/*/expression_annotation.csv" > /dev/null; then
  echo "PASS: annotation copied into bundle"
else
  echo "FAIL: annotation not copied"; fail=1
fi

echo "[10] correlation scatter (expect Pearson)"
./cli/figkit plot --recipe correlation --data example/correlation_xy.csv \
  --x gene_a --y gene_b --out "$OUT/corr"
check_bundle "$OUT/corr" "correlation"
if compgen -G "$OUT/corr/*/stats_*.csv" > /dev/null && grep -q "Pearson correlation" "$OUT"/corr/*/stats_*.csv; then
  echo "PASS: chose Pearson correlation"
else
  echo "FAIL: expected Pearson correlation"; fail=1
fi

echo "[11] clustered correlation heatmap"
./cli/figkit plot --recipe correlation_heatmap --data example/correlation_vars.csv --out "$OUT/corrhm"
check_bundle "$OUT/corrhm" "correlation_heatmap"
if compgen -G "$OUT/corrhm/*/stats_*.csv" > /dev/null && grep -q "p_adj" "$OUT"/corrhm/*/stats_*.csv; then
  echo "PASS: correlation heatmap pairwise stats"
else
  echo "FAIL: expected correlation heatmap stats"; fail=1
fi

echo "[12] UpSet plot"
./cli/figkit plot --recipe upset --data example/set_membership.csv --out "$OUT/upset"
check_bundle "$OUT/upset" "upset"
if compgen -G "$OUT/upset/*/stats_*.csv" > /dev/null && grep -q "intersection" "$OUT"/upset/*/stats_*.csv; then
  echo "PASS: UpSet intersection stats"
else
  echo "FAIL: expected UpSet stats"; fail=1
fi

echo "[13] prism theme carries through"
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --theme prism --out "$OUT/prism"
check_bundle "$OUT/prism" "prism_theme"
if compgen -G "$OUT/prism/*/plot_config.yaml" > /dev/null && grep -q "theme: prism" "$OUT"/prism/*/plot_config.yaml; then
  echo "PASS: prism theme recorded in config"
else
  echo "FAIL: expected prism theme in config"; fail=1
fi

echo "--- Python engine (all recipes) ---"
pyc() {  # $1 = short label; rest = plot args
  local label="$1"; shift
  ./cli/figkit plot --engine python "$@" --out "$OUT/py_$label" > /dev/null 2>&1
  if compgen -G "$OUT/py_$label/*/manifest_*.json" > /dev/null && \
     grep -q '"engine": "python"' "$OUT"/py_$label/*/manifest_*.json; then
    echo "PASS: python $label"
  else
    echo "FAIL: python $label"; fail=1
  fi
}
pyc tg  --recipe two_group_compare  --data example/tumor_volume.csv     --x group --y volume
pyc mg  --recipe multi_group_compare --data example/gene_expression.csv --x genotype --y expression
pyc fa2 --recipe factorial_anova    --data example/twoway_response.csv   --y response --x genotype --fill treatment
pyc fa3 --recipe factorial_anova    --data example/threeway_response.csv --y response --x genotype --fill treatment --facet sex
pyc km  --recipe survival_km        --data example/survival_trial.csv    --time time --event event --x arm
pyc cox --recipe cox_forest         --data example/survival_trial.csv    --time time --event event --covariates arm,age,sex,stage
pyc hm  --recipe heatmap            --data example/expression_matrix.csv --annotation example/expression_annotation.csv
pyc corr    --recipe correlation         --data example/correlation_xy.csv   --x gene_a --y gene_b
pyc corrhm  --recipe correlation_heatmap --data example/correlation_vars.csv
pyc upset   --recipe upset               --data example/set_membership.csv
pyc prism   --recipe two_group_compare   --data example/tumor_volume.csv     --x group --y volume --theme prism

echo
if [ "$fail" -eq 0 ]; then
  echo "ALL SMOKE TESTS PASSED"
else
  echo "SMOKE TESTS FAILED"; exit 1
fi
