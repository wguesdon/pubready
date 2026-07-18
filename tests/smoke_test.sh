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
  "REPRODUCE.md" "script_*.R" "input_*"
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

echo
if [ "$fail" -eq 0 ]; then
  echo "ALL SMOKE TESTS PASSED"
else
  echo "SMOKE TESTS FAILED"; exit 1
fi
