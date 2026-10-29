#!/usr/bin/env bash
# Benchmark test. Runs every recipe on public benchmark datasets (benchmarks/data)
# through figkit, with both engines, and asserts that each run writes a complete
# bundle. The numbers themselves are checked against benchmarks/reference_values.csv
# by the unit suites (figkit test), which compare each recipe with an independent
# call to the canonical R function. Requires the image: figkit build.
#
# Usage: tests/benchmark_test.sh [case-name-substring]
cd "$(dirname "$0")/.."

OUT="test_output/benchmarks"
FILTER="${1:-}"
rm -rf "$OUT"
fail=0
ran=0

expect_files=(
  "figure_*.pdf" "figure_*.png" "figure_*.svg"
  "stats_*.csv" "manifest_*.json" "methods_*.md"
  "data_log.md" "plot_config.yaml" "session_info.txt"
  "REPRODUCE.md" "script_*" "input_*"
  "report.qmd" "report.html" "report_provenance.json"
)

check_bundle() {  # $1 = output root, $2 = label
  local d
  d="$(find "$1" -mindepth 1 -maxdepth 1 -type d | head -1)"
  if [ -z "$d" ]; then echo "FAIL: $2 wrote no bundle"; fail=1; return 1; fi
  local missing=0 pat
  for pat in "${expect_files[@]}"; do
    if ! compgen -G "$d/$pat" > /dev/null; then
      echo "  MISSING in $2: $pat"; missing=1; fail=1
    fi
  done
  [ "$missing" -eq 0 ] && echo "PASS: $2"
}

bench() {  # $1 = case name, $2 = text the methods paragraph must hold; rest = plot args
  local name="$1" want="$2"; shift 2
  if [ -n "$FILTER" ] && [[ "$name" != *"$FILTER"* ]]; then return; fi
  local engine dir
  for engine in r python; do
    dir="$OUT/${name}_${engine}"
    ran=$((ran + 1))
    if ! ./cli/figkit plot --engine "$engine" "$@" --out "$dir" > "$dir.log" 2>&1; then
      echo "FAIL: $name ($engine) exited with an error; see $dir.log"; fail=1; continue
    fi
    check_bundle "$dir" "$name ($engine)" || continue
    if [ -n "$want" ] && ! grep -q -- "$want" "$dir"/*/methods_*.md; then
      echo "  FAIL: $name ($engine) methods do not mention '$want'"; fail=1
    fi
  done
}

mkdir -p "$OUT"
D=benchmarks/data

echo "--- Existing recipes on public datasets ---"
bench sleep_paired "paired t-test" --recipe paired_compare --data $D/sleep.csv --test t \
  --x group --y extra --id ID
bench sleep_welch "Welch" --recipe two_group_compare --data $D/sleep.csv \
  --x group --y extra --test welch_t
bench plant_growth_anova "one-way ANOVA" --recipe multi_group_compare --data $D/plant_growth.csv \
  --x group --y weight
bench insect_sprays_kruskal "Kruskal-Wallis" --recipe multi_group_compare --data $D/insect_sprays.csv \
  --x spray --y count --test kruskal
bench tooth_growth_twoway "2-way ANOVA" --recipe factorial_anova --data $D/tooth_growth.csv \
  --y len --x dose --fill supp --test anova
bench lung_km "log-rank" --recipe survival_km --data $D/lung.csv \
  --time time --event event --x sex --xlab "Days"
bench lung_cox "Cox" --recipe cox_forest --data $D/lung.csv \
  --time time --event event --covariates "age,sex,ph.ecog"
bench mtcars_correlation "Pearson" --recipe correlation --data $D/mtcars.csv \
  --x wt --y mpg --test pearson
bench mtcars_correlation_heatmap "" --recipe correlation_heatmap --data $D/mtcars.csv
bench mtcars_heatmap "" --recipe heatmap --data $D/mtcars_matrix.csv \
  --annotation $D/mtcars_annotation.csv
bench airquality_spearman "Spearman" --recipe correlation --data $D/airquality.csv \
  --x Temp --y Ozone --test spearman
bench iris_pca "PERMANOVA" --recipe pca --data $D/iris.csv --group Species
bench titanic_proportions "chi-squared" --recipe proportions --data $D/titanic.csv \
  --x Class --y Survived --test chisq
bench ucb_proportions "chi-squared" --recipe proportions --data $D/ucb_admissions.csv \
  --x Gender --y Admit --test chisq

echo "--- New recipes on public datasets ---"
bench dnase_4pl "four-parameter logistic" --recipe dose_response --data $D/dnase_run1.csv \
  --x conc --y density
bench dnase_shared "extra-sum-of-squares" --recipe dose_response --data $D/dnase.csv \
  --x conc --y density --group Run
bench asah_roc "DeLong" --recipe roc_curve --data $D/asah.csv \
  --y outcome --x s100b,ndka
bench pefr_bland_altman "Bland and Altman" --recipe bland_altman --data $D/pefr.csv \
  --x wright_1 --y mini_1
bench bcg_meta "DerSimonian-Laird" --recipe forest_meta --data $D/bcg.csv \
  --y risk_ratio --ci_low ci_low --ci_high ci_high --label trial
bench orthodont_timecourse "mixed ANOVA" --recipe timecourse --data $D/orthodont.csv \
  --x age --y distance --id Subject --group Sex
bench chick_weight_timecourse "Greenhouse-Geisser" --recipe timecourse --data $D/chick_weight.csv \
  --x Time --y weight --id Chick --group Diet
bench aml_swimmer "median duration" --recipe swimmer --data $D/aml.csv \
  --id patient --time time --event status --group x --xlab "Weeks"

echo "--- Real omics and clinical datasets ---"
bench airway_volcano "volcano" --recipe volcano --data $D/airway_de.csv --y padj
bench airway_heatmap "Welch" --recipe heatmap --data $D/airway_heatmap.csv \
  --annotation $D/airway_annotation.csv
bench airway_pca "exact" --recipe pca --data $D/airway_pca.csv --group dex
bench airway_enrichment "over-representation" --recipe enrichment_dot \
  --data $D/airway_enrichment.csv
bench airway_upset "UpSet" --recipe upset --data $D/airway_upset.csv
bench anorexia_paired "paired t-test" --recipe paired_compare --data $D/anorexia_ft.csv \
  --x condition --y weight_lbs --id subject --test t
bench veteran_km "log-rank" --recipe survival_km --data $D/veteran.csv \
  --time time --event status --x celltype --xlab "Days"
bench veteran_cox "Cox" --recipe cox_forest --data $D/veteran.csv \
  --time time --event status --covariates "trt,celltype,karno,age"
bench state_correlation_heatmap "" --recipe correlation_heatmap --data $D/state_x77.csv

echo "--- Distribution and composition ---"
bench faithful_distribution "Freedman-Diaconis" --recipe distribution --data $D/faithful.csv \
  --y eruptions
bench tooth_growth_distribution "Shapiro-Wilk" --recipe distribution --data $D/tooth_growth.csv \
  --x supp --y len
bench housing_composition "Mann-Whitney" --recipe composition --data $D/housing.csv \
  --id sample --x Sat --y Freq --group Cont
echo "--- Bench assay recipes on public datasets ---"
bench chickwts_estimation "Park-Miller" --recipe estimation_plot --data $D/chickwts.csv \
  --x feed --y weight --ref casein
bench puromycin_kinetics "Michaelis-Menten" --recipe enzyme_kinetics --data $D/puromycin.csv \
  --x conc --y rate --group state
bench dnase_standard_curve "four-parameter logistic" --recipe standard_curve \
  --data $D/dnase_standard_curve.csv --x conc --y density --id sample

echo "--- Clinical recipes on public datasets ---"
bench mgus2_cif "Gray" --recipe cumulative_incidence --data $D/mgus2.csv \
  --time etime --event event --x sex --xlab "Months"
bench pima_calibration "Hosmer-Lemeshow" --recipe calibration --data $D/pima_te.csv \
  --y type --x predicted --positive Yes
bench pbcseq_sankey "Bowker" --recipe sankey --data $D/pbcseq_edema.csv \
  --id id --x year --y edema

echo "--- Omics recipes on public datasets ---"
bench airway_ma "MA plot" --recipe ma_plot --data $D/airway_de.csv
bench gwas_manhattan "lambda_GC" --recipe manhattan --data $D/gwas_results.csv \
  --group CHR --x BP --y P --label SNP
bench iris_feature_panel_kruskal "Kruskal-Wallis" --recipe feature_panel --data $D/iris.csv \
  --x Species
bench iris_feature_panel_anova "one-way ANOVA" --recipe feature_panel --data $D/iris.csv \
  --x Species --test anova
bench iris_embedding "silhouette" --recipe embedding --data $D/iris.csv \
  --x Petal.Length --y Petal.Width --group Species

echo
if [ "$ran" -eq 0 ]; then
  echo "NO BENCHMARK MATCHED '$FILTER'"; exit 1
elif [ "$fail" -eq 0 ]; then
  echo "ALL $ran BENCHMARK RUNS PASSED"
else
  echo "BENCHMARK RUNS FAILED"; exit 1
fi
