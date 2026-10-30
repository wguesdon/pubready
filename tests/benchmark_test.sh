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

bench_on() {  # $1 = engines, $2 = case name, $3 = text the methods paragraph must hold; rest = plot args
  local engines="$1" name="$2" want="$3"; shift 3
  if [ -n "$FILTER" ] && [[ "$name" != *"$FILTER"* ]]; then return; fi
  local engine dir
  for engine in $engines; do
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

bench() { bench_on "r python" "$@"; }  # $1 = case name, $2 = methods text; rest = plot args

# A case whose method only the R engine has. The Python engine stops on it and
# routes to --engine r, which the unit suites assert.
bench_r() { bench_on "r" "$@"; }

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
bench mtcars_heatmap "Welch" --recipe heatmap --data $D/mtcars_matrix.csv \
  --annotation $D/mtcars_annotation.csv --test welch_t
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
  --annotation $D/airway_annotation.csv --test welch_t
bench airway_pca "exact" --recipe pca --data $D/airway_pca.csv --group dex
bench airway_enrichment "over-representation" --recipe enrichment_dot \
  --data $D/airway_enrichment.csv
bench airway_upset "UpSet" --recipe upset --data $D/airway_upset.csv
bench anorexia_paired "paired t-test" --recipe paired_compare --data $D/anorexia_ft.csv \
  --x condition --y weight_lbs --id subject --test t --levels pre,post
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

echo "--- Test paths that the first cases did not cover ---"
bench mtcars_fisher "Fisher" --recipe proportions --data $D/mtcars_fisher.csv \
  --x transmission --y engine --test fisher
bench infert_fisher "Fisher" --recipe proportions --data $D/infert.csv \
  --x spontaneous --y status --test fisher --levels "0,1,2+"
bench swiss_kendall "Kendall" --recipe correlation --data $D/swiss.csv \
  --x Education --y Fertility --test kendall
bench chickwts_welch "Welch" --recipe multi_group_compare --data $D/chickwts.csv \
  --x feed --y weight --test welch
bench sleep_wilcoxon "Wilcoxon" --recipe paired_compare --data $D/sleep.csv \
  --x group --y extra --id ID --test wilcoxon
bench normand_meta "DerSimonian-Laird" --recipe forest_meta --data $D/normand_stroke.csv \
  --y mean_difference --se se --label study --measure difference

echo "--- Cases of the v3 numbers audit ---"
bench cats_student "Student" --recipe two_group_compare --data $D/cats.csv \
  --x Sex --y Hwt --test student_t
bench cats_mann_whitney "Mann-Whitney" --recipe two_group_compare --data $D/cats.csv \
  --x Sex --y Hwt --test wilcoxon
bench plant_growth_mw_exact "Mann-Whitney" --recipe two_group_compare \
  --data $D/plant_growth_ctrl_trt2.csv --x group --y weight --test wilcoxon
bench sleep_two_group_paired "paired t-test" --recipe two_group_compare --data $D/sleep.csv \
  --x group --y extra --paired --test t
bench sleep_two_group_paired_auto "signed-rank" --recipe two_group_compare --data $D/sleep.csv \
  --x group --y extra --paired
bench nitrofen_two_group "Welch" --recipe two_group_compare --data $D/nitrofen_80_160.csv \
  --x conc --y total --test welch_t
bench cars93_origin "Mann-Whitney" --recipe two_group_compare --data $D/cars93.csv \
  --x Origin --y MPG.highway --test wilcoxon
bench cats_feature_panel "Welch" --recipe feature_panel --data $D/cats.csv \
  --x Sex --test welch_t
bench plant_growth_tukey "Tukey" --recipe multi_group_compare --data $D/plant_growth.csv \
  --x group --y weight --test anova
bench chickwts_games_howell "Games-Howell" --recipe multi_group_compare --data $D/chickwts.csv \
  --x feed --y weight --test welch
bench insect_sprays_dunn "Dunn" --recipe multi_group_compare --data $D/insect_sprays.csv \
  --x spray --y count --test kruskal
bench nitrofen_multi "Tukey" --recipe multi_group_compare --data $D/nitrofen.csv \
  --x conc --y total --test anova
bench quine_threeway_anova "3-way ANOVA" --recipe factorial_anova --data $D/quine.csv \
  --x Eth --fill Sex --facet Lrn --y Days --test anova
bench_r quine_threeway_art "aligned rank transform" --recipe factorial_anova --data $D/quine.csv \
  --x Eth --fill Sex --facet Lrn --y Days --test art
bench_r quine_threeway_auto "aligned rank transform" --recipe factorial_anova --data $D/quine.csv \
  --x Eth --fill Sex --facet Lrn --y Days
bench_r higgins_art "aligned rank transform" --recipe factorial_anova \
  --data $D/higgins1990_table1.csv --x Row --fill Column --y Response --test art
bench shoes_paired_estimation "Hedges" --recipe estimation_plot --data $D/shoes.csv \
  --x material --y wear --paired --id boy --ref A
bench oxboys_rm "repeated-measures" --recipe timecourse --data $D/oxboys.csv \
  --x Occasion --y height --id Subject
bench veteran_trt_km "log-rank" --recipe survival_km --data $D/veteran.csv \
  --time time --event status --x trt --xlab "Days"
bench colon_cox "Cox" --recipe cox_forest --data $D/colon_death.csv \
  --time time --event status --covariates "rx,sex,age"
bench cgd_cox "Cox" --recipe cox_forest --data $D/cgd_first.csv \
  --time time --event status --covariates "treat,inherit,age"
bench birthwt_roc "DeLong" --recipe roc_curve --data $D/birthwt.csv \
  --y low --x lwt --positive 1
bench davis_bland_altman "Bland and Altman" --recipe bland_altman --data $D/davis.csv \
  --x weight --y repwt
bench nitrendipene_ic50 "extra-sum-of-squares" --recipe dose_response --data $D/nitrendipene.csv \
  --x NIF --y activity --group Tissue
bench dnase_ec50_ci "log10 scale" --recipe dose_response --data $D/dnase_run1.csv \
  --x conc --y density
bench formaldehyde_linear "A straight line" --recipe standard_curve \
  --data $D/formaldehyde_standard_curve.csv --x carb --y optden --id sample --test linear
bench housing_kruskal "Kruskal-Wallis" --recipe composition --data $D/housing.csv \
  --id sample --x Sat --y Freq --group Infl
bench airway_pca_permanova "exact" --recipe pca --data $D/airway_pca.csv --group dex
bench iris_pca_permanova "PERMANOVA" --recipe pca --data $D/iris.csv --group Species
bench mtcars_pearson_ci "Pearson" --recipe correlation --data $D/mtcars.csv \
  --x wt --y mpg --test pearson
bench puromycin_ci "Michaelis-Menten" --recipe enzyme_kinetics --data $D/puromycin.csv \
  --x conc --y rate --group state
bench sleep_paired_effect "paired t-test" --recipe paired_compare --data $D/sleep.csv \
  --x group --y extra --id ID --test t

echo
if [ "$ran" -eq 0 ]; then
  echo "NO BENCHMARK MATCHED '$FILTER'"; exit 1
elif [ "$fail" -eq 0 ]; then
  echo "ALL $ran BENCHMARK RUNS PASSED"
else
  echo "BENCHMARK RUNS FAILED"; exit 1
fi
