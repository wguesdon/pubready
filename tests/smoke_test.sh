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
  "report.qmd" "report.html" "report_provenance.json"
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

# The string items of one array in a manifest, one per line. The array starts on
# the line of its key and ends at the first "]".
json_array() {  # $1 = file, $2 = key
  awk -v key="\"$2\":" 'index($0, key) { on = 1 } on { print } on && /\]/ { exit }' "$1" |
    grep -o '"[^"]*"' | grep -v "^\"$2\"\$"
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

echo "[14] volcano plot"
./cli/figkit plot --recipe volcano --data example/de_results.csv --out "$OUT/volcano"
check_bundle "$OUT/volcano" "volcano"
if compgen -G "$OUT/volcano/*/stats_*.csv" > /dev/null && grep -q "direction" "$OUT"/volcano/*/stats_*.csv; then
  echo "PASS: volcano labeled-hits stats"
else
  echo "FAIL: expected volcano stats"; fail=1
fi

echo "[15] enrichment dot plot"
./cli/figkit plot --recipe enrichment_dot --data example/enrichment_results.csv --out "$OUT/enrich"
check_bundle "$OUT/enrich" "enrichment_dot"

echo "[16] paired comparison"
./cli/figkit plot --recipe paired_compare --data example/paired_response.csv \
  --x condition --y value --id subject --out "$OUT/paired"
check_bundle "$OUT/paired" "paired_compare"
if compgen -G "$OUT/paired/*/stats_*.csv" > /dev/null && grep -q "n_pairs" "$OUT"/paired/*/stats_*.csv; then
  echo "PASS: paired stats with n_pairs"
else
  echo "FAIL: expected paired stats"; fail=1
fi

echo "[17] proportions (chi-square / Fisher)"
./cli/figkit plot --recipe proportions --data example/response_by_arm.csv \
  --x arm --y response --out "$OUT/prop"
check_bundle "$OUT/prop" "proportions"
if compgen -G "$OUT/prop/*/methods_*.md" > /dev/null && grep -qi "chi-squared\|Fisher" "$OUT"/prop/*/methods_*.md; then
  echo "PASS: proportions test chosen"
else
  echo "FAIL: expected proportions test"; fail=1
fi

echo "[18] PCA scatter with PERMANOVA"
./cli/figkit plot --recipe pca --data example/pca_samples.csv --group group --out "$OUT/pca"
check_bundle "$OUT/pca" "pca"
if compgen -G "$OUT/pca/*/methods_*.md" > /dev/null && grep -q "PERMANOVA" "$OUT"/pca/*/methods_*.md; then
  echo "PASS: PCA with PERMANOVA"
else
  echo "FAIL: expected PCA PERMANOVA"; fail=1
fi

echo "[19] raincloud geom"
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume --geom raincloud --out "$OUT/raincloud"
check_bundle "$OUT/raincloud" "raincloud"

echo "[20] bar geom (mean + SEM + dots)"
./cli/figkit plot --recipe multi_group_compare --data example/gene_expression.csv \
  --x genotype --y expression --geom bar --out "$OUT/bar"
check_bundle "$OUT/bar" "bar_geom"

echo "[21] spider plot (change from baseline)"
./cli/figkit plot --recipe spider_response --data example/spider_response.csv \
  --x week --y target_lesion_mm --id patient --group arm --out "$OUT/spider"
check_bundle "$OUT/spider" "spider_response"
if compgen -G "$OUT/spider/*/stats_*.csv" > /dev/null && \
   grep -q "category_from_change" "$OUT"/spider/*/stats_*.csv; then
  echo "PASS: spider per-patient categories"
else
  echo "FAIL: expected the spider category column"; fail=1
fi

echo "[22] waterfall of the best change from baseline"
./cli/figkit plot --recipe waterfall --data example/spider_response.csv \
  --x week --y target_lesion_mm --id patient --group arm --out "$OUT/waterfall"
check_bundle "$OUT/waterfall" "waterfall"
if compgen -G "$OUT/waterfall/*/stats_*.csv" > /dev/null && \
   grep -q "best_change_pct" "$OUT"/waterfall/*/stats_*.csv; then
  echo "PASS: waterfall best change per patient"
else
  echo "FAIL: expected the waterfall best change column"; fail=1
fi

echo "[23] dose-response curves (expect IC50 and the shared F test)"
./cli/figkit plot --recipe dose_response --data example/dose_response.csv \
  --x concentration_um --y viability_pct --group compound --out "$OUT/drc"
check_bundle "$OUT/drc" "dose_response"
if compgen -G "$OUT/drc/*/methods_*.md" > /dev/null && \
   grep -q "IC50" "$OUT"/drc/*/methods_*.md && grep -q "extra-sum-of-squares" "$OUT"/drc/*/methods_*.md; then
  echo "PASS: dose-response IC50 and shared F test"
else
  echo "FAIL: expected the IC50 and the shared F test"; fail=1
fi

echo "[24] ROC curves (expect DeLong)"
./cli/figkit plot --recipe roc_curve --data example/biomarker_roc.csv \
  --y status --x marker_a,marker_b --positive disease --out "$OUT/roc"
check_bundle "$OUT/roc" "roc_curve"
if compgen -G "$OUT/roc/*/stats_*.csv" > /dev/null && \
   grep -q "delong_p_vs_first" "$OUT"/roc/*/stats_*.csv && grep -q ",disease," "$OUT"/roc/*/stats_*.csv; then
  echo "PASS: ROC with the named case level and the DeLong test"
else
  echo "FAIL: expected the ROC stats with the case level"; fail=1
fi

echo "[25] Bland-Altman agreement"
./cli/figkit plot --recipe bland_altman --data example/method_agreement.csv \
  --x new_assay --y reference_assay --out "$OUT/ba"
check_bundle "$OUT/ba" "bland_altman"
if compgen -G "$OUT/ba/*/stats_*.csv" > /dev/null && grep -q "loa_high" "$OUT"/ba/*/stats_*.csv && \
   compgen -G "$OUT/ba/*/qc_normality_*.png" > /dev/null; then
  echo "PASS: Bland-Altman limits and the QC panel of the differences"
else
  echo "FAIL: expected the limits of agreement and the QC panel"; fail=1
fi

echo "[26] meta-analysis forest plot (expect DerSimonian-Laird)"
./cli/figkit plot --recipe forest_meta --data example/meta_trials.csv \
  --y hazard_ratio --ci_low ci_low --ci_high ci_high --label trial --out "$OUT/meta"
check_bundle "$OUT/meta" "forest_meta"
if compgen -G "$OUT/meta/*/stats_*.csv" > /dev/null && \
   grep -q "Random effects (DerSimonian-Laird)" "$OUT"/meta/*/stats_*.csv; then
  echo "PASS: pooled random-effects row"
else
  echo "FAIL: expected the random-effects row"; fail=1
fi

echo "[27] repeated-measures time course (expect a mixed ANOVA)"
./cli/figkit plot --recipe timecourse --data example/tumor_growth.csv \
  --x day --y volume_mm3 --id mouse --group arm --out "$OUT/tc"
check_bundle "$OUT/tc" "timecourse"
if compgen -G "$OUT/tc/*/methods_*.md" > /dev/null && grep -q "mixed ANOVA" "$OUT"/tc/*/methods_*.md && \
   grep -q "1 subject(s) without a value" "$OUT"/tc/*/methods_*.md; then
  echo "PASS: mixed ANOVA with the incomplete mouse dropped"
else
  echo "FAIL: expected the mixed ANOVA and the dropped mouse"; fail=1
fi

echo "[28] swimmer plot of the time on treatment"
./cli/figkit plot --recipe swimmer --data example/swimmer_patients.csv \
  --id patient --time months --event ended --group arm \
  --markers response_month,progression_month --out "$OUT/swimmer"
check_bundle "$OUT/swimmer" "swimmer"
if compgen -G "$OUT/swimmer/*/stats_*.csv" > /dev/null && \
   grep -q "ongoing" "$OUT"/swimmer/*/stats_*.csv && grep -q "response_month" "$OUT"/swimmer/*/stats_*.csv; then
  echo "PASS: swimmer per-patient courses with markers"
else
  echo "FAIL: expected the swimmer per-patient table"; fail=1
fi

echo "[29] distribution of one variable per group (expect Shapiro-Wilk, no comparison)"
./cli/figkit plot --recipe distribution --data example/distribution.csv \
  --x condition --y cell_area_um2 --xlab "Cell area (um^2)" --out "$OUT/dist"
check_bundle "$OUT/dist" "distribution"
if compgen -G "$OUT/dist/*/stats_*.csv" > /dev/null && grep -q "shapiro_p" "$OUT"/dist/*/stats_*.csv && \
   grep -q "Freedman-Diaconis" "$OUT"/dist/*/methods_*.md && grep -q "No test compared the groups" "$OUT"/dist/*/methods_*.md; then
  echo "PASS: distribution with Freedman-Diaconis bins and Shapiro-Wilk per group"
else
  echo "FAIL: expected the distribution summary"; fail=1
fi

echo "[30] composition of each sample (expect Mann-Whitney per category, BH-adjusted)"
./cli/figkit plot --recipe composition --data example/composition.csv \
  --id sample --x cell_type --y n_cells --group condition --out "$OUT/comp"
check_bundle "$OUT/comp" "composition"
if compgen -G "$OUT/comp/*/stats_*.csv" > /dev/null && grep -q "Mann-Whitney U test" "$OUT"/comp/*/stats_*.csv && \
   grep -q "compositional" "$OUT"/comp/*/methods_*.md; then
  echo "PASS: composition with a per-category Mann-Whitney test"
else
  echo "FAIL: expected the per-category Mann-Whitney test"; fail=1
fi

echo "[31] factorial ANOVA interaction geom (expect the same effects as the box geom)"
./cli/figkit plot --recipe factorial_anova --data example/twoway_response.csv \
  --y response --x genotype --fill treatment --ylab "Response (a.u.)" --geom interaction \
  --out "$OUT/fa_interaction"
check_bundle "$OUT/fa_interaction" "factorial interaction"
if compgen -G "$OUT/fa_interaction/*/stats_*.csv" > /dev/null && \
   grep -q "genotype:treatment" "$OUT"/fa_interaction/*/stats_*.csv && \
   grep -q "position_dodge(width = 0.25)" "$OUT"/fa_interaction/*/script_*.R; then
  echo "PASS: interaction plot with the interaction term"
else
  echo "FAIL: expected the interaction plot and term"; fail=1
fi

echo "[32] estimation plot, Cumming layout (expect Welch t-tests and BCa bootstrap CIs)"
./cli/figkit plot --recipe estimation_plot --data example/estimation_plot.csv \
  --x treatment --y speed_um_h --ref vehicle --out "$OUT/estimation"
check_bundle "$OUT/estimation" "estimation_plot"
if compgen -G "$OUT/estimation/*/stats_*.csv" > /dev/null && \
   grep -q "Welch two-sample t-test" "$OUT"/estimation/*/stats_*.csv && \
   grep -q "ci_low" "$OUT"/estimation/*/stats_*.csv && \
   grep -q "Park-Miller" "$OUT"/estimation/*/methods_*.md; then
  echo "PASS: mean differences against vehicle with BCa CIs"
else
  echo "FAIL: expected the Welch t-tests and the BCa CIs"; fail=1
fi

echo "[33] estimation plot, Gardner-Altman paired layout (expect a paired t-test)"
./cli/figkit plot --recipe estimation_plot --data example/paired_response.csv \
  --x condition --y value --id subject --paired --ref pre --out "$OUT/estimation_paired"
check_bundle "$OUT/estimation_paired" "estimation_plot paired"
if compgen -G "$OUT/estimation_paired/*/stats_*.csv" > /dev/null && \
   grep -q "paired t-test" "$OUT"/estimation_paired/*/stats_*.csv; then
  echo "PASS: paired mean difference with a BCa CI"
else
  echo "FAIL: expected the paired t-test"; fail=1
fi

echo "[34] enzyme kinetics (expect Michaelis-Menten fits and the shared Km and Vmax tests)"
./cli/figkit plot --recipe enzyme_kinetics --data example/enzyme_kinetics.csv \
  --x substrate_um --y rate --group inhibitor --out "$OUT/enzyme"
check_bundle "$OUT/enzyme" "enzyme_kinetics"
if compgen -G "$OUT/enzyme/*/methods_*.md" > /dev/null && \
   grep -q "Michaelis-Menten" "$OUT"/enzyme/*/methods_*.md && \
   grep -q "shared Km F" "$OUT"/enzyme/*/methods_*.md && grep -q "km_ci_low" "$OUT"/enzyme/*/stats_*.csv; then
  echo "PASS: Michaelis-Menten fits with the shared-parameter F tests"
else
  echo "FAIL: expected the Michaelis-Menten fits and the F tests"; fail=1
fi

echo "[35] standard curve (expect a 4PL and the out-of-range unknowns flagged)"
./cli/figkit plot --recipe standard_curve --data example/standard_curve.csv \
  --x conc_ng_ml --y absorbance --id sample --out "$OUT/stdcurve"
check_bundle "$OUT/stdcurve" "standard_curve"
if compgen -G "$OUT/stdcurve/*/stats_*.csv" > /dev/null && \
   grep -q "four-parameter logistic" "$OUT"/stdcurve/*/methods_*.md && \
   grep -q "below range" "$OUT"/stdcurve/*/stats_*.csv && grep -q "above range" "$OUT"/stdcurve/*/stats_*.csv; then
  echo "PASS: unknowns interpolated on a 4PL, out-of-range unknowns flagged"
else
  echo "FAIL: expected the 4PL and the range flags"; fail=1
fi

echo "[36] cumulative incidence of competing risks (expect Gray's test)"
./cli/figkit plot --recipe cumulative_incidence --data example/cumulative_incidence.csv \
  --time months --event status --x arm --xlab "Months after transplant" --out "$OUT/cif"
check_bundle "$OUT/cif" "cumulative_incidence"
if compgen -G "$OUT/cif/*/stats_*.csv" > /dev/null && grep -q "gray_p" "$OUT"/cif/*/stats_*.csv && \
   grep -q "Gray's k-sample test" "$OUT"/cif/*/methods_*.md; then
  echo "PASS: cumulative incidence with Gray's test"
else
  echo "FAIL: expected the cumulative incidence and Gray's test"; fail=1
fi

echo "[37] calibration of two risk models (expect Hosmer-Lemeshow)"
./cli/figkit plot --recipe calibration --data example/calibration.csv \
  --y event_1y --x risk_clinical,risk_extended --positive yes --out "$OUT/cal"
check_bundle "$OUT/cal" "calibration"
if compgen -G "$OUT/cal/*/stats_*.csv" > /dev/null && grep -q "hl_p" "$OUT"/cal/*/stats_*.csv && \
   grep -q "Hosmer-Lemeshow" "$OUT"/cal/*/methods_*.md; then
  echo "PASS: calibration with the Hosmer-Lemeshow test"
else
  echo "FAIL: expected the calibration statistics"; fail=1
fi

echo "[38] Sankey diagram of response transitions (expect descriptive, three visits)"
./cli/figkit plot --recipe sankey --data example/sankey.csv \
  --id patient --x week --y response --xlab "Week" --out "$OUT/sankey"
check_bundle "$OUT/sankey" "sankey"
if compgen -G "$OUT/sankey/*/stats_*.csv" > /dev/null && grep -q "from_state" "$OUT"/sankey/*/stats_*.csv && \
   grep -q "Sankey diagram (descriptive)" "$OUT"/sankey/*/manifest_*.json; then
  echo "PASS: Sankey transitions table, descriptive"
else
  echo "FAIL: expected the Sankey transitions table"; fail=1
fi
./cli/figkit plot --recipe sankey --data example/sankey.csv --id patient --x week --y response \
  --levels "progressive disease,stable disease,partial response,complete response" \
  --out "$OUT/sankey_levels" > /dev/null 2>&1
sk_m="$(find "$OUT/sankey_levels" -name 'manifest_*.json' | head -1)"
if [ -n "$sk_m" ] && [ "$(json_array "$sk_m" states | head -1)" = '"progressive disease"' ]; then
  echo "PASS: Sankey --levels sets the state order"
else
  echo "FAIL: expected --levels to put progressive disease first"; fail=1
fi

echo "[39] MA plot of a DESeq2 table (expect baseMean on the log10 scale)"
./cli/figkit plot --recipe ma_plot --data example/de_basemean.csv --out "$OUT/ma_plot"
check_bundle "$OUT/ma_plot" "ma_plot"
if compgen -G "$OUT/ma_plot/*/methods_*.md" > /dev/null && \
   grep -q "MA plot" "$OUT"/ma_plot/*/methods_*.md && grep -q "log10 mean of normalized counts" "$OUT"/ma_plot/*/methods_*.md; then
  echo "PASS: MA plot against log10 baseMean"
else
  echo "FAIL: expected an MA plot against log10 baseMean"; fail=1
fi

echo "[40] Manhattan and QQ plot (expect lambda_GC and lead SNPs)"
./cli/figkit plot --recipe manhattan --data example/gwas_summary.csv \
  --group chr --x pos --y p --label snp --out "$OUT/manhattan"
check_bundle "$OUT/manhattan" "manhattan"
if compgen -G "$OUT/manhattan/*/stats_*.csv" > /dev/null && \
   grep -q "lambda_gc" "$OUT"/manhattan/*/stats_*.csv && grep -q "^lead" "$OUT"/manhattan/*/stats_*.csv; then
  echo "PASS: Manhattan summary and lead SNPs"
else
  echo "FAIL: expected the Manhattan summary and lead SNPs"; fail=1
fi

echo "[41] feature panel of cytokines (expect Kruskal-Wallis with BH)"
./cli/figkit plot --recipe feature_panel --data example/cytokine_panel.csv \
  --x group --y concentration_pg_ml --facet cytokine --out "$OUT/feature_panel"
check_bundle "$OUT/feature_panel" "feature_panel"
if compgen -G "$OUT/feature_panel/*/stats_*.csv" > /dev/null && \
   grep -q "Kruskal-Wallis test" "$OUT"/feature_panel/*/stats_*.csv && grep -q "p_adj" "$OUT"/feature_panel/*/stats_*.csv; then
  echo "PASS: chose Kruskal-Wallis for every feature"
else
  echo "FAIL: expected Kruskal-Wallis with BH-adjusted p"; fail=1
fi

echo "[42] embedding coloured by cluster (expect silhouette widths)"
./cli/figkit plot --recipe embedding --data example/umap_clusters.csv \
  --x UMAP_1 --y UMAP_2 --group cell_type --out "$OUT/embedding"
check_bundle "$OUT/embedding" "embedding"
if compgen -G "$OUT/embedding/*/stats_*.csv" > /dev/null && \
   grep -q "mean_silhouette" "$OUT"/embedding/*/stats_*.csv && grep -q "CD4 T cells" "$OUT"/embedding/*/stats_*.csv; then
  echo "PASS: embedding with per-cluster silhouettes"
else
  echo "FAIL: expected the per-cluster silhouette table"; fail=1
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
pyc volcano --recipe volcano             --data example/de_results.csv
pyc enrich  --recipe enrichment_dot      --data example/enrichment_results.csv
pyc paired  --recipe paired_compare      --data example/paired_response.csv  --x condition --y value --id subject
pyc prop    --recipe proportions         --data example/response_by_arm.csv  --x arm --y response
pyc pca     --recipe pca                 --data example/pca_samples.csv      --group group
pyc rain    --recipe two_group_compare   --data example/tumor_volume.csv     --x group --y volume --geom raincloud
pyc bargeom --recipe multi_group_compare --data example/gene_expression.csv  --x genotype --y expression --geom bar
pyc spider  --recipe spider_response     --data example/spider_response.csv  --x week --y target_lesion_mm --id patient --group arm
pyc waterfall --recipe waterfall         --data example/spider_response.csv  --x week --y target_lesion_mm --id patient --group arm
pyc drc     --recipe dose_response       --data example/dose_response.csv    --x concentration_um --y viability_pct --group compound
pyc roc     --recipe roc_curve           --data example/biomarker_roc.csv    --y status --x marker_a,marker_b --positive disease
pyc ba      --recipe bland_altman        --data example/method_agreement.csv --x new_assay --y reference_assay
pyc meta    --recipe forest_meta         --data example/meta_trials.csv      --y hazard_ratio --ci_low ci_low --ci_high ci_high --label trial
pyc tc      --recipe timecourse          --data example/tumor_growth.csv     --x day --y volume_mm3 --id mouse --group arm
pyc swimmer --recipe swimmer             --data example/swimmer_patients.csv --id patient --time months --event ended --group arm --markers response_month,progression_month
pyc dist    --recipe distribution        --data example/distribution.csv     --x condition --y cell_area_um2
pyc comp    --recipe composition         --data example/composition.csv      --id sample --x cell_type --y n_cells --group condition
pyc faint   --recipe factorial_anova     --data example/twoway_response.csv  --y response --x genotype --fill treatment --geom interaction
pyc estimation --recipe estimation_plot  --data example/estimation_plot.csv  --x treatment --y speed_um_h --ref vehicle
pyc est_paired --recipe estimation_plot  --data example/paired_response.csv  --x condition --y value --id subject --paired --ref pre
pyc enzyme  --recipe enzyme_kinetics     --data example/enzyme_kinetics.csv  --x substrate_um --y rate --group inhibitor
pyc stdcurve --recipe standard_curve     --data example/standard_curve.csv   --x conc_ng_ml --y absorbance --id sample
pyc cif     --recipe cumulative_incidence --data example/cumulative_incidence.csv --time months --event status --x arm
pyc cal     --recipe calibration         --data example/calibration.csv      --y event_1y --x risk_clinical,risk_extended --positive yes
pyc sankey  --recipe sankey              --data example/sankey.csv           --id patient --x week --y response
pyc ma_plot --recipe ma_plot             --data example/de_basemean.csv
pyc manhattan --recipe manhattan         --data example/gwas_summary.csv     --group chr --x pos --y p --label snp
pyc fpanel  --recipe feature_panel       --data example/cytokine_panel.csv   --x group --y concentration_pg_ml --facet cytokine
pyc embed   --recipe embedding           --data example/umap_clusters.csv    --x UMAP_1 --y UMAP_2 --group cell_type

# Engine parity on the non-parametric paths. Both engines must report the same
# p-value and the same effect size to six significant digits.
echo
echo "[parity] R against Python on the non-parametric tests"
col_value() {  # $1 = stats csv, $2 = column name
  awk -F, -v want="$2" 'NR==1 {for (i=1;i<=NF;i++) if ($i==want) c=i; next}
                        NR==2 {printf "%.6g", $c}' "$1"
}
parity_case() {  # $1 = label; rest = plot args
  local label="$1"; shift
  ./cli/figkit plot "$@" --out "$OUT/par_${label}_r" > /dev/null 2>&1
  ./cli/figkit plot --engine python "$@" --out "$OUT/par_${label}_py" > /dev/null 2>&1
  local rf pf
  rf="$(find "$OUT/par_${label}_r" -name 'stats_*.csv' | head -1)"
  pf="$(find "$OUT/par_${label}_py" -name 'stats_*.csv' | head -1)"
  if [ -z "$rf" ] || [ -z "$pf" ]; then
    echo "FAIL: parity $label wrote no stats table"; fail=1; return
  fi
  local col ok=1
  for col in p_value effect_size; do
    if [ "$(col_value "$rf" "$col")" != "$(col_value "$pf" "$col")" ]; then
      echo "  MISMATCH: $col is $(col_value "$rf" "$col") in R and $(col_value "$pf" "$col") in Python"
      ok=0; fail=1
    fi
  done
  [ "$ok" -eq 1 ] && echo "PASS: parity $label"
}
parity_case mwu --recipe two_group_compare --data example/cytokine_pg_ml.csv \
  --x group --y concentration
parity_case signed_rank --recipe paired_compare --data example/paired_response.csv \
  --x condition --y value --id subject --test wilcoxon
parity_case composition --recipe composition --data example/composition.csv \
  --id sample --x cell_type --y n_cells --group condition

# The estimation plot bootstrap draws the same Park-Miller resamples in both
# engines, so the BCa interval must agree, not only the t-test.
echo "[parity] R against Python on the bootstrap CI of estimation_plot"
./cli/figkit plot --recipe estimation_plot --data example/estimation_plot.csv \
  --x treatment --y speed_um_h --ref vehicle --out "$OUT/par_est_r" > /dev/null 2>&1
./cli/figkit plot --engine python --recipe estimation_plot --data example/estimation_plot.csv \
  --x treatment --y speed_um_h --ref vehicle --out "$OUT/par_est_py" > /dev/null 2>&1
est_rf="$(find "$OUT/par_est_r" -name 'stats_*.csv' | head -1)"
est_pf="$(find "$OUT/par_est_py" -name 'stats_*.csv' | head -1)"
if [ -z "$est_rf" ] || [ -z "$est_pf" ]; then
  echo "FAIL: parity estimation_plot wrote no stats table"; fail=1
else
  est_ok=1
  for col in mean_diff ci_low ci_high hedges_g p_value; do
    if [ "$(col_value "$est_rf" "$col")" != "$(col_value "$est_pf" "$col")" ]; then
      echo "  MISMATCH: $col is $(col_value "$est_rf" "$col") in R and $(col_value "$est_pf" "$col") in Python"
      est_ok=0; fail=1
    fi
  done
  [ "$est_ok" -eq 1 ] && echo "PASS: parity estimation_plot bootstrap CI"
fi

# Both engines build the same tree and order each branch by the mean of the
# unscaled values (heatmap) or by the mean correlation (correlation_heatmap), so
# the rows and the columns come in one order.
order_parity() {  # $1 = label; rest = plot args
  local label="$1"; shift
  ./cli/figkit plot "$@" --out "$OUT/par_${label}_r" > /dev/null 2>&1
  ./cli/figkit plot --engine python "$@" --out "$OUT/par_${label}_py" > /dev/null 2>&1
  local rf pf key ok=1
  rf="$(find "$OUT/par_${label}_r" -name 'manifest_*.json' | head -1)"
  pf="$(find "$OUT/par_${label}_py" -name 'manifest_*.json' | head -1)"
  if [ -z "$rf" ] || [ -z "$pf" ]; then
    echo "FAIL: parity $label wrote no manifest"; fail=1; return
  fi
  for key in row_order column_order; do
    if [ -z "$(json_array "$rf" "$key")" ] || \
       [ "$(json_array "$rf" "$key")" != "$(json_array "$pf" "$key")" ]; then
      echo "  MISMATCH: $label $key differs between R and Python"
      ok=0; fail=1
    fi
  done
  [ "$ok" -eq 1 ] && echo "PASS: parity $label row and column order"
}
echo "[parity] R against Python on the heatmap row and column order"
order_parity heatmap --recipe heatmap --data example/expression_matrix.csv \
  --annotation example/expression_annotation.csv
order_parity correlation_heatmap --recipe correlation_heatmap --data example/correlation_vars.csv

# Area B, core: the read, bundle and command layer.
echo "[B1] an option that a recipe does not read stops the run before any output (F71)"
for engine in r python; do
  rc=0
  msg="$(./cli/figkit plot --engine "$engine" --recipe survival_km --data example/survival_trial.csv \
    --time time --event event --x arm --title "KM" --out "$OUT/b1_$engine" 2>&1)" || rc=$?
  if [ "$rc" -ne 0 ] && [ ! -d "$OUT/b1_$engine" ] && \
     printf '%s\n' "$msg" | grep -q -- "--title is not an option of survival_km; it is read by:"; then
    echo "PASS: $engine stops on --title for survival_km"
  else
    echo "FAIL: $engine did not stop on --title for survival_km"; fail=1
  fi
done

# Area C: the clinical recipes.
echo
echo "[C1] survival_km reads column names that hold a space, in both engines"
sed '1s/.*/"follow up months","death event","treatment arm","age","sex","stage"/' \
  example/survival_trial.csv > "$OUT/km_space_names.csv"
for eng in r python; do
  ./cli/figkit plot --engine "$eng" --recipe survival_km --data "$OUT/km_space_names.csv" \
    --time "follow up months" --event "death event" --x "treatment arm" \
    --out "$OUT/km_space_$eng" > /dev/null 2>&1
  check_bundle "$OUT/km_space_$eng" "survival_km with spaced names ($eng)"
done

# Area A: the comparison recipes.
echo
echo "[A1] column names with a space run in both engines (F44)"
SP="$OUT/space_cols"; mkdir -p "$SP"
sed '1s/.*/"treatment group","tumor volume"/' example/tumor_volume.csv > "$SP/tumor_space.csv"
sed '1s/.*/"mouse genotype","expression level"/' example/gene_expression.csv > "$SP/gene_space.csv"
sed '1s/.*/"mouse genotype","drug treatment","response value"/' example/twoway_response.csv \
  > "$SP/twoway_space.csv"
for engine in r python; do
  ./cli/figkit plot --engine "$engine" --recipe two_group_compare --data "$SP/tumor_space.csv" \
    --x "treatment group" --y "tumor volume" --out "$SP/tg_$engine" > /dev/null 2>&1
  check_bundle "$SP/tg_$engine" "two_group_compare with spaces in the names ($engine)"
  ./cli/figkit plot --engine "$engine" --recipe multi_group_compare --data "$SP/gene_space.csv" \
    --x "mouse genotype" --y "expression level" --test kruskal --out "$SP/mg_$engine" > /dev/null 2>&1
  check_bundle "$SP/mg_$engine" "multi_group_compare with spaces in the names ($engine)"
done
./cli/figkit plot --recipe factorial_anova --data "$SP/twoway_space.csv" --x "mouse genotype" \
  --fill "drug treatment" --y "response value" --test art --out "$SP/fa_art_r" > /dev/null 2>&1
check_bundle "$SP/fa_art_r" "factorial_anova --test art with spaces in the names (r)"

echo
if [ "$fail" -eq 0 ]; then
  echo "ALL SMOKE TESTS PASSED"
else
  echo "SMOKE TESTS FAILED"; exit 1
fi
