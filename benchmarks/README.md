# Benchmarks

This folder holds public benchmark datasets and the reference values that each
recipe must reproduce. Each reference value comes from an independent call to
the canonical R function, for example `stats::t.test`, `survival::coxph`,
`pROC::roc` or `metafor::rma`. A pubready recipe never computes a reference
value. A match is therefore evidence that the recipe computes the correct number,
in each engine.

## Contents

| Path | Contents |
|---|---|
| `data/*.csv` | The public datasets, one CSV file for each dataset. |
| `reference_values.csv` | One row for each checked quantity: the case, the quantity, the value, the relative tolerance and the function that computed the value. |

## How to run the benchmarks

1. Run `./cli/figkit test`. The command mounts this folder, and the R suite
   (`test-benchmarks.R`) and the Python suite (`test_benchmarks.py`) compare
   each recipe with `reference_values.csv`.
2. Run `./tests/benchmark_test.sh`. The script runs every case through
   `figkit plot` with both engines and asserts that each run writes a complete
   bundle. Give a case name as the first argument to run only the matching cases.

To export the datasets and compute the reference values again, run the exporter
in the pinned image from the repository root:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work localhost/pubready:<version> \
  Rscript tools/make_benchmarks/make_benchmark_data.R
```

## Cases

| Case | File | Recipe | Reference function |
|---|---|---|---|
| `sleep_paired` | `sleep.csv` | `paired_compare --test t` | `stats::t.test(paired = TRUE)` |
| `sleep_welch` | `sleep.csv` | `two_group_compare --test welch_t` | `stats::t.test(var.equal = FALSE)` |
| `plant_growth_anova` | `plant_growth.csv` | `multi_group_compare` | `stats::anova(lm())`, `stats::TukeyHSD` |
| `insect_sprays_kruskal` | `insect_sprays.csv` | `multi_group_compare --test kruskal` | `stats::kruskal.test` |
| `tooth_growth_twoway` | `tooth_growth.csv` | `factorial_anova --test anova` | `car::Anova(type = 2)` |
| `lung_km` | `lung.csv` | `survival_km` | `survival::survdiff` |
| `lung_cox` | `lung.csv` | `cox_forest` | `survival::coxph`, `survival::cox.zph` |
| `mtcars_correlation` | `mtcars.csv` | `correlation --test pearson` | `stats::cor.test` |
| `mtcars_correlation_heatmap` | `mtcars.csv` | `correlation_heatmap` | `stats::cor` |
| `mtcars_heatmap` | `mtcars_matrix.csv`, `mtcars_annotation.csv` | `heatmap` | `stats::t.test` for each feature |
| `airquality_spearman` | `airquality.csv` | `correlation --test spearman` | `stats::cor.test(method = "spearman")` |
| `iris_pca` | `iris.csv` | `pca` | `stats::prcomp(scale. = TRUE)` |
| `titanic_proportions` | `titanic.csv` | `proportions --test chisq` | `stats::chisq.test` |
| `ucb_proportions` | `ucb_admissions.csv` | `proportions --test chisq` | `stats::chisq.test`, cross-product ratio |
| `dnase_4pl` | `dnase_run1.csv` | `dose_response` | `stats::nls` with `SSfpl` |
| `dnase_shared` | `dnase.csv` | `dose_response --group Run` | `stats::nls` with `SSfpl`, extra sum of squares |
| `asah_roc` | `asah.csv` | `roc_curve` | `pROC::roc`, `pROC::ci.auc`, `pROC::roc.test` |
| `pefr_bland_altman` | `pefr.csv` | `bland_altman` | the published bias and SD, the mean and SD of the differences |
| `bcg_meta` | `bcg.csv` | `forest_meta` | `metafor::rma(method = "DL")` and `"FE"` |
| `orthodont_timecourse` | `orthodont.csv` | `timecourse` | `stats::aov` with `Error()`, `car::Anova(idesign)` |
| `chick_weight_timecourse` | `chick_weight.csv` | `timecourse` | `stats::aov` with `Error()`, `car::Anova(idesign)` |
| `aml_swimmer` | `aml.csv` | `swimmer` | `stats::median` and the count of censored courses for each arm |
| `airway_volcano` | `airway_de.csv` | `volcano --y padj` | the DESeq2 genes at padj of 0.05 or less and an absolute log2FC of 1 or more |
| `airway_heatmap` | `airway_heatmap.csv`, `airway_annotation.csv` | `heatmap` | `stats::t.test` for each gene |
| `airway_pca` | `airway_pca.csv` | `pca --group dex` | `stats::prcomp(scale. = TRUE)` |
| `airway_enrichment` | `airway_enrichment.csv` | `enrichment_dot` | `stats::phyper`, `stats::p.adjust(method = "BH")` |
| `airway_upset` | `airway_upset.csv` | `upset` | `base::table` of the membership patterns |
| `anorexia_paired` | `anorexia_ft.csv` | `paired_compare --test t` | `stats::t.test(paired = TRUE)` |
| `veteran_km` | `veteran.csv` | `survival_km --x celltype` | `survival::survdiff` |
| `veteran_cox` | `veteran.csv` | `cox_forest` | `survival::coxph`, `survival::cox.zph` |
| `state_correlation_heatmap` | `state_x77.csv` | `correlation_heatmap` | `stats::cor`, `stats::cor.test`, `stats::p.adjust(method = "BH")` |
| `faithful_distribution` | `faithful.csv` | `distribution --y eruptions` | `stats::shapiro.test`, `mean`, `stats::sd`, `stats::bw.nrd0`, `stats::IQR` |
| `tooth_growth_distribution` | `tooth_growth.csv` | `distribution --x supp --y len` | `stats::shapiro.test`, `mean`, `stats::sd`, `stats::bw.nrd0` for each supplement |
| `housing_composition` | `housing.csv` | `composition --id sample --x Sat --y Freq --group Cont` | `stats::wilcox.test` for each satisfaction level, `stats::p.adjust(method = "BH")` |
| `chickwts_estimation` | `chickwts.csv` | `estimation_plot --ref casein` | the mean difference, `stats::t.test(var.equal = FALSE)` for each feed |
| `puromycin_kinetics` | `puromycin.csv` | `enzyme_kinetics --group state` | `stats::nls` for each state, `stats::anova` of the shared-K and the shared-Vm model |
| `dnase_standard_curve` | `dnase_standard_curve.csv` | `standard_curve --id sample` | the inverse of the `stats::nls` fit with `SSfpl` to the run 1 standards |
| `mgus2_cif` | `mgus2.csv` | `cumulative_incidence --x sex` | `cmprsk::cuminc` (Gray's test), `cmprsk::timepoints` |
| `pima_calibration` | `pima_te.csv` | `calibration --positive Yes` | `ResourceSelection::hoslem.test(g = 10)`, `stats::glm` |
| `pbcseq_sankey` | `pbcseq_edema.csv` | `sankey --x year --y edema` | `stats::mcnemar.test(correct = FALSE)` |
| `airway_ma` | `airway_de.csv` | `ma_plot` | the DESeq2 genes at padj of 0.05 or less and an absolute log2FC of 1 or more |
| `gwas_manhattan` | `gwas_results.csv` | `manhattan --group CHR --x BP --y P --label SNP` | the median of `qchisq(1 - p, 1)` over `qchisq(0.5, 1)`, the counts below 5e-8 and 1e-5, the smallest p |
| `iris_feature_panel_kruskal` | `iris.csv` | `feature_panel --x Species` | `stats::kruskal.test` for each measurement, `stats::p.adjust(method = "BH")` |
| `iris_feature_panel_anova` | `iris.csv` | `feature_panel --x Species --test anova` | `stats::oneway.test(var.equal = TRUE)` for each measurement, `stats::p.adjust(method = "BH")` |
| `iris_embedding` | `iris.csv` | `embedding --x Petal.Length --y Petal.Width --group Species` | `cluster::silhouette` on the two petal measurements |

## Tolerances

Most reference values carry a relative tolerance of 1e-6. The cases below are
wider, and the reason follows:

- **The 4PL parameters of `dnase_4pl` and `dnase_shared`.** Two optimisers stop
  at their own convergence criteria on a flat optimum. The fitted parameters
  therefore carry 1e-4. The bottom plateau is near zero, so it carries 1e-2. The
  residual sum of squares stays at 1e-6.
- **The published PEFR values.** The paper rounds the bias to -2.1 and the SD to
  38.8, so these two checks carry the tolerance of the rounding.
- **The interpolated concentrations of `dnase_standard_curve`.** Each value is
  read off a 4PL fit, so it carries the 1e-4 of the `dnase_4pl` parameters. The
  counts of the unknowns in range, below range and above range stay at 1e-6.

## Data preparation

The exporter keeps each dataset as its package distributes it. The changes
below make a dataset readable by a recipe. Each change is in
`tools/make_benchmarks/make_benchmark_data.R`.

- **`lung.csv`.** The survival package codes `status` as 1 = censored and 2 =
  dead. pubready reads 1 = event, so the exporter adds `event = status - 1`. The
  exporter also writes `sex` as `male` and `female`.
- **`mtcars_matrix.csv` and `mtcars_annotation.csv`.** These files hold six
  variables of `mtcars` as a features-by-samples matrix for the heatmap recipe.
  The annotation gives the transmission of each car.
- **`titanic.csv` and `ucb_admissions.csv`.** These files expand the
  contingency tables to one row for each person.
- **`dnase_run1.csv`.** This file holds run 1 of `DNase`, which is the example of
  the `SSfpl` help page.
- **`pefr.csv`.** This file is a transcription of Table 1 of Bland and Altman
  (1986), checked against the table on the page of Martin Bland at
  <https://www-users.york.ac.uk/~mb55/meas/ba.htm>.
- **`bcg.csv`.** The exporter converts the trial counts to a risk ratio and its
  95% CI with `metafor::escalc`, because a paper reports a trial in that form.
- **`chick_weight.csv`.** Five chicks died before day 21. The file keeps them,
  and the `timecourse` recipe drops them, because a repeated-measures ANOVA needs
  complete subjects.
- **`airway_de.csv`.** The exporter runs DESeq2 on the `airway` data with the
  design `~ cell + dex`, so the test blocks the donor. The file keeps the 12,625
  genes that have an adjusted p and a mean count of 50 or more. This limit keeps
  the file at 544 KB. The file holds the `baseMean` of each gene, so the
  `ma_plot` case reads the same table as the `volcano` case.
- **`airway_heatmap.csv` and `airway_pca.csv`.** These files hold the
  variance-stabilised counts of `DESeq2::vst`. The heatmap file holds the 40 genes
  with the smallest adjusted p and an absolute log2FC of 1 or more. The PCA file
  holds the 500 most variable genes, one row for each sample.
- **`airway_enrichment.csv`.** The exporter tests each GO biological-process
  term of 15 to 500 genes for over-representation of the 596 annotated DE genes,
  with the hypergeometric test and the Benjamini-Hochberg adjustment. The file
  keeps the 30 strongest terms in the layout of clusterProfiler.
- **`airway_upset.csv`.** This file holds the DE genes of five enriched terms.
  The exporter goes down the ranking and skips a term that overlaps a chosen term
  with a Jaccard index of 0.5 or more, so the five sets are distinct.
- **`anorexia_ft.csv`.** This file holds the 17 patients of the family-therapy
  group of `MASS::anorexia`, in long form with a subject id. The weights are in
  pounds. The original source states kilograms in error.
- **`veteran.csv`.** The survival package codes `status` as 1 = dead, which is
  the pubready code. The exporter writes `trt` as `standard` and `test`.
- **`state_x77.csv`.** The file keeps the column names `Life Exp` and `HS Grad`
  as R distributes them. The case therefore also tests a column name with a space.
- **`faithful.csv`.** The file keeps the 272 eruptions of `datasets::faithful`
  as R distributes them. The `distribution` case uses the bimodal `eruptions`
  column.
- **`tooth_growth.csv`.** The `tooth_growth_distribution` case reads the file of
  the factorial case. It describes `len` within each level of `supp`.
- **`housing.csv`.** The exporter writes the 72 rows of `MASS::housing` with a
  `sample` column. A sample is one combination of `Infl`, `Type` and `Cont`, so
  the file holds 24 samples. The categories are the `Sat` levels, the counts are
  `Freq`, and the group is `Cont`. The reference proportions are `Freq` divided
  by the sample total.
- **`chickwts.csv`.** The file keeps the `feed` and `weight` columns. Casein sorts
  first, but the command names it with `--ref casein`, so the case also tests
  the option.
- **`puromycin.csv`.** The file keeps `conc`, `rate` and `state`. The reference
  fits use a convergence tolerance of 1e-8, so each reference value is the
  optimum and not a point near it.
- **`dnase_standard_curve.csv`.** Run 1 of `DNase` gives the 16 standards. The
  16 densities of run 2, with an empty concentration, are the unknowns. The
  `sample` column names each row.
- **`mgus2.csv`.** The file follows the survival vignette "compete". The time is
  the time to progression for a patient who progressed, else the follow-up time.
  The event is 1 for a progression to a plasma cell malignancy, 2 for a death
  without progression and 0 for a censored patient.
- **`pima_te.csv`.** The exporter fits the logistic model
  `type ~ npreg + glu + bmi + ped + age` to `MASS::Pima.tr`. The file holds
  `MASS::Pima.te` with the predicted risk of each woman in the column `predicted`.
- **`pbcseq_edema.csv`.** The file holds the edema score of each patient at the
  first visit, as `year = 0`, and at the visit closest to day 730 among the visits
  from day 365 to day 1095, as `year = 2`. A patient without such a visit has only
  the first row, so the test uses the 249 patients with both visits.
- **`gwas_results.csv`.** This file is the `gwasResults` data of the qqman
  package as the package distributes it. The positions count the SNPs on each
  chromosome, because the data is simulated.
- **`iris.csv` for `feature_panel` and `embedding`.** The feature panel reads
  the file in wide form, so each of the four measurements is a feature. The
  embedding case uses the petal length and the petal width as the two
  coordinates. A silhouette on these two columns describes the separation of the
  species on two measurements and not on a computed embedding.

## Sources

Each dataset comes from the R package or the publication in this table. The
license of that package applies to the data.

| File | Package | Source |
|---|---|---|
| `sleep.csv` | datasets | Cushny, A. R. and Peebles, A. R. (1905) The action of optical isomers: II hyoscines. The Journal of Physiology 32, 501-510. Student (1908) The probable error of the mean. Biometrika 6, 20. |
| `plant_growth.csv` | datasets | Dobson, A. J. (1983) An Introduction to Statistical Modelling. London: Chapman and Hall. |
| `insect_sprays.csv` | datasets | Beall, G. (1942) The transformation of data from entomological field experiments. Biometrika 29, 243-262. |
| `tooth_growth.csv` | datasets | Bliss, C. I. (1952) The Statistics of Bioassay. Academic Press. Crampton, E. W. (1947) The Journal of Nutrition 33(5), 491-504. |
| `lung.csv` | survival | Loprinzi, C. L. et al. (1994) Prospective evaluation of prognostic variables from patient-completed questionnaires. Journal of Clinical Oncology 12(3), 601-607. |
| `mtcars.csv` | datasets | Henderson and Velleman (1981) Building multiple regression models interactively. Biometrics 37, 391-411. |
| `airquality.csv` | datasets | New York State Department of Conservation and the National Weather Service. Chambers, J. M. et al. (1983) Graphical Methods for Data Analysis. Wadsworth. |
| `iris.csv` | datasets | Anderson, E. (1935) The irises of the Gaspe Peninsula. Bulletin of the American Iris Society 59, 2-5. Fisher, R. A. (1936) Annals of Eugenics 7, 179-188. |
| `titanic.csv` | datasets | Dawson, R. J. MacG. (1995) The "Unusual Episode" data revisited. Journal of Statistics Education 3. |
| `ucb_admissions.csv` | datasets | Bickel, P. J., Hammel, E. A. and O'Connell, J. W. (1975) Sex bias in graduate admissions: data from Berkeley. Science 187, 398-403. |
| `dnase.csv`, `dnase_run1.csv` | datasets | Davidian, M. and Giltinan, D. M. (1995) Nonlinear Models for Repeated Measurement Data. Chapman and Hall, section 5.2.4. |
| `asah.csv` | pROC | Turck, N. et al. (2010) A multiparameter panel method for outcome prediction following aneurysmal subarachnoid hemorrhage. Intensive Care Medicine 36(1), 107-115. |
| `pefr.csv` | none | Bland, J. M. and Altman, D. G. (1986) Statistical methods for assessing agreement between two methods of clinical measurement. Lancet 327, 307-310. |
| `bcg.csv` | metadat | Colditz, G. A. et al. (1994) Efficacy of BCG vaccine in the prevention of tuberculosis: meta-analysis of the published literature. JAMA 271(9), 698-702. |
| `orthodont.csv` | nlme | Potthoff, R. F. and Roy, S. N. (1964) A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika 51, 313-326. |
| `chick_weight.csv` | datasets | Crowder, M. and Hand, D. (1990) Analysis of Repeated Measures. Chapman and Hall, example 5.3. |
| `aml.csv` | survival | Miller, R. G. (1997) Survival Analysis. John Wiley and Sons. |
| `airway_*.csv` | airway, DESeq2, org.Hs.eg.db, GO.db | Himes, B. E. et al. (2014) RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells. PLoS ONE 9(6), e99625. GEO GSE52778. The gene symbols and the GO annotation come from org.Hs.eg.db 3.20.0 and GO.db 3.20.0. |
| `anorexia_ft.csv` | MASS | Hand, D. J., Daly, F., McConway, K., Lunn, D. and Ostrowski, E. eds (1993) A Handbook of Small Data Sets. Chapman and Hall, data set 285. |
| `veteran.csv` | survival | Kalbfleisch, J. D. and Prentice, R. L. (1980) The Statistical Analysis of Failure Time Data. Wiley. |
| `state_x77.csv` | datasets | U.S. Department of Commerce, Bureau of the Census (1977) Statistical Abstract of the United States. |
| `faithful.csv` | datasets | Azzalini, A. and Bowman, A. W. (1990) A look at some data on the Old Faithful geyser. Applied Statistics 39, 357-365. |
| `housing.csv` | MASS | Madsen, M. (1976) Statistical analysis of multiple contingency tables. Two examples. Scandinavian Journal of Statistics 3, 97-106. |
| `chickwts.csv` | datasets | Anonymous (1948) Biometrika 35, 214. McNeil, D. R. (1977) Interactive Data Analysis. Wiley. |
| `puromycin.csv` | datasets | Treloar, M. A. (1974) Effects of Puromycin on Galactosyltransferase in Golgi Membranes. M.Sc. Thesis, University of Toronto. Bates, D. M. and Watts, D. G. (1988) Nonlinear Regression Analysis and Its Applications. Wiley. |
| `dnase_standard_curve.csv` | datasets | Davidian, M. and Giltinan, D. M. (1995) Nonlinear Models for Repeated Measurement Data. Chapman and Hall, section 5.2.4. |
| `mgus2.csv` | survival | Kyle, R. A. et al. (2002) A long-term study of prognosis in monoclonal gammopathy of undetermined significance. New England Journal of Medicine 346, 564-569. |
| `pima_te.csv` | MASS | Smith, J. W. et al. (1988) Using the ADAP learning algorithm to forecast the onset of diabetes mellitus. Proceedings of the Annual Symposium on Computer Application in Medical Care, 261-265. |
| `pbcseq_edema.csv` | survival | Murtaugh, P. A. et al. (1994) Primary biliary cirrhosis: prediction of short-term survival based on repeated patient visits. Hepatology 20(1), 126-134. |
| `gwas_results.csv` | qqman | Turner, S. D. (2018) qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software 3(25), 731. |
