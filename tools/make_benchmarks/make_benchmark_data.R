#!/usr/bin/env Rscript
# Export the public benchmark datasets to benchmarks/data/ and compute the
# reference values that each recipe must reproduce. Every reference value comes
# from an independent call to the canonical R function (t.test, aov, survdiff,
# coxph, prcomp, nls with SSfpl, pROC, metafor), never from a pubready recipe, so
# a match is evidence that the recipe computes the right number.
#
# Run inside the pinned image from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:<version> \
#     Rscript tools/make_benchmarks/make_benchmark_data.R
suppressPackageStartupMessages({
  library(survival)
  library(pROC)
  library(metafor)
  library(DESeq2)
  library(org.Hs.eg.db)
})

out_dir <- file.path("benchmarks", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

refs <- list()

# Append one reference value. rel_tol is the relative tolerance the engine tests
# apply; a published value rounded in the paper gets a wider tolerance.
add_ref <- function(case, quantity, value, source, rel_tol = 1e-6) {
  refs[[length(refs) + 1]] <<- data.frame(
    case = case, quantity = quantity, value = signif(value, 12),
    rel_tol = rel_tol, source = source, stringsAsFactors = FALSE)
}

write_data <- function(df, name) {
  utils::write.csv(df, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  # Read back, so each reference is computed on exactly the bytes the recipes read.
  utils::read.csv(file.path(out_dir, paste0(name, ".csv")), check.names = FALSE,
                  stringsAsFactors = FALSE)
}

# Jaccard index between a candidate gene set and an already chosen set.
jaccard_with <- function(id, genes, members) {
  other <- members[[id]]
  length(intersect(genes, other)) / length(union(genes, other))
}

# Expand a contingency table into one row per observation.
expand_table <- function(tab) {
  d <- as.data.frame(tab, stringsAsFactors = FALSE)
  d[rep(seq_len(nrow(d)), d$Freq), setdiff(names(d), "Freq"), drop = FALSE]
}

# 1. sleep: Cushny and Peebles (1905), the data of Student (1908).
sleep_df <- write_data(datasets::sleep[, c("ID", "group", "extra")], "sleep")
g1 <- sleep_df$extra[sleep_df$group == 1]
g2 <- sleep_df$extra[sleep_df$group == 2]
tt <- stats::t.test(g2, g1, paired = TRUE)
add_ref("sleep_paired", "statistic", tt$statistic, "stats::t.test(paired = TRUE)")
add_ref("sleep_paired", "df", tt$parameter, "stats::t.test(paired = TRUE)")
add_ref("sleep_paired", "p_value", tt$p.value, "stats::t.test(paired = TRUE)")
tw <- stats::t.test(g1, g2, var.equal = FALSE)
add_ref("sleep_welch", "statistic", tw$statistic, "stats::t.test(var.equal = FALSE)")
add_ref("sleep_welch", "df", tw$parameter, "stats::t.test(var.equal = FALSE)")
add_ref("sleep_welch", "p_value", tw$p.value, "stats::t.test(var.equal = FALSE)")

# 2. PlantGrowth: Dobson (1983), one-way design with three groups.
pg <- write_data(datasets::PlantGrowth, "plant_growth")
pa <- stats::anova(stats::lm(weight ~ group, data = pg))
add_ref("plant_growth_anova", "statistic", pa$`F value`[1], "stats::anova(lm())")
add_ref("plant_growth_anova", "p_value", pa$`Pr(>F)`[1], "stats::anova(lm())")
tk <- stats::TukeyHSD(stats::aov(weight ~ group, data = pg))$group
add_ref("plant_growth_anova", "tukey_p_trt1_trt2", tk["trt2-trt1", "p adj"], "stats::TukeyHSD()")

# 3. InsectSprays: Beall (1942), counts with unequal variances.
ins <- write_data(datasets::InsectSprays, "insect_sprays")
kt <- stats::kruskal.test(count ~ spray, data = ins)
add_ref("insect_sprays_kruskal", "statistic", kt$statistic, "stats::kruskal.test()")
add_ref("insect_sprays_kruskal", "p_value", kt$p.value, "stats::kruskal.test()")

# 4. ToothGrowth: Bliss (1952) and Crampton (1947), a balanced 2 x 3 factorial design.
tg <- write_data(datasets::ToothGrowth, "tooth_growth")
tg$dose <- factor(tg$dose); tg$supp <- factor(tg$supp)
ta <- car::Anova(stats::lm(len ~ dose * supp, data = tg), type = 2)
for (term in c("dose", "supp", "dose:supp")) {
  add_ref("tooth_growth_twoway", paste0("F_", term), ta[term, "F value"], "car::Anova(type = 2)")
  add_ref("tooth_growth_twoway", paste0("p_", term), ta[term, "Pr(>F)"], "car::Anova(type = 2)")
}

# 5. lung: the NCCTG advanced lung cancer trial, Loprinzi et al. (1994).
# survival codes status 1 = censored, 2 = dead; pubready reads 1 = event.
lung <- survival::lung
lung$event <- as.integer(lung$status == 2)
lung$sex <- ifelse(lung$sex == 1, "male", "female")
lung_df <- write_data(lung[, c("inst", "time", "event", "age", "sex", "ph.ecog",
                               "ph.karno", "wt.loss")], "lung")
sd <- survival::survdiff(Surv(time, event) ~ sex, data = lung_df)
add_ref("lung_km", "log_rank_chisq", sd$chisq, "survival::survdiff()")
add_ref("lung_km", "log_rank_p",
        stats::pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE), "survival::survdiff()")
cx <- survival::coxph(Surv(time, event) ~ age + sex + ph.ecog, data = lung_df)
cs <- summary(cx)$coefficients
for (term in rownames(cs)) {
  add_ref("lung_cox", paste0("HR_", term), cs[term, "exp(coef)"], "survival::coxph()")
  add_ref("lung_cox", paste0("p_", term), cs[term, "Pr(>|z|)"], "survival::coxph()")
}

# 6. mtcars: Henderson and Velleman (1981), Motor Trend 1974.
mt <- cbind(car = rownames(datasets::mtcars), datasets::mtcars)
mt_df <- write_data(mt, "mtcars")
ct <- stats::cor.test(mt_df$wt, mt_df$mpg, method = "pearson")
add_ref("mtcars_correlation", "estimate", ct$estimate, "stats::cor.test(pearson)")
add_ref("mtcars_correlation", "p_value", ct$p.value, "stats::cor.test(pearson)")
cm <- stats::cor(mt_df[, -1])
add_ref("mtcars_correlation_heatmap", "r_mpg_wt", cm["mpg", "wt"], "stats::cor()")
add_ref("mtcars_correlation_heatmap", "r_disp_cyl", cm["disp", "cyl"], "stats::cor()")

# The same cars as a features-by-samples matrix for the heatmap recipe, with the
# transmission as the sample annotation.
feats <- c("mpg", "disp", "hp", "drat", "wt", "qsec")
mmat <- data.frame(feature = feats, t(as.matrix(datasets::mtcars[, feats])),
                   check.names = FALSE)
mmat_df <- write_data(mmat, "mtcars_matrix")
mann <- data.frame(car = rownames(datasets::mtcars),
                   transmission = ifelse(datasets::mtcars$am == 1, "manual", "automatic"))
mann_df <- write_data(mann, "mtcars_annotation")
for (f in feats) {
  v <- as.numeric(mmat_df[mmat_df$feature == f, -1])
  g <- mann_df$transmission[match(names(mmat_df)[-1], mann_df$car)]
  add_ref("mtcars_heatmap", paste0("p_", f),
          stats::t.test(v[g == "automatic"], v[g == "manual"])$p.value,
          "stats::t.test(var.equal = FALSE)")
}

# 7. airquality: New York, May to September 1973, with missing ozone values.
aq <- write_data(datasets::airquality, "airquality")
cs2 <- suppressWarnings(stats::cor.test(aq$Temp, aq$Ozone, method = "spearman", exact = FALSE))
add_ref("airquality_spearman", "estimate", cs2$estimate, "stats::cor.test(spearman)")
add_ref("airquality_spearman", "n", sum(stats::complete.cases(aq$Temp, aq$Ozone)), "stats::complete.cases()")

# 8. iris: Anderson (1935), Fisher (1936).
ir <- write_data(datasets::iris, "iris")
pc <- stats::prcomp(ir[, 1:4], center = TRUE, scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2) * 100
add_ref("iris_pca", "pc1_variance_pct", ve[1], "stats::prcomp(scale. = TRUE)")
add_ref("iris_pca", "pc2_variance_pct", ve[2], "stats::prcomp(scale. = TRUE)")

# 9. Titanic: passenger class and survival, Dawson (1995).
ti <- write_data(expand_table(datasets::Titanic), "titanic")
tct <- stats::chisq.test(table(ti$Class, ti$Survived), correct = FALSE)
add_ref("titanic_proportions", "statistic", tct$statistic, "stats::chisq.test(correct = FALSE)")
add_ref("titanic_proportions", "p_value", tct$p.value, "stats::chisq.test(correct = FALSE)")

# 10. UCBAdmissions: Bickel et al. (1975), a 2 x 2 table of gender by admission.
ucb <- write_data(expand_table(datasets::UCBAdmissions), "ucb_admissions")
utab <- table(factor(ucb$Gender, levels = unique(ucb$Gender)),
              factor(ucb$Admit, levels = unique(ucb$Admit)))
uct <- stats::chisq.test(utab, correct = FALSE)
add_ref("ucb_proportions", "statistic", uct$statistic, "stats::chisq.test(correct = FALSE)")
add_ref("ucb_proportions", "p_value", uct$p.value, "stats::chisq.test(correct = FALSE)")
add_ref("ucb_proportions", "odds_ratio",
        (utab[1, 1] * utab[2, 2]) / (utab[1, 2] * utab[2, 1]), "cross-product ratio")

# 11. DNase: ELISA of recombinant DNase in rat serum, Davidian and Giltinan (1995).
dn <- as.data.frame(datasets::DNase)
dn$Run <- as.integer(as.character(dn$Run))
dn <- dn[order(dn$Run, dn$conc), ]
dn_df <- write_data(dn, "dnase")
dn1_df <- write_data(dn[dn$Run == 1, ], "dnase_run1")
# SSfpl is the self-starting four-parameter logistic of the stats package. Its
# xmid is the natural log of the EC50 and its scal is the inverse Hill slope.
# Two optimisers stop at their own tolerance on a flat optimum, so the fitted
# parameters carry a relative tolerance of 1e-4, and the bottom plateau, which
# sits near zero, 1e-2. The residual sum of squares is the sharp check.
fpl <- stats::nls(density ~ SSfpl(log(conc), A, B, xmid, scal), data = dn1_df)
cf <- stats::coef(fpl)
add_ref("dnase_4pl", "bottom", cf[["A"]], "stats::nls(SSfpl)", rel_tol = 1e-2)
add_ref("dnase_4pl", "top", cf[["B"]], "stats::nls(SSfpl)", rel_tol = 1e-4)
add_ref("dnase_4pl", "ec50", exp(cf[["xmid"]]), "stats::nls(SSfpl)", rel_tol = 1e-4)
add_ref("dnase_4pl", "hill", 1 / cf[["scal"]], "stats::nls(SSfpl)", rel_tol = 1e-4)
add_ref("dnase_4pl", "rss", sum(stats::residuals(fpl)^2), "stats::nls(SSfpl)")

# The shared-EC50 F test over the eleven runs, with the SSfpl parameterisation:
# one natural-log midpoint for all runs, and the plateaus and scales per run.
dn_df$run_i <- match(dn_df$Run, sort(unique(dn_df$Run)))
sep_rss <- 0
sep_start <- NULL
for (r in sort(unique(dn_df$Run))) {
  fr <- stats::nls(density ~ SSfpl(log(conc), A, B, xmid, scal), data = dn_df[dn_df$Run == r, ])
  sep_rss <- sep_rss + sum(stats::residuals(fr)^2)
  sep_start <- rbind(sep_start, stats::coef(fr))
}
shared_fit <- stats::nls(
  density ~ A[run_i] + (B[run_i] - A[run_i]) / (1 + exp((xmid - log(conc)) / scal[run_i])),
  data = dn_df,
  start = list(A = sep_start[, "A"], B = sep_start[, "B"], scal = sep_start[, "scal"],
               xmid = mean(sep_start[, "xmid"])))
n_runs <- nrow(sep_start)
df_sep <- nrow(dn_df) - 4 * n_runs
f_shared <- ((sum(stats::residuals(shared_fit)^2) - sep_rss) / (n_runs - 1)) / (sep_rss / df_sep)
add_ref("dnase_shared", "statistic", f_shared, "stats::nls(SSfpl), extra sum of squares", rel_tol = 1e-4)
add_ref("dnase_shared", "p_value", stats::pf(f_shared, n_runs - 1, df_sep, lower.tail = FALSE),
        "stats::nls(SSfpl), extra sum of squares", rel_tol = 1e-4)

# 12. aSAH: aneurysmal subarachnoid haemorrhage, Turck et al. (2010), in pROC.
as_df <- write_data(as.data.frame(pROC::aSAH), "asah")
r1 <- pROC::roc(as_df$outcome, as_df$s100b, levels = c("Good", "Poor"), quiet = TRUE)
r2 <- pROC::roc(as_df$outcome, as_df$ndka, levels = c("Good", "Poor"), quiet = TRUE)
c1 <- pROC::ci.auc(r1, method = "delong")
add_ref("asah_roc", "auc_s100b", as.numeric(pROC::auc(r1)), "pROC::auc()")
add_ref("asah_roc", "auc_low_s100b", c1[1], "pROC::ci.auc(delong)")
add_ref("asah_roc", "auc_high_s100b", c1[3], "pROC::ci.auc(delong)")
add_ref("asah_roc", "auc_ndka", as.numeric(pROC::auc(r2)), "pROC::auc()")
rt <- pROC::roc.test(r1, r2, method = "delong")
add_ref("asah_roc", "delong_p", rt$p.value, "pROC::roc.test(delong)")

# 13. PEFR: Bland and Altman (1986), Lancet 327:307-310, Table 1. Wright and
# mini Wright peak flow meters on 17 subjects, two readings each (l/min).
# Transcribed from https://www-users.york.ac.uk/~mb55/meas/ba.htm.
pefr <- data.frame(
  subject  = 1:17,
  wright_1 = c(494, 395, 516, 434, 476, 557, 413, 442, 650, 433, 417, 656, 267, 478, 178, 423, 427),
  wright_2 = c(490, 397, 512, 401, 470, 611, 415, 431, 638, 429, 420, 633, 275, 492, 165, 372, 421),
  mini_1   = c(512, 430, 520, 428, 500, 600, 364, 380, 658, 445, 432, 626, 260, 477, 259, 350, 451),
  mini_2   = c(525, 415, 508, 444, 500, 625, 460, 390, 642, 432, 420, 605, 227, 467, 268, 370, 443))
pefr_df <- write_data(pefr, "pefr")
dd <- pefr_df$wright_1 - pefr_df$mini_1
add_ref("pefr_bland_altman", "bias_published", -2.1, "Bland and Altman (1986), rounded", rel_tol = 0.03)
add_ref("pefr_bland_altman", "sd_published", 38.8, "Bland and Altman (1986), rounded", rel_tol = 0.002)
add_ref("pefr_bland_altman", "bias", mean(dd), "mean of the differences")
add_ref("pefr_bland_altman", "loa_low", mean(dd) - stats::qnorm(0.975) * stats::sd(dd), "mean - 1.96 SD")
add_ref("pefr_bland_altman", "loa_high", mean(dd) + stats::qnorm(0.975) * stats::sd(dd), "mean + 1.96 SD")

# 14. BCG: thirteen trials of the BCG vaccine against tuberculosis, Colditz et
# al. (1994), as distributed in metadat (a metafor dependency). The recipe input
# is the risk ratio and its 95% CI, which is how a paper reports a trial.
bcg <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                       data = metadat::dat.bcg)
z <- stats::qnorm(0.975)
bcg_out <- data.frame(
  trial = paste0(bcg$author, " (", bcg$year, ")"),
  risk_ratio = exp(bcg$yi),
  ci_low = exp(bcg$yi - z * sqrt(bcg$vi)),
  ci_high = exp(bcg$yi + z * sqrt(bcg$vi)),
  ablat = bcg$ablat, alloc = bcg$alloc, stringsAsFactors = FALSE)
bcg_df <- write_data(bcg_out, "bcg")
yi <- log(bcg_df$risk_ratio)
vi <- ((log(bcg_df$ci_high) - log(bcg_df$ci_low)) / (2 * z))^2
dl <- metafor::rma(yi = yi, vi = vi, method = "DL")
fe <- metafor::rma(yi = yi, vi = vi, method = "FE")
add_ref("bcg_meta", "random_estimate", exp(dl$b[1]), "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "random_ci_low", exp(dl$ci.lb), "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "random_ci_high", exp(dl$ci.ub), "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "tau2", dl$tau2, "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "i2", dl$I2, "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "q", dl$QE, "metafor::rma(method = 'DL')")
add_ref("bcg_meta", "fixed_estimate", exp(fe$b[1]), "metafor::rma(method = 'FE')")

# 15. Orthodont: Potthoff and Roy (1964), the distance from the pituitary to the
# pterygomaxillary fissure in 27 children at ages 8, 10, 12 and 14.
or <- as.data.frame(nlme::Orthodont)
or$Subject <- as.character(or$Subject); or$Sex <- as.character(or$Sex)
or_df <- write_data(or[, c("Subject", "Sex", "age", "distance")], "orthodont")

# Mixed ANOVA references for a between factor and a within factor. The F tests
# come from the classical split-plot aov with an Error stratum per subject. The
# Greenhouse-Geisser epsilon and Mauchly's test come from car::Anova on the
# multivariate linear model, which is the route rstatix takes, at full precision.
mixed_refs <- function(case, d, id, grp, time, y) {
  d <- d[, c(id, grp, time, y)]
  names(d) <- c("id", "grp", "time", "y")
  d$id <- factor(d$id); d$grp <- factor(d$grp); d$time_f <- factor(d$time)
  oa <- summary(stats::aov(y ~ grp * time_f + Error(id / time_f), data = d))
  between <- oa[["Error: id"]][[1]]
  within <- oa[["Error: id:time_f"]][[1]]
  rn_b <- trimws(rownames(between)); rn_w <- trimws(rownames(within))
  add_ref(case, "F_group", between[rn_b == "grp", "F value"], "stats::aov(Error())")
  add_ref(case, "p_group", between[rn_b == "grp", "Pr(>F)"], "stats::aov(Error())")
  add_ref(case, "F_time", within[rn_w == "time_f", "F value"], "stats::aov(Error())")
  add_ref(case, "F_interaction", within[rn_w == "grp:time_f", "F value"], "stats::aov(Error())")
  add_ref(case, "p_interaction", within[rn_w == "grp:time_f", "Pr(>F)"], "stats::aov(Error())")
  wide <- stats::reshape(d[, c("id", "grp", "time", "y")], idvar = c("id", "grp"),
                         timevar = "time", direction = "wide")
  ycols <- grep("^y\\.", names(wide), value = TRUE)
  mlm <- stats::lm(as.matrix(wide[, ycols]) ~ grp, data = wide)
  idata <- data.frame(time_f = factor(sub("^y\\.", "", ycols), levels = sub("^y\\.", "", ycols)))
  av <- car::Anova(mlm, idata = idata, idesign = ~time_f, type = 2)
  sv <- summary(av, multivariate = FALSE)
  add_ref(case, "gg_epsilon", sv$pval.adjustments["time_f", "GG eps"], "car::Anova(idesign)")
  add_ref(case, "mauchly_p", sv$sphericity.tests["time_f", "p-value"], "car::Anova(idesign)")
}
mixed_refs("orthodont_timecourse", or_df, "Subject", "Sex", "age", "distance")

# 16. ChickWeight: Crowder and Hand (1990), body weight of chicks on four diets.
# Five chicks died before day 21, so the mixed ANOVA uses the 45 complete chicks.
cw <- as.data.frame(datasets::ChickWeight)
cw$Chick <- as.integer(as.character(cw$Chick)); cw$Diet <- as.integer(as.character(cw$Diet))
cw_df <- write_data(cw[order(cw$Chick, cw$Time), c("Chick", "Diet", "Time", "weight")], "chick_weight")
complete_chicks <- names(which(table(cw_df$Chick) == length(unique(cw_df$Time))))
add_ref("chick_weight_timecourse", "n_complete_chicks", length(complete_chicks),
        "chicks seen at all 12 times")
mixed_refs("chick_weight_timecourse", cw_df[cw_df$Chick %in% complete_chicks, ],
           "Chick", "Diet", "Time", "weight")

# 17. aml: acute myelogenous leukaemia, maintenance chemotherapy or not, Miller
# (1997). status 1 = the event, 0 = censored, which is the swimmer convention of
# 1 = the course ended and 0 = ongoing.
aml <- as.data.frame(survival::aml)
aml$patient <- sprintf("AML-%02d", seq_len(nrow(aml)))
aml_df <- write_data(aml[, c("patient", "x", "time", "status")], "aml")
for (arm in sort(unique(aml_df$x))) {
  key <- tolower(arm)
  add_ref("aml_swimmer", paste0("median_", key), stats::median(aml_df$time[aml_df$x == arm]),
          "stats::median()")
  add_ref("aml_swimmer", paste0("ongoing_", key), sum(aml_df$status[aml_df$x == arm] == 0),
          "count of status 0")
}

# 18. airway: Himes et al. (2014), GSE52778. Airway smooth muscle cells from four
# donors, treated with dexamethasone or not. DESeq2 blocks the donor and tests
# treated against untreated. The table keeps the genes that DESeq2 tested (an
# adjusted p exists) with a mean count of 50 or more, so the file stays small.
data("airway", package = "airway")
air <- airway
air$dex <- stats::relevel(air$dex, "untrt")
dds <- DESeqDataSet(air, design = ~ cell + dex)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("dex", "trt", "untrt")))
res$ensembl <- rownames(res)
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = res$ensembl, column = "SYMBOL",
                             keytype = "ENSEMBL", multiVals = "first")
res$gene <- ifelse(is.na(sym), res$ensembl, sym)
res <- res[!is.na(res$padj) & res$baseMean >= 50, ]
res <- res[order(res$padj, res$ensembl), ]
de_df <- write_data(data.frame(gene = res$gene, log2FoldChange = signif(res$log2FoldChange, 5),
                               pvalue = signif(res$pvalue, 5), padj = signif(res$padj, 5)),
                    "airway_de")
add_ref("airway_volcano", "n_up", sum(de_df$log2FoldChange >= 1 & de_df$padj <= 0.05),
        "DESeq2 results, padj <= 0.05 and log2FC >= 1")
add_ref("airway_volcano", "n_down", sum(de_df$log2FoldChange <= -1 & de_df$padj <= 0.05),
        "DESeq2 results, padj <= 0.05 and log2FC <= -1")
add_ref("airway_volcano", "n_total", nrow(de_df), "DESeq2 results, rows tested")

# The variance-stabilised counts give the heatmap of the 40 strongest
# differentially expressed genes and the PCA of the 500 most variable genes.
vsd <- assay(vst(dds, blind = FALSE))
ann <- data.frame(sample = colnames(vsd), dex = as.character(air$dex), cell = as.character(air$cell))
ann_df <- write_data(ann, "airway_annotation")
top40 <- res$ensembl[abs(res$log2FoldChange) >= 1][1:40]
hm <- data.frame(feature = make.unique(res$gene[match(top40, res$ensembl)]),
                 round(vsd[top40, ], 4), check.names = FALSE)
hm_df <- write_data(hm, "airway_heatmap")
for (i in 1:5) {
  v <- as.numeric(hm_df[i, -1])
  g <- ann_df$dex[match(names(hm_df)[-1], ann_df$sample)]
  add_ref("airway_heatmap", paste0("p_", hm_df$feature[i]),
          stats::t.test(v[g == "untrt"], v[g == "trt"])$p.value, "stats::t.test(var.equal = FALSE)")
}
var_rank <- order(-apply(vsd, 1, stats::var), rownames(vsd))[1:500]
pcx <- t(round(vsd[var_rank, ], 4))
colnames(pcx) <- make.unique(ifelse(is.na(sym[rownames(vsd)[var_rank]]), rownames(vsd)[var_rank],
                                    sym[rownames(vsd)[var_rank]]))
pca_df <- write_data(data.frame(sample = rownames(pcx), dex = ann$dex, cell = ann$cell, pcx,
                                check.names = FALSE), "airway_pca")
pca_air <- stats::prcomp(pca_df[, -(1:3)], center = TRUE, scale. = TRUE)
ve_air <- pca_air$sdev^2 / sum(pca_air$sdev^2) * 100
add_ref("airway_pca", "pc1_variance_pct", ve_air[1], "stats::prcomp(scale. = TRUE)")
add_ref("airway_pca", "pc2_variance_pct", ve_air[2], "stats::prcomp(scale. = TRUE)")

# GO biological-process over-representation of the differentially expressed genes
# (padj < 0.05 and |log2FC| >= 1) against the tested genes, by the hypergeometric
# test with Benjamini-Hochberg adjustment, in the layout of clusterProfiler.
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = res$ensembl,
                                             columns = c("GOALL", "ONTOLOGYALL"), keytype = "ENSEMBL"))
go <- unique(go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", c("ENSEMBL", "GOALL")])
universe <- unique(go$ENSEMBL)
de_genes <- intersect(res$ensembl[res$padj < 0.05 & abs(res$log2FoldChange) >= 1], universe)
term_size <- table(go$GOALL)
term_hits <- table(go$GOALL[go$ENSEMBL %in% de_genes])
terms <- names(term_size)[term_size >= 15 & term_size <= 500]
hits <- as.integer(term_hits[terms]); hits[is.na(hits)] <- 0L
size <- as.integer(term_size[terms])
pv <- stats::phyper(hits - 1, size, length(universe) - size, length(de_genes), lower.tail = FALSE)
ora <- data.frame(ID = terms, Count = hits, size = size, pvalue = pv,
                  p.adjust = stats::p.adjust(pv, "BH"), stringsAsFactors = FALSE)
ora <- ora[ora$Count > 0, ]
ora <- ora[order(ora$p.adjust, ora$pvalue, ora$ID), ]
ora$Description <- AnnotationDbi::Term(GO.db::GOTERM[ora$ID])
ora$GeneRatio <- sprintf("%d/%d", ora$Count, length(de_genes))
ora$BgRatio <- sprintf("%d/%d", ora$size, length(universe))
enr_df <- write_data(data.frame(ID = ora$ID, Description = ora$Description, GeneRatio = ora$GeneRatio,
                                BgRatio = ora$BgRatio, pvalue = signif(ora$pvalue, 5),
                                p.adjust = signif(ora$p.adjust, 5), Count = ora$Count)[1:30, ],
                     "airway_enrichment")
add_ref("airway_enrichment", "top_count", enr_df$Count[1], "stats::phyper(), top term")
add_ref("airway_enrichment", "top_p_adjust", enr_df$p.adjust[1], "stats::p.adjust(BH), top term")
add_ref("airway_enrichment", "top_gene_ratio", enr_df$Count[1] / length(de_genes), "Count / DE genes")

# Gene-set membership of the DE genes in five enriched GO terms, chosen down the
# ranking and skipping a term whose genes overlap an already chosen term by a
# Jaccard index of 0.5 or more, so that the UpSet plot shows distinct sets.
members <- split(go$ENSEMBL[go$ENSEMBL %in% de_genes], go$GOALL[go$ENSEMBL %in% de_genes])
chosen <- character(0)
for (id in ora$ID) {
  overlap <- vapply(chosen, jaccard_with, numeric(1), genes = members[[id]], members = members)
  if (all(overlap < 0.5)) chosen <- c(chosen, id)
  if (length(chosen) == 5) break
}
in_any <- sort(unique(unlist(members[chosen])))
up <- data.frame(gene = make.unique(ifelse(is.na(sym[in_any]), in_any, sym[in_any])),
                 check.names = FALSE)
for (id in chosen) up[[AnnotationDbi::Term(GO.db::GOTERM[[id]])]] <- as.integer(in_any %in% members[[id]])
up_df <- write_data(up, "airway_upset")
codes <- apply(as.matrix(up_df[, -1]), 1, paste, collapse = "")
add_ref("airway_upset", "n_elements", nrow(up_df), "genes in at least one set")
add_ref("airway_upset", "largest_intersection", max(table(codes)), "base::table() of membership patterns")
add_ref("airway_upset", "n_intersections", length(unique(codes)), "distinct membership patterns")

# 19. anorexia: Hand et al. (1993), weight of young women before and after
# family therapy (the FT group of MASS::anorexia), in long form with a subject id.
an <- MASS::anorexia[MASS::anorexia$Treat == "FT", ]
an_long <- data.frame(subject = rep(sprintf("FT%02d", seq_len(nrow(an))), 2),
                      condition = rep(c("pre", "post"), each = nrow(an)),
                      weight_lbs = c(an$Prewt, an$Postwt))
an_df <- write_data(an_long, "anorexia_ft")
an_t <- stats::t.test(an_df$weight_lbs[an_df$condition == "post"],
                      an_df$weight_lbs[an_df$condition == "pre"], paired = TRUE)
add_ref("anorexia_paired", "statistic", an_t$statistic, "stats::t.test(paired = TRUE)")
add_ref("anorexia_paired", "p_value", an_t$p.value, "stats::t.test(paired = TRUE)")

# 20. veteran: the Veterans' Administration lung cancer trial, Kalbfleisch and
# Prentice (1980). status is already 1 = dead, 0 = censored; trt is relabelled.
vet <- survival::veteran
vet$trt <- ifelse(vet$trt == 1, "standard", "test")
vet_df <- write_data(vet, "veteran")
vsd_lr <- survival::survdiff(Surv(time, status) ~ celltype, data = vet_df)
add_ref("veteran_km", "log_rank_chisq", vsd_lr$chisq, "survival::survdiff()")
add_ref("veteran_km", "log_rank_p",
        stats::pchisq(vsd_lr$chisq, length(vsd_lr$n) - 1, lower.tail = FALSE), "survival::survdiff()")
vcx <- summary(survival::coxph(Surv(time, status) ~ trt + celltype + karno + age, data = vet_df))$coefficients
for (term in rownames(vcx)) {
  add_ref("veteran_cox", paste0("HR_", term), vcx[term, "exp(coef)"], "survival::coxph()")
  add_ref("veteran_cox", paste0("p_", term), vcx[term, "Pr(>|z|)"], "survival::coxph()")
}

# 21. state.x77: facts about the 50 US states, U.S. Department of Commerce (1977).
# Two column names hold a space ("Life Exp", "HS Grad"), as R distributes them.
st_df <- write_data(data.frame(state = rownames(datasets::state.x77), datasets::state.x77,
                               check.names = FALSE), "state_x77")
st_num <- st_df[, -1]
st_p <- NULL
for (i in 1:(ncol(st_num) - 1)) for (j in (i + 1):ncol(st_num)) {
  st_p <- rbind(st_p, data.frame(a = names(st_num)[i], b = names(st_num)[j],
                                 p = stats::cor.test(st_num[[i]], st_num[[j]])$p.value))
}
st_p$padj <- stats::p.adjust(st_p$p, "BH")
add_ref("state_correlation_heatmap", "r_murder_life_exp", stats::cor(st_num$Murder, st_num$`Life Exp`),
        "stats::cor()")
add_ref("state_correlation_heatmap", "padj_murder_life_exp",
        st_p$padj[st_p$a == "Life Exp" & st_p$b == "Murder"], "stats::cor.test(), p.adjust(BH)")
add_ref("state_correlation_heatmap", "r_illiteracy_frost", stats::cor(st_num$Illiteracy, st_num$Frost),
        "stats::cor()")

# 22. faithful: the eruption durations of the Old Faithful geyser, Azzalini and
# Bowman (1990). The durations are bimodal, so the normality test rejects.
fa_df <- write_data(datasets::faithful, "faithful")
er <- fa_df$eruptions
er_sw <- stats::shapiro.test(er)
add_ref("faithful_distribution", "shapiro_w", er_sw$statistic, "stats::shapiro.test()")
add_ref("faithful_distribution", "shapiro_p", er_sw$p.value, "stats::shapiro.test()")
add_ref("faithful_distribution", "mean", mean(er), "base::mean()")
add_ref("faithful_distribution", "sd", stats::sd(er), "stats::sd()")
add_ref("faithful_distribution", "bandwidth", stats::bw.nrd0(er), "stats::bw.nrd0()")
add_ref("faithful_distribution", "bin_width", 2 * stats::IQR(er) * length(er)^(-1 / 3),
        "2 * stats::IQR() * n^(-1/3)")

# 23. ToothGrowth (case 4): the tooth length within each supplement type.
for (s in c("OJ", "VC")) {
  tg_v <- tg$len[tg$supp == s]
  tg_sw <- stats::shapiro.test(tg_v)
  add_ref("tooth_growth_distribution", paste0("shapiro_w_", s), tg_sw$statistic, "stats::shapiro.test()")
  add_ref("tooth_growth_distribution", paste0("shapiro_p_", s), tg_sw$p.value, "stats::shapiro.test()")
  add_ref("tooth_growth_distribution", paste0("mean_", s), mean(tg_v), "base::mean()")
  add_ref("tooth_growth_distribution", paste0("sd_", s), stats::sd(tg_v), "stats::sd()")
  add_ref("tooth_growth_distribution", paste0("bandwidth_", s), stats::bw.nrd0(tg_v), "stats::bw.nrd0()")
}

# 24. housing: satisfaction with housing in Copenhagen, Madsen (1976), in MASS.
# A sample is one combination of influence, housing type and contact; the
# categories are the satisfaction levels and the group is the contact level.
hs <- MASS::housing
hs_df <- write_data(data.frame(sample = paste(hs$Infl, hs$Type, hs$Cont, sep = "_"),
                               Infl = as.character(hs$Infl), Type = as.character(hs$Type),
                               Cont = as.character(hs$Cont), Sat = as.character(hs$Sat),
                               Freq = hs$Freq, stringsAsFactors = FALSE), "housing")
hs_tot <- tapply(hs_df$Freq, hs_df$sample, sum)
hs_df$prop <- hs_df$Freq / as.numeric(hs_tot[hs_df$sample])
hs_p <- NULL
for (s in c("Low", "Medium", "High")) {
  hs_s <- hs_df[hs_df$Sat == s, ]
  hs_wt <- suppressWarnings(stats::wilcox.test(hs_s$prop[hs_s$Cont == "High"],
                                               hs_s$prop[hs_s$Cont == "Low"]))
  add_ref("housing_composition", paste0("W_", s), hs_wt$statistic, "stats::wilcox.test()")
  add_ref("housing_composition", paste0("p_", s), hs_wt$p.value, "stats::wilcox.test()")
  hs_p <- c(hs_p, hs_wt$p.value)
}
hs_adj <- stats::p.adjust(hs_p, method = "BH")
for (k in seq_along(hs_adj)) {
  add_ref("housing_composition", paste0("p_adj_", c("Low", "Medium", "High")[k]), hs_adj[k],
          "stats::p.adjust(method = \"BH\")")
}

# 22. chickwts: chick weight after six weeks on six feed supplements, Anonymous
# (1948), in datasets. Each feed is compared with casein, the estimation_plot
# reference: the mean difference and the Welch two-sample t-test.
cw_feed <- write_data(datasets::chickwts[, c("feed", "weight")], "chickwts")
casein <- cw_feed$weight[cw_feed$feed == "casein"]
for (fd in setdiff(sort(unique(cw_feed$feed)), "casein")) {
  v <- cw_feed$weight[cw_feed$feed == fd]
  wt <- stats::t.test(v, casein, var.equal = FALSE)
  add_ref("chickwts_estimation", paste0("mean_diff_", fd), mean(v) - mean(casein), "mean difference")
  add_ref("chickwts_estimation", paste0("statistic_", fd), wt$statistic, "stats::t.test(var.equal = FALSE)")
  add_ref("chickwts_estimation", paste0("p_", fd), wt$p.value, "stats::t.test(var.equal = FALSE)")
}

# 23. Puromycin: the rate of an enzymatic reaction in cells treated or not with
# puromycin, Treloar (1974), in datasets. Michaelis-Menten fits per state, and
# the nested shared-K and shared-Vm models against the separate fits by anova().
# The fits use a tight convergence tolerance, so the reference is the optimum.
pu <- write_data(datasets::Puromycin, "puromycin")
pu_ctrl <- stats::nls.control(maxiter = 200, tol = 1e-8, minFactor = 1e-10)
for (st in c("treated", "untreated")) {
  fp <- stats::nls(rate ~ Vm * conc / (K + conc), data = pu[pu$state == st, ],
                   start = list(Vm = 200, K = 0.1), control = pu_ctrl)
  cs <- summary(fp)$coefficients
  add_ref("puromycin_kinetics", paste0("vmax_", st), cs["Vm", "Estimate"], "stats::nls()")
  add_ref("puromycin_kinetics", paste0("vmax_se_", st), cs["Vm", "Std. Error"], "stats::nls()")
  add_ref("puromycin_kinetics", paste0("km_", st), cs["K", "Estimate"], "stats::nls()")
  add_ref("puromycin_kinetics", paste0("km_se_", st), cs["K", "Std. Error"], "stats::nls()")
}
pu$si <- match(pu$state, c("treated", "untreated"))
pu_sep <- stats::nls(rate ~ Vm[si] * conc / (K[si] + conc), data = pu,
                     start = list(Vm = c(200, 160), K = c(0.1, 0.1)), control = pu_ctrl)
pu_shk <- stats::nls(rate ~ Vm[si] * conc / (K + conc), data = pu,
                     start = list(Vm = c(200, 160), K = 0.1), control = pu_ctrl)
pu_shv <- stats::nls(rate ~ Vm * conc / (K[si] + conc), data = pu,
                     start = list(Vm = 200, K = c(0.1, 0.1)), control = pu_ctrl)
ak <- stats::anova(pu_shk, pu_sep)
av <- stats::anova(pu_shv, pu_sep)
add_ref("puromycin_kinetics", "shared_km_F", ak$`F value`[2], "stats::anova(nls)")
add_ref("puromycin_kinetics", "shared_km_p", ak$`Pr(>F)`[2], "stats::anova(nls)")
add_ref("puromycin_kinetics", "shared_vmax_F", av$`F value`[2], "stats::anova(nls)")
add_ref("puromycin_kinetics", "shared_vmax_p", av$`Pr(>F)`[2], "stats::anova(nls)")

# 24. DNase as a standard curve: run 1 (conc, density) is the set of standards,
# and the densities of run 2, with the concentration removed, are the unknowns.
# The reference curve is the self-starting SSfpl fit to the standards, and each
# unknown is the inverse of that curve when its density lies within the curve's
# range over the standards. Two optimisers stop at their own tolerance, so the
# concentrations carry the 1e-4 of dnase_4pl.
r1 <- dn[dn$Run == 1, ]; r2 <- dn[dn$Run == 2, ]
sc_df <- write_data(rbind(
  data.frame(sample = sprintf("run1_%02d", seq_len(nrow(r1))), conc = r1$conc, density = r1$density),
  data.frame(sample = sprintf("run2_%02d", seq_len(nrow(r2))), conc = NA, density = r2$density)),
  "dnase_standard_curve")
sc_std <- sc_df[!is.na(sc_df$conc), ]
sc_unk <- sc_df[is.na(sc_df$conc), ]
sc_fit <- stats::coef(stats::nls(density ~ SSfpl(log(conc), A, B, xmid, scal), data = sc_std))
sc_f <- sc_fit[["A"]] + (sc_fit[["B"]] - sc_fit[["A"]]) /
  (1 + exp((sc_fit[["xmid"]] - log(range(sc_std$conc))) / sc_fit[["scal"]]))
sc_in <- sc_unk$density >= min(sc_f) & sc_unk$density <= max(sc_f)
sc_conc <- exp(sc_fit[["xmid"]] - sc_fit[["scal"]] *
                 log((sc_fit[["B"]] - sc_fit[["A"]]) / (sc_unk$density - sc_fit[["A"]]) - 1))
add_ref("dnase_standard_curve", "n_in_range", sum(sc_in), "inverse of stats::nls(SSfpl)")
add_ref("dnase_standard_curve", "n_below", sum(!sc_in & sc_unk$density < min(sc_f)),
        "inverse of stats::nls(SSfpl)")
add_ref("dnase_standard_curve", "n_above", sum(!sc_in & sc_unk$density > max(sc_f)),
        "inverse of stats::nls(SSfpl)")
for (i in which(sc_in)[c(1, ceiling(sum(sc_in) / 2), sum(sc_in))]) {
  add_ref("dnase_standard_curve", paste0("conc_", sc_unk$sample[i]), sc_conc[i],
          "inverse of stats::nls(SSfpl)", rel_tol = 1e-4)
}

# 22. mgus2: monoclonal gammopathy of undetermined significance, Kyle et al.
# (2002). Competing risks of progression to a plasma cell malignancy (event 1)
# and death without progression (event 2), as in the survival vignette "compete".
mg <- survival::mgus2
mg$etime <- ifelse(mg$pstat == 0, mg$futime, mg$ptime)
mg$event <- ifelse(mg$pstat == 0, 2 * mg$death, 1)
mg$sex <- as.character(mg$sex)
mg <- mg[stats::complete.cases(mg[, c("etime", "event", "sex")]), ]
mg_df <- write_data(mg[, c("id", "age", "sex", "etime", "event")], "mgus2")
mg_fit <- cmprsk::cuminc(mg_df$etime, mg_df$event, mg_df$sex)
for (cause in c("1", "2")) {
  add_ref("mgus2_cif", paste0("gray_stat_", cause), mg_fit$Tests[cause, "stat"], "cmprsk::cuminc()$Tests")
  add_ref("mgus2_cif", paste0("gray_p_", cause), mg_fit$Tests[cause, "pv"], "cmprsk::cuminc()$Tests")
}
mg_tp <- cmprsk::timepoints(mg_fit, c(120, 240))$est
for (curve in rownames(mg_tp)) {
  for (tt in colnames(mg_tp)) {
    add_ref("mgus2_cif", sprintf("cif_%s_%s", gsub(" ", "_", curve), tt), mg_tp[curve, tt],
            "cmprsk::timepoints()")
  }
}

# 23. Pima.tr and Pima.te: diabetes in Pima Indian women, Smith et al. (1988).
# A logistic model fitted on Pima.tr predicts the risk of each woman in Pima.te,
# which is then the validation set of the calibration recipe.
pima_fit <- stats::glm(type ~ npreg + glu + bmi + ped + age, family = stats::binomial(),
                       data = MASS::Pima.tr)
pima_te <- MASS::Pima.te
pima_te$type <- as.character(pima_te$type)
pima_te$predicted <- stats::predict(pima_fit, newdata = MASS::Pima.te, type = "response")
pima_df <- write_data(pima_te, "pima_te")
pima_y <- as.integer(pima_df$type == "Yes")
pima_p <- pima_df$predicted
pima_hl <- ResourceSelection::hoslem.test(pima_y, pima_p, g = 10)
add_ref("pima_calibration", "hl_statistic", pima_hl$statistic, "ResourceSelection::hoslem.test(g = 10)")
add_ref("pima_calibration", "hl_df", pima_hl$parameter, "ResourceSelection::hoslem.test(g = 10)")
add_ref("pima_calibration", "hl_p", pima_hl$p.value, "ResourceSelection::hoslem.test(g = 10)")
add_ref("pima_calibration", "brier", mean((pima_p - pima_y)^2), "mean((p - y)^2)")
pima_lp <- stats::qlogis(pima_p)
add_ref("pima_calibration", "slope",
        stats::coef(stats::glm(pima_y ~ pima_lp, family = stats::binomial()))[["pima_lp"]],
        "stats::glm(y ~ logit(p))")
add_ref("pima_calibration", "intercept",
        stats::coef(stats::glm(pima_y ~ 1 + offset(pima_lp), family = stats::binomial()))[[1]],
        "stats::glm(y ~ 1 + offset(logit(p)))")

# 24. pbcseq: the Mayo Clinic primary biliary cholangitis trial with repeated
# visits, Murtaugh et al. (1994). Edema at the first visit and at the visit
# closest to day 730 among the visits between day 365 and day 1095.
pbs <- survival::pbcseq
pbs <- pbs[order(pbs$id, pbs$day), ]
pb_first <- pbs[!duplicated(pbs$id), c("id", "day", "edema")]
pb_win <- pbs[pbs$day >= 365 & pbs$day <= 1095, c("id", "day", "edema")]
pb_win <- pb_win[order(pb_win$id, abs(pb_win$day - 730), pb_win$day), ]
pb_second <- pb_win[!duplicated(pb_win$id), ]
pb_long <- rbind(data.frame(id = pb_first$id, year = 0, day = pb_first$day, edema = pb_first$edema),
                 data.frame(id = pb_second$id, year = 2, day = pb_second$day, edema = pb_second$edema))
pb_df <- write_data(pb_long[order(pb_long$id, pb_long$year), ], "pbcseq_edema")
pb_pair <- merge(pb_df[pb_df$year == 0, ], pb_df[pb_df$year == 2, ], by = "id")
pb_lv <- sort(unique(c(pb_pair$edema.x, pb_pair$edema.y)))
pb_tab <- table(factor(pb_pair$edema.x, levels = pb_lv), factor(pb_pair$edema.y, levels = pb_lv))
pb_mt <- stats::mcnemar.test(pb_tab, correct = FALSE)
add_ref("pbcseq_sankey", "statistic", pb_mt$statistic, "stats::mcnemar.test(correct = FALSE)")
add_ref("pbcseq_sankey", "df", pb_mt$parameter, "stats::mcnemar.test(correct = FALSE)")
add_ref("pbcseq_sankey", "p_value", pb_mt$p.value, "stats::mcnemar.test(correct = FALSE)")
add_ref("pbcseq_sankey", "n_paired", nrow(pb_pair), "patients with both visits")
add_ref("pbcseq_sankey", "n_0_to_0.5", pb_tab["0", "0.5"], "base::table()")
add_ref("pbcseq_sankey", "n_0.5_to_0", pb_tab["0.5", "0"], "base::table()")

ref_df <- do.call(rbind, refs)
utils::write.csv(ref_df, file.path("benchmarks", "reference_values.csv"), row.names = FALSE)
cat(sprintf("Wrote %d datasets and %d reference values.\n",
            length(list.files(out_dir, pattern = "\\.csv$")), nrow(ref_df)))
