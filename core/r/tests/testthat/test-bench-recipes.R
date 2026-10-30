# Behaviour of the bench-assay recipes (estimation_plot, enzyme_kinetics and
# standard_curve) on small fixed inputs. The numbers on public data are checked
# in test-benchmarks.R.

bench_recipe_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

fixed_unpaired <- function() {
  data.frame(group = rep(c("ctrl", "trt"), c(7, 6)),
             value = c(10.2, 11.5, 9.8, 12.1, 10.9, 11.3, 9.4, 13.4, 12.8, 14.1, 12.2, 13.9, 11.7))
}

fixed_paired <- function() {
  data.frame(id = rep(paste0("s", 1:5), 2), cond = rep(c("pre", "post"), each = 5),
             value = c(5.1, 4.8, 6.2, 5.5, 5.9, 5.9, 5.1, 7.0, 6.1, 6.8))
}

fpl_signal <- function(conc, bottom, top, ec50, hill) {
  bottom + (top - bottom) / (1 + 10^((log10(ec50) - log10(conc)) * hill))
}

test_that("the Park-Miller generator reaches the published state after 10000 draws", {
  # Park, Miller and Stockmeyer (1993): seed 1, multiplier 48271, state 399268537.
  g <- .est_uniforms(10000, 1, .est_block_multipliers(.est_block))
  expect_identical(g$state, 399268537)
  expect_true(all(g$u > 0 & g$u < 1))
})

test_that("estimation_plot gives the BCa interval of the Python engine (unpaired)", {
  res <- recipe_estimation_plot(fixed_unpaired(), bench_recipe_spec("estimation_plot", x = "group", y = "value"))
  # The Python engine computes these two values from the same Park-Miller resamples.
  expect_equal(res$stats$ci_low, 1.3148939612822759, tolerance = 1e-9)
  expect_equal(res$stats$ci_high, 3.2523809523809533, tolerance = 1e-9)
  tt <- stats::t.test(fixed_unpaired()$value[7 + 1:6], fixed_unpaired()$value[1:7])
  expect_equal(res$stats$statistic, unname(tt$statistic), tolerance = 1e-12)
  expect_equal(res$stats$p_value, tt$p.value, tolerance = 1e-12)
  expect_equal(res$stats$mean_diff, mean(fixed_unpaired()$value[8:13]) - mean(fixed_unpaired()$value[1:7]))
  expect_equal(res$test_meta$reference, "ctrl")
})

test_that("estimation_plot gives the BCa interval of the Python engine (paired)", {
  res <- recipe_estimation_plot(fixed_paired(),
                                bench_recipe_spec("estimation_plot", x = "cond", y = "value", id = "id",
                                                  paired = TRUE, ref = "pre"))
  expect_equal(res$stats$ci_low, 0.4, tolerance = 1e-9)
  expect_equal(res$stats$ci_high, 0.8199999999999996, tolerance = 1e-9)
  expect_equal(res$stats$test, "paired t-test")
  expect_equal(res$stats$n_ref, 5)
})

test_that("estimation_plot puts --ref first and compares the other levels in sorted order", {
  df <- data.frame(g = rep(c("b", "c", "a"), each = 4), v = c(1:4, 5:8, 3:6))
  res <- recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x = "g", y = "v", ref = "c"))
  expect_equal(res$stats$reference, c("c", "c"))
  expect_equal(res$stats$group, c("a", "b"))
  expect_error(recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x = "g", y = "v", ref = "z")),
               "is not a level")
})

test_that("estimation_plot --paired leaves out a subject without both values", {
  df <- fixed_paired()[-10, ]
  res <- recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x = "cond", y = "value", id = "id",
                                                      paired = TRUE, ref = "pre"))
  expect_equal(res$stats$n_group, 4)
  expect_true(any(grepl("left out 1 subject", res$clean_steps)))
})

test_that("enzyme_kinetics recovers Vmax and Km and their Prism-style intervals", {
  s <- rep(c(0.5, 1, 2, 4, 8, 16, 32), each = 2)
  wobble <- rep(c(1.01, 0.99), 7)
  df <- data.frame(s = s, v = 10 * s / (2 + s) * wobble)
  res <- recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x = "s", y = "v"))
  r <- res$stats
  expect_equal(r$vmax, 10, tolerance = 0.01)
  expect_equal(r$km, 2, tolerance = 0.02)
  tq <- stats::qt(0.975, nrow(df) - 2)
  expect_equal(r$km_ci_high - r$km, tq * r$km_se, tolerance = 1e-10)
  expect_null(res$test_meta$shared_km_test)
})

test_that("enzyme_kinetics tests a shared Km and a shared Vmax with groups", {
  s <- rep(c(0.5, 1, 2, 4, 8, 16, 32), 2)
  df <- data.frame(s = s, grp = rep(c("none", "inhibitor"), each = 7),
                   v = c(10 * s[1:7] / (2 + s[1:7]), 10 * s[8:14] / (6 + s[8:14])) * rep(c(1.01, 0.99, 1), length.out = 14))
  res <- recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x = "s", y = "v", group = "grp"))
  expect_equal(res$stats$group, c("none", "inhibitor"))
  expect_lt(res$test_meta$shared_km_test$p_value, 0.001)
  expect_gt(res$test_meta$shared_vmax_test$p_value, 0.05)
  expect_equal(res$test_meta$shared_km_test$df1, 1)
})

test_that("enzyme_kinetics needs three substrate concentrations per curve", {
  df <- data.frame(s = c(1, 2, 1, 2), v = c(1, 2, 1.1, 2.1))
  expect_error(recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x = "s", y = "v")),
               "at least 3 distinct positive substrate")
})

test_that("standard_curve picks the 4PL from five levels and flags the out-of-range unknowns", {
  conc <- rep(c(0, 1, 3, 10, 30, 100, 300), each = 2)
  sig <- c(0.05, 0.05, fpl_signal(conc[-(1:2)], 0.05, 2.5, 25, 1.1)) * rep(c(1.01, 0.99), 7)
  unk <- c(0.02, fpl_signal(c(2, 20, 200), 0.05, 2.5, 25, 1.1), 2.6)
  df <- data.frame(id = c(paste0("std", seq_along(conc)), paste0("u", 1:5)),
                   conc = c(conc, rep(NA, 5)), od = c(sig, unk))
  res <- recipe_standard_curve(df, bench_recipe_spec("standard_curve", x = "conc", y = "od", id = "id"))
  expect_equal(res$test_meta$model, "four-parameter logistic")
  expect_equal(res$stats$flag, c("below range", "in range", "in range", "in range", "above range"))
  expect_true(all(is.na(res$stats$concentration[c(1, 5)])))
  expect_equal(res$stats$concentration[2:4], c(2, 20, 200), tolerance = 0.03)
  rec <- res$test_meta$recovery
  expect_equal(rec$n_levels, 6)
  # Near the plateaus a 1% error in signal is a large error in concentration, so
  # only the levels on the steep part of the curve are held to a tight recovery.
  mid <- vapply(rec$levels[3:4], `[[`, numeric(1), "recovery_pct")
  expect_true(all(abs(mid - 100) < 3))
  expect_true(any(grepl("concentration of zero out of the four-parameter", res$clean_steps)))
})

test_that("standard_curve falls back to a line with four levels and names unknowns by row", {
  df <- data.frame(conc = c(0, 0, 10, 10, 20, 20, 40, 40, NA, NA),
                   signal = c(0.02, 0.03, 0.21, 0.23, 0.41, 0.44, 0.8, 0.83, 0.33, 0.95))
  res <- recipe_standard_curve(df, bench_recipe_spec("standard_curve", x = "conc", y = "signal"))
  expect_equal(res$test_meta$model, "linear")
  expect_equal(res$stats$id, c("row 9", "row 10"))
  fit <- stats::lm(signal ~ conc, data = df[1:8, ])
  expect_equal(res$stats$concentration[1], (0.33 - stats::coef(fit)[[1]]) / stats::coef(fit)[[2]],
               tolerance = 1e-12)
  expect_equal(res$stats$flag[2], "above range")
  expect_error(recipe_standard_curve(df, bench_recipe_spec("standard_curve", x = "conc", y = "signal",
                                                           test = "4pl")),
               "needs at least 4 distinct positive")
})

test_that("standard_curve reads a falling curve the right way round", {
  conc <- c(1, 3, 10, 30, 100, 300)
  df <- data.frame(conc = c(conc, NA, NA),
                   od = c(fpl_signal(conc, 0.1, 2, 20, -1) * c(1.01, 0.99, 1.02, 0.98, 1.01, 0.99), 2.05, 0.05))
  res <- recipe_standard_curve(df, bench_recipe_spec("standard_curve", x = "conc", y = "od"))
  # A competitive assay: a high signal means a low concentration.
  expect_equal(res$stats$flag, c("below range", "above range"))
})

test_that("the category labels wrap after spaces and after underscores, as in Python", {
  expect_equal(.est_wrap("Heterozygous knockout vehicle"), "Heterozygous\nknockout\nvehicle")
  expect_equal(.est_wrap("stimulated_plus_inhibitor"), "stimulated_\nplus_\ninhibitor")
})

test_that("estimation_plot says that its P values carry no adjustment", {
  expect_match(.est_p_text("Welch two-sample t-test", 5, "casein"),
               "no adjustment for the 5 comparisons against casein", fixed = TRUE)
  expect_equal(.est_p_text("paired t-test", 1, "pre"), "A paired t-test gives the P value. ")
})

test_that("the standard curve recovery counts only the levels read back", {
  rec <- list(n_levels = 8, n_read_back = 6, min_pct = 94.9, max_pct = 111.9, n_within_80_120 = 6)
  txt <- .stc_recovery_text(rec)
  expect_match(txt, "6 of the 8 standard levels", fixed = TRUE)
  expect_match(txt, "6 of these 6 levels fell within 80% to 120%", fixed = TRUE)
  expect_match(txt, "The other 2 levels could not be read back", fixed = TRUE)
  rec$n_read_back <- 8; rec$n_within_80_120 <- 8
  expect_false(grepl("could not be read back", .stc_recovery_text(rec)))
})

# A Michaelis-Menten fit whose Km lies above the highest substrate, or whose CI
# holds 0, is reported as not determined and stays out of the F tests (F50). The
# Python suite uses the same data.
ek_two_curves <- function() {
  s <- c(1, 2, 5, 10, 20, 50)
  noise <- c(0.8, -1.1, 0.6, -0.4, 0.9, -0.7)
  data.frame(s = rep(s, 2), group = rep(c("mm", "linear"), each = 6),
             v = c(100 * s / (5 + s) + noise, 0.1 * s + noise / 10))
}

test_that("enzyme_kinetics reports a fit above the substrate range as not determined", {
  res <- recipe_enzyme_kinetics(ek_two_curves(),
                                bench_recipe_spec("enzyme_kinetics", x = "s", y = "v", group = "group"))
  st <- res$stats
  expect_equal(st$fit_status, c("fitted", "not determined"))
  expect_equal(st$status_reason[2], "Km lies above the highest substrate concentration")
  expect_true(is.na(st$km[2]))
  expect_equal(st$km[1], 5, tolerance = 0.1)
  expect_null(res$test_meta$shared_km_test)
  expect_match(res$methods, "fewer than two curves were determined", fixed = TRUE)
  expect_match(res$plot$labels$subtitle, "Not determined, left out of the tests: linear", fixed = TRUE)
})
