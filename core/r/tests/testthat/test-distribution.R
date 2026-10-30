# Behaviour of the distribution recipe on small synthetic inputs. The numbers
# themselves are checked against public benchmarks in test-benchmarks.R.

dist_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "distribution", data = "d.csv", test = "auto"),
                                  list(...)))
}

test_that("distribution uses the Freedman-Diaconis width anchored at a multiple of it", {
  set.seed(3)
  v <- stats::rnorm(120, 50, 10)
  res <- recipe_distribution(data.frame(value = v), dist_spec(y = "value"))
  w <- 2 * stats::IQR(v) * length(v)^(-1 / 3)
  expect_equal(res$stats$bin_width, w, tolerance = 1e-12)
  expect_equal(res$test_meta$bins$start, floor(min(v) / w) * w)
  expect_equal(res$test_meta$bins$rule, "Freedman-Diaconis")
  expect_equal(res$stats$bandwidth, stats::bw.nrd0(v), tolerance = 1e-12)
  expect_equal(res$stats$group, "all")
  expect_equal(res$stats$shapiro_p, stats::shapiro.test(v)$p.value, tolerance = 1e-12)
})

test_that("distribution counts every value once in left-closed bins", {
  v <- c(0, 1, 1, 2, 3, 3, 3, 4, 5, 8)
  bins <- .dist_bins(v)
  counts <- .dist_counts(v, bins)
  expect_equal(sum(counts), length(v))
  expect_equal(counts[1], sum(v >= bins$start & v < bins$start + bins$width))
})

test_that("distribution falls back to Sturges when the IQR is zero", {
  v <- c(rep(5, 20), 1, 9)
  res <- recipe_distribution(data.frame(value = v), dist_spec(y = "value"))
  expect_equal(res$test_meta$bins$rule, "Sturges")
  expect_equal(res$stats$bin_width, 8 / ceiling(log2(22) + 1))
})

test_that("distribution reports each group and leaves Shapiro-Wilk empty below 3 values", {
  set.seed(4)
  df <- data.frame(g = c(rep("a", 30), rep("b", 2)), v = c(stats::rnorm(30), 0.1, 0.4))
  res <- recipe_distribution(df, dist_spec(x = "g", y = "v"))
  expect_equal(res$stats$group, c("a", "b"))
  expect_equal(res$stats$n, c(30, 2))
  expect_true(is.na(res$stats$shapiro_p[2]))
  expect_equal(length(unique(res$stats$bin_width)), 1)
  expect_named(res$qc$panels, c("a", "b"))
  expect_s3_class(ggplot2::ggplot_build(res$plot), "ggplot_built")
})

# A variance that is not finite at double precision stops the run with a clear
# message, in distribution and correlation_heatmap (F83).
test_that("values near 1e300 or 1e-300 stop distribution and correlation_heatmap", {
  big <- c(-7.858e299, 4.285e299, 9.2e299, 1.1e299, -3.3e299)
  expect_error(recipe_distribution(data.frame(value = big), dist_spec(y = "value")),
               "the variance of column 'value' is not finite at double precision")
  tiny <- c(1.2e-300, 3.4e-300, 2.2e-300, 5.1e-300, 0.7e-300)
  expect_error(recipe_distribution(data.frame(value = tiny), dist_spec(y = "value")),
               "not finite at double precision")
  df <- data.frame(big = big, b = c(1, 3, 2, 5, 4), c = c(2, 1, 4, 3, 6))
  expect_error(recipe_correlation_heatmap(df, spec_from_opt(list(recipe = "correlation_heatmap",
                                                                 data = "v.csv"))),
               "correlation_heatmap: the variance of column 'big'")
})

# A group with no spread gets no density curve and no bandwidth, so it does not
# stretch the shared x axis (F125).
test_that("a group with no spread gets no density curve", {
  df <- data.frame(condition = c(rep("ctrl", 6), "pair", "pair"),
                   value = c(11.4, 13.1, 12.2, 10.9, 12.8, 11.7, 12.0, 12.0))
  res <- recipe_distribution(df, dist_spec(y = "value", x = "condition"))
  st <- res$stats
  expect_true(is.na(st$bandwidth[st$group == "pair"]))
  expect_false(is.na(st$bandwidth[st$group == "ctrl"]))
  expect_match(res$methods, "The group(s) pair have no spread", fixed = TRUE)
  expect_null(.dist_group_curves("pair", data.frame(group = df$condition, value = df$value),
                                 .dist_bins(df$value)))
})

# A short number such as 1.3 has no padding in the subtitle (F127).
test_that("the distribution subtitle has one space before the bin width", {
  v <- (1:8) * 1.3 / 3.5
  df <- data.frame(condition = rep(c("a", "b"), 4), value = v)
  res <- recipe_distribution(df, dist_spec(y = "value", x = "condition"))
  expect_equal(res$stats$bin_width[1], 1.3, tolerance = 1e-9)
  expect_match(res$plot$labels$subtitle, "bin width 1.3$")
  expect_false(grepl("  ", res$methods, fixed = TRUE))
})
