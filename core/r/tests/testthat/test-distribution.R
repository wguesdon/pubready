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
