# The correlation recipe: the method resolver and an end-to-end run that must
# report the forced method, the right n, and a strong, significant coefficient.
test_that("resolve_correlation_method picks the method from normality", {
  expect_equal(resolve_correlation_method("auto", TRUE, TRUE)$method, "pearson")
  expect_equal(resolve_correlation_method("auto", TRUE, FALSE)$method, "spearman")
  expect_equal(resolve_correlation_method("auto", FALSE, TRUE)$method, "spearman")
  expect_equal(resolve_correlation_method("spearman", TRUE, TRUE)$method, "spearman")
  expect_error(resolve_correlation_method("nonsense", TRUE, TRUE),
               "unknown correlation method")
})

test_that("recipe_correlation runs end to end", {
  set.seed(1)
  x <- seq(0, 10, length.out = 40)
  y <- 2 * x + rnorm(40)
  df <- data.frame(x = x, y = y)
  opt <- list(recipe = "correlation", data = "toy.csv", x = "x", y = "y",
              test = "pearson")
  spec <- spec_from_opt(opt)

  res <- recipe_correlation(df, spec)

  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$stats$method, "Pearson correlation")
  expect_equal(res$stats$n, 40)
  expect_true(res$stats$estimate > 0.9)
  expect_true(res$stats$p_value < 1e-6)
})

test_that("the correlation subtitle names the method, the coefficient and P as Python does", {
  expect_equal(cor_label("Pearson correlation", "r", 0.7123, 2.1e-10),
               "Pearson correlation, r = 0.71, P = 2.1e-10")
  expect_equal(cor_label("Spearman rank correlation", "rho", -0.456, 0.0312),
               "Spearman rank correlation, rho = -0.46, P = 0.031")
  df <- data.frame(x = 1:20, y = (1:20) + sin(1:20))
  spec <- spec_from_opt(list(recipe = "correlation", data = "toy.csv", x = "x", y = "y",
                             test = "kendall"))
  res <- recipe_correlation(df, spec)
  expect_match(res$plot$labels$subtitle, "^Kendall rank correlation, tau = ")
})

test_that("a rank correlation scatter gets no fitted line", {
  set.seed(4)
  df <- data.frame(a = 1:30, b = (1:30)^2 + stats::rnorm(30))
  base <- list(recipe = "correlation", data = "c.csv", x = "a", y = "b")
  geoms <- function(res) vapply(res$plot$layers, function(l) class(l$geom)[1], character(1))
  sp <- recipe_correlation(df, spec_from_opt(c(base, test = "spearman")))
  expect_false("GeomSmooth" %in% geoms(sp))
  expect_match(sp$methods, "No line is drawn", fixed = TRUE)
  pe <- recipe_correlation(df, spec_from_opt(c(base, test = "pearson")))
  expect_true("GeomSmooth" %in% geoms(pe))
})
