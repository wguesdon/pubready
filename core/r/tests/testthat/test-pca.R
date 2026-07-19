# The PCA recipe: two well-separated groups give a significant PERMANOVA.
test_that("recipe_pca separates two groups", {
  set.seed(0)
  a <- matrix(rnorm(50, 0, 1), 10, 5)
  b <- matrix(rnorm(50, 6, 1), 10, 5)
  mat <- rbind(a, b); colnames(mat) <- paste0("f", 1:5)
  df <- data.frame(mat, group = rep(c("A", "B"), each = 10), check.names = FALSE)
  opt <- list(recipe = "pca", data = "pca.csv", group = "group")
  spec <- spec_from_opt(opt)

  res <- recipe_pca(df, spec)

  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$test_meta$n_samples, 20)
  expect_equal(res$test_meta$n_features, 5)
  expect_true(res$test_meta$permanova$p_value < 0.05)
})
