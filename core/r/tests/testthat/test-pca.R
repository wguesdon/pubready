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

# The features leave out --id and --label, --features sets them, and the data
# log names each feature (F47).
pca_id_data <- function() {
  set.seed(2)
  data.frame(sample_id = 1:12, group = rep(c("A", "B", "C"), each = 4),
             f1 = rnorm(12), f2 = rnorm(12), f3 = rnorm(12), f4 = rnorm(12))
}

test_that("pca leaves out --id, takes --features and names the features", {
  df <- pca_id_data()
  plain <- recipe_pca(df, spec_from_opt(list(recipe = "pca", data = "p.csv", group = "group")))
  expect_equal(plain$test_meta$n_features, 5)
  expect_match(plain$clean_steps[1], "'sample_id', 'f1'", fixed = TRUE)
  res <- recipe_pca(df, spec_from_opt(list(recipe = "pca", data = "p.csv", group = "group",
                                           id = "sample_id")))
  expect_equal(res$test_meta$n_features, 4)
  expect_match(res$clean_steps[1], "4 numeric features: 'f1', 'f2', 'f3', 'f4'.", fixed = TRUE)
  expect_true(any(grepl("Left out the --id and --label column(s) 'sample_id'", res$clean_steps,
                        fixed = TRUE)))
  sel <- recipe_pca(df, spec_from_opt(list(recipe = "pca", data = "p.csv", group = "group",
                                           features = "f1,f3,f4")))
  expect_equal(sel$test_meta$n_features, 3)
  expect_error(recipe_pca(df, spec_from_opt(list(recipe = "pca", data = "p.csv", group = "group",
                                                 features = "f1,f9"))), "'f9'")
})
