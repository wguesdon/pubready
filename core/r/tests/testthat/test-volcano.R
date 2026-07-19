# The volcano recipe: counts of up/down/significant features from the cutoffs.
test_that("recipe_volcano counts up/down significant genes", {
  df <- data.frame(
    gene = paste0("g", 1:6),
    log2FoldChange = c(3, -3, 0.1, 2.5, -0.2, 0),
    pvalue = c(1e-8, 1e-9, 0.5, 1e-6, 0.3, 0.9),
    stringsAsFactors = FALSE)
  opt <- list(recipe = "volcano", data = "de.csv", fc_cutoff = 1, p_cutoff = 0.05, top_n = 15)
  spec <- spec_from_opt(opt)

  res <- recipe_volcano(df, spec)

  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$test_meta$n_up, 2)
  expect_equal(res$test_meta$n_down, 1)
  expect_equal(res$test_meta$n_significant, 3)
})
