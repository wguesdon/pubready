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

# A class with no feature gets no legend key, as in the Python engine.
test_that("the volcano legend lists only the classes that hold a feature", {
  cat <- factor(c("NS", "FC and p", "p"), levels = c("NS", "FC", "p", "FC and p"))
  expect_equal(.volcano_present(cat), c("NS", "p", "FC and p"))
})

# The methods text names the p column of the cutoff and the label rule.
test_that("the volcano methods name the p column and the label rule", {
  df <- data.frame(gene = paste0("g", 1:6),
                   log2FoldChange = c(3, -3, 0.1, 2.5, -0.2, 0),
                   padj = c(1e-8, 1e-9, 0.5, 1e-6, 0.3, 0.9), stringsAsFactors = FALSE)
  spec <- spec_from_opt(list(recipe = "volcano", data = "de.csv", y = "padj", top_n = 15))
  res <- recipe_volcano(df, spec)
  expect_match(res$methods, "|log2 fold change| >= 1 and padj <= 0.05", fixed = TRUE)
  expect_match(res$methods, "The 3 significant features with the largest product", fixed = TRUE)
})
