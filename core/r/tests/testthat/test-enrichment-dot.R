# The enrichment dot plot recipe: ORA mode when there is no NES column.
test_that("recipe_enrichment_dot builds an ORA dot plot", {
  df <- data.frame(
    Description = paste0("term", 1:5),
    GeneRatio = c("10/100", "8/100", "6/100", "4/100", "2/100"),
    Count = c(10, 8, 6, 4, 2),
    p.adjust = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
    stringsAsFactors = FALSE, check.names = FALSE)
  opt <- list(recipe = "enrichment_dot", data = "e.csv", top_n = 15)
  spec <- spec_from_opt(opt)

  res <- recipe_enrichment_dot(df, spec)

  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$resolved$label, "ORA dot plot")
  expect_equal(nrow(res$stats), 5)
})
