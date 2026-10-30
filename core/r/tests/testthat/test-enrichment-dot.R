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

test_that("the stats table lists the terms in the order of the adjusted p", {
  # The gene ratio runs against the adjusted p here, so the plot order and the
  # table order differ; the Python engine writes the same table order.
  df <- data.frame(
    Description = paste0("term", 1:4),
    GeneRatio = c("8/100", "6/100", "4/100", "2/100"),
    Count = c(8, 6, 4, 2),
    p.adjust = c(1e-6, 1e-5, 1e-4, 1e-3),
    stringsAsFactors = FALSE, check.names = FALSE)
  spec <- spec_from_opt(list(recipe = "enrichment_dot", data = "e.csv", top_n = 15))

  res <- recipe_enrichment_dot(df, spec)

  expect_equal(res$stats$term, paste0("term", 1:4))
  expect_equal(res$stats$p_adjust, c(1e-6, 1e-5, 1e-4, 1e-3))
})

# Long GO terms wrap into balanced lines of at most 45 characters; the Python
# engine has the same cases in test_enrichment_dot.py.
test_that(".enrich_wrap balances long labels", {
  expect_equal(.enrich_wrap("glycolytic process"), "glycolytic process")
  expect_equal(.enrich_wrap("regulation of transcription by RNA polymerase II"),
               "regulation of transcription\nby RNA polymerase II")
  expect_equal(.enrich_wrap(paste("negative regulation of extrinsic apoptotic signaling",
                                  "pathway via death domain receptors")),
               paste0("negative regulation of extrinsic apoptotic\n",
                      "signaling pathway via death domain receptors"))
  size <- .enrich_size(c("a\nb", "short"))
  expect_equal(size$height, 4.0)
})

# The enrichment table of the area D tests: one term with an adjusted p of 0,
# one with no adjusted p and one with no name. The Python suite uses the same.
enrich_area_d <- function() {
  data.frame(Description = c("Pathway 0", "Pathway 1", NA, "Pathway 3", "Pathway 4"),
             GeneRatio = c("5/200", "6/200", "7/200", "8/200", "9/200"),
             Count = c(5, 6, 7, 8, 9), p.adjust = c(0, NA, 1e-4, 1e-5, 1e-6),
             check.names = FALSE, stringsAsFactors = FALSE)
}

# An adjusted p of 0 is kept and marked, and each dropped row is logged (F51).
test_that("enrichment_dot keeps an adjusted p of 0 and logs the dropped rows", {
  res <- recipe_enrichment_dot(enrich_area_d(),
                               spec_from_opt(list(recipe = "enrichment_dot", data = "e.csv")))
  expect_equal(res$test_meta$n_terms_total, 3)
  expect_equal(res$test_meta$n_zero_p, 1)
  expect_equal(res$stats$term[1], "Pathway 0")
  expect_equal(res$stats$p_adjust[1], .Machine$double.xmin)
  steps <- paste(res$clean_steps, collapse = " | ")
  expect_match(steps, "Dropped 2 row(s) with a missing term, adjusted p or x value (5 -> 3).", fixed = TRUE)
  expect_match(steps, "Replaced 1 adjusted p value(s) of exactly 0", fixed = TRUE)
  expect_match(res$methods, "drawn as triangles", fixed = TRUE)
})
