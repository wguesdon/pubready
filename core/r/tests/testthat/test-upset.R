# The UpSet recipe: elements in no set are dropped, and the stats table carries
# both set sizes and intersection sizes.
test_that("recipe_upset summarizes set and intersection sizes", {
  df <- data.frame(
    gene = paste0("g", 1:6),
    A = c(1, 1, 0, 1, 0, 0),
    B = c(1, 0, 1, 1, 0, 0),
    C = c(0, 0, 1, 1, 1, 0)
  )
  opt <- list(recipe = "upset", data = "sets.csv")
  spec <- spec_from_opt(opt)

  res <- recipe_upset(df, spec)

  expect_true(is.function(res$draw))
  expect_null(res$plot)
  sets <- res$stats[res$stats$type == "set", ]
  expect_equal(sets$size[sets$sets == "A"], 3)
  # g6 belongs to no set, so it is dropped: 5 elements remain.
  expect_equal(res$test_meta$n_elements, 5)
  expect_true(any(res$stats$type == "intersection"))
})
