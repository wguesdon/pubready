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

# The figure shows the top_n largest intersections, and ties sort by degree and
# then by set name, the same order as the Python engine (test_upset.py).
test_that("recipe_upset caps the drawn intersections at top_n", {
  df <- data.frame(
    gene = paste0("g", 1:9),
    A = c(1, 1, 1, 0, 0, 1, 0, 1, 0),
    B = c(0, 0, 0, 1, 1, 1, 0, 1, 0),
    C = c(0, 0, 0, 0, 0, 0, 1, 1, 1))
  spec <- spec_from_opt(list(recipe = "upset", data = "sets.csv", top_n = 3))
  res <- recipe_upset(df, spec)
  inter <- res$stats[res$stats$type == "intersection", ]
  expect_equal(inter$sets, c("A", "B", "C", "A & B", "A & B & C"))
  expect_equal(inter$size, c(3L, 2L, 2L, 1L, 1L))
  expect_equal(res$test_meta$n_intersections, 5)
  expect_equal(res$test_meta$n_intersections_shown, 3)
  expect_match(.up_subtitle(3, 5), "The 3 largest of 5 intersections; 2 smaller ones")
  expect_null(.up_subtitle(5, 5))
})
