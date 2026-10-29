# The heatmap leaf order is the same in both engines. The Python suite
# (test_heatmap.py) uses the same matrix and expects the same order.
#
# g1 and g3 form one branch and g2 and g4 the other. The branch g2, g4 has the
# higher mean, and inside each branch the gene with the higher mean comes first.
# The columns s1 and s3 are equal, as are s2 and s4, so the order inside each
# column branch falls to the tie rule: the earlier input column first.

hm_matrix <- function() {
  data.frame(gene = c("g1", "g2", "g3", "g4"),
             s1 = c(1, 10, 1.5, 10.2), s2 = c(2, 11, 2.5, 11.2),
             s3 = c(1, 10, 1.5, 10.2), s4 = c(2, 11, 2.5, 11.2))
}

hm_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "heatmap", data = "m.csv", scale = "none"),
                                  list(...)))
}

test_that("the heatmap orders each branch by the mean of its unscaled values", {
  res <- recipe_heatmap(hm_matrix(), hm_spec())
  expect_equal(unlist(res$test_meta$row_order), c("g4", "g2", "g3", "g1"))
  expect_equal(unlist(res$test_meta$column_order), c("s2", "s4", "s1", "s3"))
  grDevices::pdf(NULL)
  drawn <- res$draw()
  grDevices::dev.off()
  expect_equal(ComplexHeatmap::row_order(drawn), c(4L, 2L, 3L, 1L))
  expect_equal(ComplexHeatmap::column_order(drawn), c(2L, 4L, 1L, 3L))
})

test_that("an axis that is not clustered keeps the input order", {
  res <- recipe_heatmap(hm_matrix(), hm_spec(cluster = "rows"))
  expect_equal(unlist(res$test_meta$row_order), c("g4", "g2", "g3", "g1"))
  expect_equal(unlist(res$test_meta$column_order), c("s1", "s2", "s3", "s4"))
})

test_that("the emitted script carries the order rule and runs", {
  df <- hm_matrix()
  res <- recipe_heatmap(df, hm_spec())
  script <- call_build_script(res$build_script, "input.csv", "figure")
  expect_true(any(grepl(".hm_order_tree <- function", script, fixed = TRUE)))
  dir <- tempfile("heatmap_script_")
  dir.create(dir)
  utils::write.csv(df, file.path(dir, "input.csv"), row.names = FALSE)
  writeLines(script, file.path(dir, "script.R"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "script.R",
                                  stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  expect_true(file.exists(file.path(dir, "figure.pdf")))
})

test_that("the heatmap methods name the groups, the independence and the count", {
  txt <- .hm_methods("row", "both", 40, 8, c("trt", "untrt"), 38)
  expect_match(txt, "compared between trt and untrt with a Welch two-sample t-test", fixed = TRUE)
  expect_match(txt, "treats the samples as independent", fixed = TRUE)
  expect_match(txt, "38 of the 40 features had an adjusted P below 0.05", fixed = TRUE)
  expect_false(grepl("feature(s)", .hm_methods("row", "both", 40, 8, NULL, NA), fixed = TRUE))
})
