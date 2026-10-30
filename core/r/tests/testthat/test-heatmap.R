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

# A null result has no star to draw. The heatmap gets no star column and draws,
# where anno_text stopped on a column of empty strings (F46).
hm_null_matrix <- function() {
  data.frame(gene = paste0("g", 1:4),
             a1 = c(5.1, 3.0, 7.2, 1.1), a2 = c(4.9, 3.3, 6.8, 0.9), a3 = c(5.3, 2.8, 7.0, 1.2),
             b1 = c(5.0, 3.1, 7.1, 1.0), b2 = c(5.2, 2.9, 6.9, 1.3), b3 = c(4.8, 3.2, 7.3, 0.8))
}

hm_null_annotation <- function() {
  path <- tempfile("hm_ann_", fileext = ".csv")
  utils::write.csv(data.frame(sample = c("a1", "a2", "a3", "b1", "b2", "b3"),
                              group = rep(c("a", "b"), each = 3)), path, row.names = FALSE)
  path
}

test_that("a heatmap with no significant feature draws without a star column", {
  res <- recipe_heatmap(hm_null_matrix(), hm_spec(annotation = hm_null_annotation(),
                                                  test = "welch_t"))
  expect_equal(res$test_meta$n_significant, 0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(res$draw())
})

# The Welch test runs only with --test welch_t (F97).
test_that("the heatmap runs no test by default and checks --test welch_t", {
  ann <- hm_null_annotation()
  res <- recipe_heatmap(hm_null_matrix(), hm_spec(annotation = ann))
  expect_equal(names(res$stats), c("feature", "mean", "sd"))
  expect_null(res$test_meta$per_feature_test)
  expect_match(res$methods, "No statistical test was applied.", fixed = TRUE)
  three <- utils::read.csv(ann)
  three$group[1] <- "c"
  utils::write.csv(three, ann, row.names = FALSE)
  expect_error(recipe_heatmap(hm_null_matrix(), hm_spec(annotation = ann, test = "welch_t")),
               "compares two groups")
  expect_error(recipe_heatmap(hm_null_matrix(), hm_spec(annotation = ann, test = "anova")),
               "unknown --test 'anova'")
})

# One rule for a degenerate Welch test in both engines: a feature with no
# variance in either group gets P = NA and stays out of the BH count. A sample
# that the annotation misses, or a group of one sample, stops the run (F87).
test_that("the heatmap leaves a degenerate feature out of the BH count", {
  m <- hm_null_matrix()
  m <- rbind(m, data.frame(gene = "g_floor", a1 = 0, a2 = 0, a3 = 0, b1 = 3, b2 = 3, b3 = 3))
  ann <- hm_null_annotation()
  res <- recipe_heatmap(m, hm_spec(annotation = ann, test = "welch_t"))
  st <- res$stats
  expect_true(is.na(st$p_value[st$feature == "g_floor"]))
  expect_true(is.na(st$p_adj[st$feature == "g_floor"]))
  ok <- st$feature != "g_floor"
  expect_equal(st$p_adj[ok], stats::p.adjust(st$p_value[ok], "BH"))
  expect_equal(res$test_meta$n_untested, 1)
  expect_match(res$methods, "1 feature(s) had no variance in either group", fixed = TRUE)
  short <- utils::read.csv(ann)[-6, ]
  utils::write.csv(short, ann, row.names = FALSE)
  expect_error(recipe_heatmap(m, hm_spec(annotation = ann)), "has no row for the sample\\(s\\) 'b3'")
  one <- data.frame(sample = c("a1", "a2", "a3", "b1", "b2", "b3"),
                    group = c("a", "a", "a", "a", "a", "b"))
  utils::write.csv(one, ann, row.names = FALSE)
  expect_error(recipe_heatmap(m, hm_spec(annotation = ann, test = "welch_t")),
               "group 'b' has 1")
})

# The palette continues across the annotation tracks, so two tracks do not
# share a colour (F121). The Python suite expects the same colours.
test_that("the annotation tracks do not share a colour", {
  ann <- data.frame(group = c("b", "a", "a", "b"), batch = c(1, 2, 1, 2))
  ann$batch <- level_labels(ann$batch)
  cols <- .hm_track_colours(ann)
  expect_equal(unname(cols$group), pubready_palette[1:2])
  expect_equal(names(cols$group), c("a", "b"))
  expect_equal(unname(cols$batch), pubready_palette[3:4])
  expect_equal(names(cols$batch), c("1", "2"))
})
