# The spider recipe: every value becomes a percentage change from the patient's
# own baseline, and the two RECIST thresholds put each patient in a category.

sample_spider <- function() {
  data.frame(
    patient = rep(c("p1", "p2", "p3"), each = 3),
    arm     = rep(c("treated", "control", "treated"), each = 3),
    week    = rep(c(0, 8, 16), times = 3),
    size    = c(100, 80, 50,     # p1 falls to -50%: partial response
                100, 110, 130,   # p2 rises to +30%: progressive disease
                100, 95, 105),   # p3 stays inside both lines: stable
    stringsAsFactors = FALSE)
}

spider_spec <- function(...) {
  spec_from_opt(utils::modifyList(
    list(recipe = "spider_response", data = "spider.csv",
         x = "week", y = "size", id = "patient"), list(...)))
}

test_that("the change is measured from the first visit of each patient", {
  res <- recipe_spider_response(sample_spider(), spider_spec())

  first <- res$df_used[res$df_used$week == 0, ]
  expect_equal(unique(first$change), 0)
  expect_equal(res$test_meta$n_patients, 3)
  expect_equal(res$test_meta$n_measurements, 9)
})

test_that("the two thresholds put each patient in a category", {
  res <- recipe_spider_response(sample_spider(), spider_spec())

  got <- setNames(res$stats$category_from_change, res$stats$patient)
  expect_equal(unname(got["p1"]), "partial response")
  expect_equal(unname(got["p2"]), "progressive disease")
  expect_equal(unname(got["p3"]), "stable disease")
  expect_equal(res$test_meta$partial_response, 1)
  expect_equal(res$test_meta$progressive_disease, 1)
  expect_equal(res$test_meta$stable_disease, 1)
})

test_that("a threshold the reader moves changes the category", {
  # p3 reaches +5% and -5%. A pair of thresholds inside that range makes it
  # progressive disease, which is what says the parameters are read.
  res <- recipe_spider_response(sample_spider(),
                                spider_spec(pd_threshold = 4, pr_threshold = -4))

  got <- setNames(res$stats$category_from_change, res$stats$patient)
  expect_equal(unname(got["p3"]), "partial response")
  expect_equal(res$test_meta$pd_threshold, 4)
  expect_equal(res$test_meta$pr_threshold, -4)
})

test_that("a patient with a zero baseline is dropped and recorded", {
  df <- rbind(sample_spider(),
              data.frame(patient = "p4", arm = "control", week = c(0, 8, 16),
                         size = c(0, 10, 20), stringsAsFactors = FALSE))

  res <- recipe_spider_response(df, spider_spec())

  expect_equal(res$test_meta$n_patients, 3)
  expect_equal(res$test_meta$n_dropped_patients, 1)
  expect_false("p4" %in% res$stats$patient)
  expect_true(any(grepl("non-zero baseline", res$clean_steps)))
})

test_that("the recipe names the column it cannot find", {
  expect_error(recipe_spider_response(sample_spider(), spider_spec(id = "subject")),
               "subject")
  expect_error(recipe_spider_response(sample_spider(),
                                      spec_from_opt(list(recipe = "spider_response",
                                                         data = "spider.csv",
                                                         x = "week", y = "size"))),
               "--id")
})

test_that("the group column colours the lines and reaches the table", {
  res <- recipe_spider_response(sample_spider(), spider_spec(group = "arm"))

  expect_true("group" %in% names(res$stats))
  expect_equal(sort(unique(res$stats$group)), c("control", "treated"))
})

test_that("a patient who only grows carries its smallest growth", {
  # The baseline reads zero for everyone, so it cannot be the best change.
  df <- data.frame(patient = rep("p1", 3), week = c(0, 8, 16),
                   size = c(100, 105, 112), stringsAsFactors = FALSE)

  res <- recipe_spider_response(df, spider_spec())

  expect_equal(res$stats$best_change_pct[1], 5)
  expect_equal(res$stats$worst_change_pct[1], 12)
  expect_equal(res$stats$category_from_change[1], "stable disease")
})

test_that("a patient with no visit after the baseline is dropped", {
  df <- rbind(sample_spider(),
              data.frame(patient = "p5", arm = "treated", week = 0,
                         size = 100, stringsAsFactors = FALSE))

  res <- recipe_spider_response(df, spider_spec())

  expect_equal(res$test_meta$n_without_follow_up, 1)
  expect_false("p5" %in% res$stats$patient)
  expect_true(any(grepl("after the baseline", res$clean_steps)))
})
