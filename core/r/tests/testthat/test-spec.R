test_that("spec_from_opt builds the expected spec with house defaults", {
  opt <- list(recipe = "two_group_compare", data = "/path/to/tumor_volume.csv",
              x = "group", y = "volume", test = "auto", paired = FALSE,
              p_adjust = "none")
  spec <- spec_from_opt(opt)

  expect_equal(spec$engine, "r")
  expect_equal(spec$recipe, "two_group_compare")
  expect_equal(spec$data$file, "tumor_volume.csv")
  expect_equal(spec$data$x, "group")
  expect_equal(spec$data$y, "volume")
  expect_equal(spec$test$method, "auto")
  expect_false(spec$test$paired)
  expect_equal(spec$appearance$geom, "box")
  expect_true(spec$appearance$show_points)
})

test_that("spec_from_opt splits palette and covariates on commas", {
  opt <- list(recipe = "cox_forest", data = "trial.csv", palette = "#3B6DB3, #C1432B",
              covariates = "age, sex, stage", test = "auto")
  spec <- spec_from_opt(opt)
  expect_equal(spec$appearance$palette, c("#3B6DB3", "#C1432B"))
  expect_equal(spec$data$covariates, c("age", "sex", "stage"))
})

test_that("spec_from_opt carries the theme and defaults to the house style", {
  opt <- list(recipe = "two_group_compare", data = "d.csv", x = "g", y = "v",
              theme = "prism")
  expect_equal(spec_from_opt(opt)$appearance$theme, "prism")
  opt2 <- list(recipe = "two_group_compare", data = "d.csv", x = "g", y = "v")
  expect_equal(spec_from_opt(opt2)$appearance$theme, "pubplot_house")
})

test_that("fix_axis_keys restores YAML-coerced boolean axis keys", {
  d <- list("TRUE" = "volume", x = "group")
  fixed <- fix_axis_keys(d)
  expect_true("y" %in% names(fixed))
  expect_false("TRUE" %in% names(fixed))
  expect_equal(fixed$y, "volume")
})
