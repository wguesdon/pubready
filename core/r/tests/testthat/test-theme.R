# The shared label and axis helpers. The Python suite (test_theme.py) expects
# the same strings.

test_that("wrap_label breaks a long label at the spaces, as the Python engine does", {
  expect_equal(wrap_label("Vehicle control (DMSO 0.1%)"), "Vehicle\ncontrol\n(DMSO 0.1%)")
  expect_equal(wrap_label("Compound 7b, 10 uM (72 h)"), "Compound 7b,\n10 uM (72 h)")
  expect_equal(wrap_label("Overexpression"), "Overexpression")
  expect_equal(wrap_label("WT"), "WT")
  expect_equal(wrap_labels(c("KO + rescue", "High dose + adjuvant")),
               c("KO + rescue", "High dose +\nadjuvant"))
  expect_equal(longest_line("High dose +\nadjuvant"), 11L)
})

test_that("plain_number_labels prints plain numbers with a thousands separator", {
  expect_equal(plain_number_labels(c(-2e5, 0, 2e5, 4e5)), c("-200,000", "0", "200,000", "400,000"))
  expect_equal(plain_number_labels(c(0, 0.25, 0.5)), c("0.00", "0.25", "0.50"))
  expect_equal(plain_number_labels(c(1000, 1500, 2000)), c("1,000", "1,500", "2,000"))
  expect_equal(plain_number_labels(c(10, 20, NA)), c("10", "20", NA))
})

test_that("category_width grows with the number of groups", {
  expect_equal(category_width(4), 5.4)
  expect_equal(category_width(8), 9.8)
  expect_equal(category_width(2), 4.2)
})

test_that("category_width widens the slots for a long unbreakable word", {
  expect_equal(category_width(8, wrap_labels(c("WT", "Overexpression"))), 0.085 * 14 * 8 + 1)
})

test_that("rain_offsets puts the raincloud points left of the box, as Python does", {
  off <- rain_offsets(3)
  expect_equal(off, c(-0.1058359, -0.1516718, -0.0775078), tolerance = 1e-6)
  expect_true(all(rain_offsets(500) >= -0.18 & rain_offsets(500) <= -0.06))
})

test_that("geom_top reaches past the data for an untrimmed violin", {
  v <- c(1, 2, 3, 4, 10)
  expect_equal(geom_top(v, "box"), 10)
  expect_equal(geom_top(v, "violin"), 10 + 3 * stats::bw.nrd0(v))
  expect_equal(geom_top(c(1, 1, 1, 1), "bar"), 1)
})

.is_point_layer <- function(layer) inherits(layer$geom, "GeomPoint")

test_that("the comparison points sit at fixed offsets, the same as in Python", {
  df <- data.frame(group = rep(c("a", "b"), each = 3), value = c(1, 2, 3, 4, 5, 6))
  spec <- spec_from_opt(list(recipe = "two_group_compare", data = "d.csv", x = "group",
                             y = "value", test = "welch_t"))
  res <- recipe_two_group_compare(df, spec)
  i <- which(vapply(res$plot$layers, .is_point_layer, logical(1)))[1]
  x <- as.numeric(ggplot2::ggplot_build(res$plot)$data[[i]]$x)
  expect_equal(x, c(1.0283282, 0.9366563, 1.0849845, 2.0283282, 1.9366563, 2.0849845),
               tolerance = 1e-6)
})
