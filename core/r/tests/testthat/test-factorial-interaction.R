# The interaction geom of factorial_anova and the label wrapping of its x axis.

fa_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "factorial_anova", data = "d.csv", test = "anova"),
                                  list(...)))
}

fa_data <- function() {
  set.seed(9)
  df <- expand.grid(rep = 1:5, dose = c("Vehicle control", "High dose (10 mg/kg)"),
                    strain = c("C57BL/6J wild type", "NOD scid gamma (NSG)", "FVB/N"),
                    stringsAsFactors = FALSE)
  df$y <- 10 + stats::rnorm(nrow(df)) + ifelse(df$dose == "Vehicle control", 0, 3)
  df
}

test_that("the interaction geom builds and keeps the statistics of the box geom", {
  df <- fa_data()
  box <- recipe_factorial_anova(df, fa_spec(x = "strain", fill = "dose", y = "y"))
  int <- recipe_factorial_anova(df, fa_spec(x = "strain", fill = "dose", y = "y", geom = "interaction"))
  expect_equal(int$stats, box$stats)
  built <- ggplot2::ggplot_build(int$plot)
  expect_equal(nrow(built$data[[3]]), 6)
  expect_true(inherits(int$plot$layers[[1]]$geom, "GeomLine"))
  script <- call_build_script(int$build_script, "input.csv", "figure")
  expect_true(any(grepl("position_dodge(width = 0.25)", script, fixed = TRUE)))
})

test_that("the cell summary gives the mean and SEM of each cell", {
  df <- fa_data()
  df$strain <- factor(df$strain); df$dose <- factor(df$dose)
  s <- .fa_cell_summary(df, "y", c("strain", "dose"))
  cell <- df$y[df$strain == "FVB/N" & df$dose == "Vehicle control"]
  r <- s[s$strain == "FVB/N" & s$dose == "Vehicle control", ]
  expect_equal(r$mean, mean(cell))
  expect_equal(r$sem, stats::sd(cell) / sqrt(5))
})

test_that("long x levels wrap at 12 characters and the figure grows with the design", {
  expect_equal(.fa_wrap_labels("NOD scid gamma (NSG)"), "NOD scid\ngamma (NSG)")
  expect_false(.fa_needs_wrap(c("wt", "ko")))
  small <- .fa_fig_size(c("wt", "ko"), c("vehicle", "low", "high"), "treatment", 1, FALSE)
  expect_equal(small$width, 5.2)
  expect_equal(small$height, 4.4)
  big <- .fa_fig_size(sprintf("strain %d", 1:4), c("a", "b", "c"), "treatment", 1, FALSE)
  expect_gt(big$width, 5.2)
})
