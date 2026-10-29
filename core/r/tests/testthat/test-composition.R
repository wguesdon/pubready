# Behaviour of the composition recipe on small synthetic inputs. The numbers
# themselves are checked against public benchmarks in test-benchmarks.R.

comp_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "composition", data = "d.csv", test = "auto"),
                                  list(...)))
}

comp_data <- function() {
  set.seed(7)
  samples <- sprintf("S%d", 1:10)
  cond <- rep(c("control", "treated"), each = 5)
  cells <- c("T", "B", "NK")
  df <- expand.grid(cell = cells, sample = samples, stringsAsFactors = FALSE)
  df$condition <- cond[match(df$sample, samples)]
  df$n <- stats::rpois(nrow(df), ifelse(df$cell == "B" & df$condition == "treated", 80, 40))
  df
}

test_that("composition turns counts into proportions and compares two groups", {
  df <- comp_data()
  res <- recipe_composition(df, comp_spec(id = "sample", x = "cell", y = "n", group = "condition"))
  long <- res$df_used
  expect_equal(as.numeric(tapply(long$proportion, long$sample, sum)), rep(1, 10), tolerance = 1e-12)
  props <- long$proportion[long$category == "B"]
  grp <- long$group[long$category == "B"]
  wt <- suppressWarnings(stats::wilcox.test(props[grp == "control"], props[grp == "treated"]))
  b <- res$stats[res$stats$category == "B", ]
  expect_equal(b$p_value[1], wt$p.value, tolerance = 1e-12)
  expect_equal(b$statistic[1], unname(wt$statistic))
  expect_equal(unique(res$stats$test), "Mann-Whitney U test")
  p <- res$stats$p_value[!duplicated(res$stats$category)]
  expect_equal(res$stats$p_adj[!duplicated(res$stats$category)], stats::p.adjust(p, "BH"))
  expect_true(grepl("compositional", res$methods))
})

test_that("composition orders the samples by group and then by natural name", {
  df <- comp_data()
  df$sample[df$sample == "S1"] <- "S10b"
  res <- recipe_composition(df, comp_spec(id = "sample", x = "cell", y = "n", group = "condition"))
  ids <- unique(res$df_used$sample)
  expect_equal(ids[1:4], c("S2", "S3", "S4", "S5"))
  expect_equal(.comp_natural_order(c("S10", "S2", "S1")), c(3L, 2L, 1L))
})

test_that("composition counts rows without --y and uses Kruskal-Wallis for three groups", {
  set.seed(8)
  df <- data.frame(sample = rep(sprintf("s%02d", 1:9), each = 30),
                   cell = sample(c("x", "y", "z"), 270, replace = TRUE))
  df$grp <- rep(c("a", "b", "c"), each = 90)
  res <- recipe_composition(df, comp_spec(id = "sample", x = "cell", group = "grp"))
  expect_equal(unique(res$stats$test), "Kruskal-Wallis test")
  expect_equal(sum(res$df_used$count), 270)
  expect_true(any(grepl("counts as one observation", res$clean_steps)))
})

test_that("composition without --group is descriptive and refuses a sample in two groups", {
  df <- comp_data()
  res <- recipe_composition(df, comp_spec(id = "sample", x = "cell", y = "n"))
  expect_true(all(is.na(res$stats$p_value)))
  expect_equal(res$test_meta$name, "descriptive (no test)")
  df$condition[1] <- "treated"
  expect_error(recipe_composition(df, comp_spec(id = "sample", x = "cell", y = "n", group = "condition")),
               "one group")
})
