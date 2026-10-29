# multi_group_compare: the bracket rule and the compact letter display. The
# Python suite (test_multi_group.py) uses the same inputs and expects the same
# letters, so the two engines label the groups identically.

cld_pairs <- function() {
  g <- utils::combn(c("A", "B", "C", "D", "E"), 2)
  sig <- c("A-C", "A-D", "A-E", "B-D", "B-E", "C-E")
  data.frame(group1 = g[1, ], group2 = g[2, ],
             p = ifelse(paste(g[1, ], g[2, ], sep = "-") %in% sig, 0.001, 0.4),
             stringsAsFactors = FALSE)
}

test_that("cld_letters follows the insert-and-absorb algorithm", {
  pr <- cld_pairs()
  expect_equal(unname(cld_letters(LETTERS[1:5], pr$group1, pr$group2, pr$p)),
               c("a", "ab", "bc", "cd", "d"))
  expect_equal(unname(cld_letters(LETTERS[1:5], pr$group1, pr$group2, rep(0.5, 10))),
               rep("a", 5))
  expect_equal(unname(cld_letters(c("x", "y", "z"), c("x", "x", "y"), c("y", "z", "z"),
                                  c(0.01, 0.01, 0.01))),
               c("a", "b", "c"))
})

# Six groups in three pairs of equal means, with the same fixed noise in both
# engines: 12 pairs differ, so the figure shows letters.
six_groups <- function() {
  g <- rep(sprintf("G%d", 1:6), each = 6)
  mu <- rep(c(0, 0, 5, 5, 10, 10), each = 6)
  data.frame(group = g, value = mu + 0.3 * sin(seq_along(g) * 1.3), stringsAsFactors = FALSE)
}

test_that("multi_group_compare switches to letters above six significant pairs", {
  spec <- spec_from_opt(list(recipe = "multi_group_compare", data = "d.csv", x = "group",
                             y = "value", test = "anova"))
  res <- recipe_multi_group_compare(six_groups(), spec)
  expect_equal(res$test_meta$n_significant, 12L)
  expect_equal(res$test_meta$display, "compact letter display")
  expect_equal(unlist(res$test_meta$letters),
               c(G1 = "a", G2 = "a", G3 = "b", G4 = "b", G5 = "c", G6 = "c"))
  expect_match(res$methods, "compact letter display")
  expect_equal(res$width, 7.6)
})

test_that("multi_group_compare keeps brackets for a few pairs, one level each", {
  df <- six_groups()[1:18, ]
  spec <- spec_from_opt(list(recipe = "multi_group_compare", data = "d.csv", x = "group",
                             y = "value", test = "anova"))
  res <- recipe_multi_group_compare(df, spec)
  expect_equal(res$test_meta$display, "brackets")
  rows <- .mg_bracket_rows(data.frame(group1 = c("G1", "G1"), group2 = c("G3", "G2"),
                                      p.adj = c(1e-5, 0.005), stringsAsFactors = FALSE),
                           c("G1", "G2", "G3"), top_all = 10, span = 10, bracket_label = "p.signif")
  expect_equal(rows$group2, c("G2", "G3"))
  expect_equal(rows$y.position, c(10.8, 12.1))
  expect_equal(rows$lab, c("**", "****"))
})

test_that("multi_group_compare names the direction of each estimate", {
  spec <- spec_from_opt(list(recipe = "multi_group_compare", data = "d.csv", x = "group",
                             y = "value", test = "anova"))
  res <- recipe_multi_group_compare(six_groups(), spec)
  expect_equal(res$stats$contrast, paste(res$stats$group2, "-", res$stats$group1))
})
