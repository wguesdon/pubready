# Behaviour of the four omics recipes (ma_plot, manhattan, feature_panel,
# embedding) on small synthetic inputs. The numbers themselves are checked
# against public benchmarks in test-benchmarks.R.

omics_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

test_that("ma_plot reads a limma table, keeps a missing p as not significant and sorts by p", {
  df <- data.frame(ID = paste0("g", 1:6), logFC = c(2, -3, 0.2, 1.5, -1.2, 4),
                   AveExpr = c(5, 6, 7, 8, 9, 10),
                   adj.P.Val = c(1e-4, 1e-6, 0.5, 0.2, 0.01, NA), stringsAsFactors = FALSE)
  res <- recipe_ma_plot(df, omics_spec("ma_plot"))
  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$test_meta$mean_column, "AveExpr")
  expect_equal(res$test_meta$p_column, "adj.P.Val")
  expect_false(res$test_meta$mean_log10)
  expect_equal(res$test_meta$n_up, 1)
  expect_equal(res$test_meta$n_down, 2)
  expect_equal(res$test_meta$n_not_significant, 3)
  expect_equal(res$test_meta$n_missing_p, 1)
  expect_equal(res$stats$gene, c("g2", "g1", "g5"))
  expect_equal(res$stats$A, c(6, 5, 9))
})

test_that("ma_plot draws log10 of a DESeq2 baseMean and drops a mean of zero", {
  df <- data.frame(gene = c("a", "b", "c"), baseMean = c(100, 0, 1000),
                   log2FoldChange = c(1.5, 2, -0.1), padj = c(0.01, 0.01, 0.9))
  res <- recipe_ma_plot(df, omics_spec("ma_plot"))
  expect_true(res$test_meta$mean_log10)
  expect_equal(res$test_meta$n_total, 2)
  expect_equal(res$stats$A, 2)
  expect_true(any(grepl("baseMean of zero", res$clean_steps)))
})

test_that("ma_plot needs a mean column", {
  df <- data.frame(gene = "a", log2FoldChange = 1, padj = 0.01)
  expect_error(recipe_ma_plot(df, omics_spec("ma_plot")), "no mean-expression column")
})

test_that("manhattan puts the chromosomes in natural order and labels one lead per chromosome", {
  df <- data.frame(snp = paste0("rs", 1:8),
                   chr = c("chrX", "chr10", "chr2", "chr2", "chrMT", "chr1", "chr10", "chr22"),
                   pos = c(5, 10, 20, 30, 1, 7, 40, 3),
                   p = c(0.5, 1e-9, 2e-10, 3e-12, 0.9, 0.04, 2e-9, 0.3))
  res <- recipe_manhattan(df, omics_spec("manhattan", group = "chr", x = "pos", y = "p", label = "snp"))
  expect_equal(res$test_meta$n_chromosomes, 6)
  leads <- res$stats[res$stats$row_type == "lead", ]
  expect_equal(leads$snp, c("rs4", "rs2"))
  expect_equal(leads$chr, c("2", "10"))
  expect_equal(.mh_chr_order(c("X", "10", "2", "MT", "1", "22", "Y")),
               c("1", "2", "10", "22", "X", "Y", "MT"))
  expect_equal(res$test_meta$n_genome_wide, 4)
  expect_true(any(grepl("'chr' prefix", res$clean_steps)))
})

test_that("manhattan sets a p of exactly 0 to the smallest double and rejects p above 1", {
  df <- data.frame(chr = c(1, 1, 2, 2), pos = 1:4, p = c(0, 0.5, 0.2, 0.9))
  res <- recipe_manhattan(df, omics_spec("manhattan", group = "chr", x = "pos", y = "p"))
  expect_equal(res$test_meta$n_zero_p, 1)
  expect_equal(res$test_meta$min_p, .Machine$double.xmin)
  expect_true(any(grepl("exactly 0", res$clean_steps)))
  bad <- data.frame(chr = 1, pos = 1:2, p = c(0.5, 1.5))
  expect_error(recipe_manhattan(bad, omics_spec("manhattan", group = "chr", x = "pos", y = "p")),
               "outside 0 to 1")
})

test_that("manhattan thins only the points with p >= 0.01 and computes lambda on all", {
  set.seed(7)
  n <- 60000
  df <- data.frame(chr = rep(1:3, length.out = n), pos = seq_len(n), p = stats::runif(n))
  res <- recipe_manhattan(df, omics_spec("manhattan", group = "chr", x = "pos", y = "p"))
  keep <- .mh_thin(df$p)
  expect_true(all(keep$keep[df$p < 0.01]))
  expect_gt(res$test_meta$thinning_step, 1)
  expect_lte(res$test_meta$n_drawn, 51000)
  expect_equal(res$test_meta$lambda_gc,
               stats::median(stats::qchisq(1 - df$p, 1)) / stats::qchisq(0.5, 1), tolerance = 1e-10)
})

test_that("feature_panel runs the Welch t-test on two normal groups and adjusts with BH", {
  set.seed(11)
  df <- data.frame(group = rep(c("ctrl", "trt"), each = 12),
                   f1 = stats::rnorm(24), f2 = stats::rnorm(24) + rep(c(0, 2), each = 12),
                   f3 = stats::rnorm(24))
  res <- recipe_feature_panel(df, omics_spec("feature_panel", x = "group", test = "welch_t"))
  tt <- stats::t.test(df$f2[1:12], df$f2[13:24])
  r <- res$stats[res$stats$feature == "f2", ]
  expect_equal(r$statistic, unname(tt$statistic), tolerance = 1e-10)
  expect_equal(r$df1, unname(tt$parameter), tolerance = 1e-10)
  expect_equal(res$stats$p_adj, stats::p.adjust(res$stats$p_value, "BH"))
  expect_equal(res$stats[["n_ctrl"]], c(12L, 12L, 12L))
  expect_equal(r$effect_size_name, "Cohen's d")
})

test_that("feature_panel reads long data and applies R's wilcox rule", {
  df <- data.frame(gene = rep(c("A", "B"), each = 10), arm = rep(rep(c("x", "y"), each = 5), 2),
                   value = c(1:5, 6:10, c(1, 3, 5, 7, 9), c(2, 4, 6, 8, 10)))
  res <- recipe_feature_panel(df, omics_spec("feature_panel", x = "arm", y = "value", facet = "gene",
                                             test = "wilcoxon"))
  wt <- stats::wilcox.test(c(1, 3, 5, 7, 9), c(2, 4, 6, 8, 10))
  expect_equal(res$stats$p_value[res$stats$feature == "B"], wt$p.value)
  expect_equal(res$stats$feature, c("A", "B"))
})

test_that("feature_panel refuses more than 36 features and a two-group test on three groups", {
  wide <- as.data.frame(matrix(stats::rnorm(37 * 6), nrow = 6))
  wide$group <- rep(c("a", "b"), 3)
  expect_error(recipe_feature_panel(wide, omics_spec("feature_panel", x = "group")), "heatmap")
  three <- data.frame(group = rep(c("a", "b", "c"), each = 4), v = stats::rnorm(12))
  expect_error(recipe_feature_panel(three, omics_spec("feature_panel", x = "group", test = "welch_t")),
               "exactly 2 groups")
})

test_that("embedding silhouettes match cluster::silhouette on a 2D layout", {
  set.seed(5)
  df <- data.frame(u1 = c(stats::rnorm(40), stats::rnorm(30, 4)), u2 = c(stats::rnorm(40), stats::rnorm(30, 3)),
                   cl = rep(c("b", "a"), c(40, 30)))
  res <- recipe_embedding(df, omics_spec("embedding", x = "u1", y = "u2", group = "cl"))
  s <- summary(cluster::silhouette(as.integer(factor(df$cl)), stats::dist(df[, 1:2])))
  expect_equal(res$stats$cluster, c("a", "b"))
  expect_equal(res$stats$mean_silhouette, as.numeric(s$clus.avg.widths), tolerance = 1e-12)
  expect_equal(res$test_meta$mean_silhouette, s$avg.width, tolerance = 1e-12)
  expect_equal(sum(res$stats$percent), 100)
})

test_that("embedding treats whole-number cluster ids as clusters and a score as continuous", {
  set.seed(6)
  ids <- data.frame(x = stats::rnorm(200), y = stats::rnorm(200), k = rep(0:19, 10))
  res <- recipe_embedding(ids, omics_spec("embedding", x = "x", y = "y", group = "k"))
  expect_false(res$test_meta$continuous)
  expect_equal(res$stats$cluster[1:3], c("0", "1", "2"))
  score <- data.frame(x = stats::rnorm(200), y = stats::rnorm(200), s = stats::runif(200))
  res2 <- recipe_embedding(score, omics_spec("embedding", x = "x", y = "y", group = "s"))
  expect_true(res2$test_meta$continuous)
  expect_true(is.na(res2$stats$mean_silhouette))
})

test_that("embedding computes the silhouette on every k-th cell above 5000 cells", {
  set.seed(8)
  n <- 12001
  df <- data.frame(x = stats::rnorm(n), y = stats::rnorm(n), g = sample(c("p", "q"), n, TRUE))
  res <- recipe_embedding(df, omics_spec("embedding", x = "x", y = "y", group = "g"))
  expect_equal(res$test_meta$silhouette_step, 3)
  expect_equal(res$test_meta$silhouette_n, length(seq(1, n, by = 3)))
  expect_match(res$methods, "one in every 3 cells")
})

test_that("feature_panel reports the Kruskal eta-squared (H) of rstatix", {
  df <- data.frame(group = rep(c("a", "b", "c"), each = 8),
                   f1 = c(1:8, (1:8) + 3, (1:8) + 6), f2 = c(8:1, 1:8, c(2, 4, 1, 7, 3, 8, 5, 6)))
  res <- recipe_feature_panel(df, omics_spec("feature_panel", x = "group", test = "kruskal"))
  long <- data.frame(group = df$group, v = df$f1)
  ref <- unname(rstatix::kruskal_effsize(long, v ~ group)$effsize)
  expect_equal(res$stats$effect_size[res$stats$feature == "f1"], ref, tolerance = 1e-10)
  expect_equal(unique(res$stats$effect_size_name), "eta-squared (H)")
})

# A numeric chromosome column with an empty cell reads "1", not "1.0" (F64).
test_that("manhattan writes a numeric chromosome as a whole number", {
  expect_equal(.mh_chr_name(c(1, 10, NA, 2)), c("1", "10", NA, "2"))
  expect_equal(.mh_chr_name(c("chr1", "chrX")), c("1", "X"))
})
