# Small shared helpers.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Shapiro-Wilk normality on a numeric vector, guarded for small and large n.
shapiro_safe <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 3)   return(list(test = "shapiro-wilk", n = n, p = NA_real_, normal = NA))
  if (n > 5000) return(list(test = "shapiro-wilk", n = n, p = NA_real_, normal = TRUE))
  p <- tryCatch(stats::shapiro.test(v)$p.value, error = function(e) NA_real_)
  list(test = "shapiro-wilk", n = n, p = p, normal = isTRUE(p > 0.05))
}

# Rank-biserial correlation for two independent groups, from the Mann-Whitney U
# of the first group. Negative when the first group ranks lower. This is
# computed here rather than with rstatix::wilcox_effsize, which needs the coin
# package and returns the magnitude without the sign. The Python engine uses the
# same formula.
rank_biserial_unpaired <- function(u, n1, n2) {
  2 * as.numeric(u) / (n1 * n2) - 1
}

# Rank-biserial correlation for paired differences (Kerby's formula). Zero
# differences are dropped, as in the signed-rank test itself.
rank_biserial_paired <- function(d) {
  d <- d[is.finite(d) & d != 0]
  if (length(d) == 0) return(NA_real_)
  rk <- rank(abs(d))
  (sum(rk[d > 0]) - sum(rk[d < 0])) / sum(rk)
}

# Resolve which two-group test to run from the request and the assumption checks.
resolve_two_group_method <- function(method, all_normal, equal_var, paired) {
  method <- tolower(method %||% "auto")
  if (method == "auto") {
    if (isTRUE(all_normal)) {
      family <- "t"; var_equal <- isTRUE(equal_var)
    } else {
      family <- "wilcoxon"; var_equal <- NA
    }
  } else if (method %in% c("t", "t_test", "student_t", "welch_t")) {
    family <- "t"
    var_equal <- if (method == "student_t") TRUE
                 else if (method == "welch_t") FALSE
                 else isTRUE(equal_var)
  } else if (method %in% c("wilcoxon", "mann_whitney", "wilcox", "mwu")) {
    family <- "wilcoxon"; var_equal <- NA
  } else {
    stop(sprintf("unknown test method '%s'", method))
  }
  label <- if (family == "t") {
    if (paired) "paired t-test"
    else if (isTRUE(var_equal)) "Student's two-sample t-test"
    else "Welch two-sample t-test"
  } else {
    if (paired) "Wilcoxon signed-rank test" else "Mann-Whitney U test"
  }
  list(family = family, var_equal = var_equal, label = label, paired = paired)
}

# Installed versions of the named packages, as a named list of strings.
pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) as.character(utils::packageVersion(p)))
  names(out) <- pkgs
  out
}

# Significance stars with the cutpoints of rstatix::add_significance, as
# util.p_stars of the Python engine: **** below 1e-4, *** below 0.001, ** below
# 0.01, * below 0.05, and ns otherwise or for a missing P. Vectorised.
p_stars <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep("ns", length(p))
  ok <- !is.na(p)
  out[ok & p < 0.05] <- "*"
  out[ok & p < 0.01] <- "**"
  out[ok & p < 1e-3] <- "***"
  out[ok & p < 1e-4] <- "****"
  out
}

# The sort key of one letter column: "0" for a member, "1" for a non-member, so
# that an ascending sort puts the column that starts at the first group first.
.cld_key <- function(col) paste(ifelse(col, "0", "1"), collapse = "")

# The absorb step: drop a column that is a proper subset of another column, a
# column equal to an earlier column, and an empty column.
.cld_absorb <- function(m) {
  k <- ncol(m)
  keep <- rep(TRUE, k)
  for (a in seq_len(k)) {
    if (!any(m[, a])) {
      keep[a] <- FALSE
      next
    }
    for (b in seq_len(k)) {
      if (a == b || !all(m[m[, a], b])) next
      if (any(m[!m[, a], b]) || b < a) {
        keep[a] <- FALSE
        break
      }
    }
  }
  m[, keep, drop = FALSE]
}

# Compact letter display by the insert-and-absorb algorithm of Piepho (2004),
# the algorithm of multcompView::multcompLetters. Each pair that differs
# (p < alpha) splits every letter column that holds both groups into two, one
# without each group; then the redundant columns are absorbed. The pairs are
# taken in the order of the levels, the columns are sorted so that "a" starts
# at the first group, and the letters run a-z then A-Z. Groups that share a
# letter do not differ. The Python engine has the same steps (util.cld_letters).
cld_letters <- function(levels, group1, group2, p, alpha = 0.05) {
  levels <- as.character(levels)
  n <- length(levels)
  m <- matrix(TRUE, nrow = n, ncol = 1)
  i1 <- match(as.character(group1), levels)
  i2 <- match(as.character(group2), levels)
  sig <- which(!is.na(p) & p < alpha & !is.na(i1) & !is.na(i2))
  lo <- pmin(i1[sig], i2[sig]); hi <- pmax(i1[sig], i2[sig])
  for (k in order(lo, hi)) {
    i <- lo[k]; j <- hi[k]
    both <- which(m[i, ] & m[j, ])
    if (length(both) == 0) next
    dup <- m[, both, drop = FALSE]
    dup[i, ] <- FALSE
    m[j, both] <- FALSE
    m <- .cld_absorb(cbind(m, dup))
  }
  keys <- apply(m, 2, .cld_key)
  m <- m[, order(keys, method = "radix"), drop = FALSE]
  pool <- c(letters, LETTERS)
  if (ncol(m) > length(pool)) pool <- c(pool, sprintf("z%d", seq_len(ncol(m) - length(pool))))
  out <- vapply(seq_len(n), .cld_row, character(1), m = m, pool = pool)
  stats::setNames(out, levels)
}

.cld_row <- function(i, m, pool) paste(pool[which(m[i, ])], collapse = "")

cap_first <- function(s) {
  if (length(s) == 0 || is.na(s) || nchar(s) == 0) return(s)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}
