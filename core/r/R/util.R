# Small shared helpers.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Condition handlers for tryCatch, at module level so no recipe defines a
# function inside another one.
na_on_error <- function(e) NA_real_
null_on_error <- function(e) NULL

# Whether one shapiro_safe() result says the values look normal.
is_normal_check <- function(n) isTRUE(n$normal)

# The installed version of one package, as a string.
pkg_version_string <- function(p) as.character(utils::packageVersion(p))

# Shapiro-Wilk normality on a numeric vector, guarded for small and large n. A
# sample with zero range (every value the same) is not normal and has no P, in
# both engines (util.shapiro_normal).
shapiro_safe <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 3)   return(list(test = "shapiro-wilk", n = n, p = NA_real_, normal = NA))
  if (max(v) == min(v)) return(list(test = "shapiro-wilk", n = n, p = NA_real_, normal = FALSE))
  if (n > 5000) return(list(test = "shapiro-wilk", n = n, p = NA_real_, normal = TRUE))
  p <- tryCatch(stats::shapiro.test(v)$p.value, error = na_on_error)
  list(test = "shapiro-wilk", n = n, p = p, normal = isTRUE(p > 0.05))
}

# Whether a sample has zero variance: two or more finite values, all equal.
is_zero_variance <- function(v) {
  v <- v[is.finite(v)]
  length(v) >= 2 && max(v) == min(v)
}

# The names of the groups of a named list that have zero variance.
zero_variance_groups <- function(groups) {
  names(groups)[vapply(groups, is_zero_variance, logical(1))]
}

# Stop on a --test value that a recipe does not know, and list the values that
# it knows. The Python engine writes the same message (util.check_test_method).
check_test_method <- function(recipe, method, allowed) {
  m <- tolower(method %||% "auto")
  if (nzchar(m) && !m %in% allowed) {
    stop(sprintf("%s does not know --test '%s'; use one of: %s", recipe, m,
                 paste(allowed, collapse = ", ")))
  }
  invisible(m)
}

# The quoted names of a list, joined as text: 'a'; 'a' and 'b'; 'a', 'b' and 'c'.
quote_join <- function(x) {
  q <- sprintf("'%s'", x)
  if (length(q) <= 1) return(paste(q, collapse = ""))
  paste(paste(q[-length(q)], collapse = ", "), "and", q[length(q)])
}

# The message that stops a test that divides by a variance of zero. `subject`
# names the samples, for example "the group 'a' has" or "the within-pair
# differences have"; `hint` names what to run instead. The Python engine writes
# the same words (util.zero_variance_message).
zero_variance_message <- function(test, subject, hint) {
  sprintf("%s cannot run because %s zero variance (all values are equal); %s",
          test, subject, hint)
}

# The subject of zero_variance_message() for a set of groups.
zero_variance_subject <- function(groups) {
  if (length(groups) == 1) sprintf("the group %s has", quote_join(groups))
  else sprintf("the groups %s have", quote_join(groups))
}

# Rank-biserial correlation for two independent groups, from the Mann-Whitney U
# of the first group. Negative when the first group ranks lower. This is
# computed here rather than with rstatix::wilcox_effsize, which needs the coin
# package and returns the magnitude without the sign. The Python engine uses the
# same formula.
rank_biserial_unpaired <- function(u, n1, n2) {
  2 * as.numeric(u) / (n1 * n2) - 1
}

# Signed Cohen's d of two groups, first minus second. Two independent groups use
# the pooled standard deviation; paired data use dz, the mean of the
# differences over their standard deviation. The Python engine uses the same
# formulas (util.cohens_d_r).
cohens_d_signed <- function(a, b, paired = FALSE) {
  if (paired) {
    d <- a - b
    return(mean(d) / stats::sd(d))
  }
  na <- length(a); nb <- length(b)
  pooled <- sqrt(((na - 1) * stats::var(a) + (nb - 1) * stats::var(b)) / (na + nb - 2))
  (mean(a) - mean(b)) / pooled
}

# Rank-biserial correlation for paired differences (Kerby's formula). Zero
# differences are dropped, as in the signed-rank test itself.
rank_biserial_paired <- function(d) {
  d <- d[is.finite(d) & d != 0]
  if (length(d) == 0) return(NA_real_)
  rk <- rank(abs(d))
  (sum(rk[d > 0]) - sum(rk[d < 0])) / sum(rk)
}

# The smallest P that a rank test can reach at the given group sizes: the same
# test function on data with no tie and complete separation. Independent groups
# take the ranks 1 to N in group order ("mann_whitney" for two groups, "kruskal"
# for the Kruskal-Wallis test); paired data take n distinct positive
# differences ("signed_rank"). NA when a group is empty. The Python engine has
# the same rule (util.rank_test_min_p).
rank_test_min_p <- function(sizes, test) {
  sizes <- as.integer(sizes)
  if (length(sizes) == 0 || any(is.na(sizes)) || any(sizes < 1)) return(NA_real_)
  if (test == "signed_rank") {
    return(suppressWarnings(stats::wilcox.test(as.numeric(seq_len(sizes[1])))$p.value))
  }
  if (length(sizes) < 2) return(NA_real_)
  v <- as.numeric(seq_len(sum(sizes)))
  g <- rep(seq_along(sizes), sizes)
  if (test == "mann_whitney") {
    return(suppressWarnings(stats::wilcox.test(v[g == 1], v[g == 2])$p.value))
  }
  stats::kruskal.test(v, factor(g))$p.value
}

# Whether the rank test cannot reach P < 0.05: then --test auto keeps the
# parametric test after a failed normality check.
rank_test_cannot_reach <- function(min_p) isTRUE(min_p > 0.05)

# The sentence of a methods text on the attainable P of the rank test, empty
# when the rank test can reach P < 0.05. `kept` says that the parametric test
# was kept for this reason. The Python engine writes the same words.
rank_floor_text <- function(min_p, kept) {
  if (!rank_test_cannot_reach(min_p)) return("")
  txt <- sprintf("The rank test cannot reach P < 0.05 at these group sizes (smallest attainable P = %s)",
                 sig_text(min_p, 2))
  if (isTRUE(kept)) paste0(txt, ", so the parametric test was kept.") else paste0(txt, ".")
}

# Resolve which two-group test to run from the request and the assumption checks.
# The t-test of two independent groups is Welch's t-test, which does not assume
# equal variances; --test student_t forces Student's t-test. Levene's test
# (equal_var) still runs for the stats table, but it does not choose the test.
# When the normality check fails but the rank test cannot reach P < 0.05 at the
# group sizes (rank_min_p), --test auto keeps the parametric test, unless the
# t-test cannot run (t_can_run, a variance of zero).
resolve_two_group_method <- function(method, all_normal, equal_var, paired,
                                     rank_min_p = NA_real_, t_can_run = TRUE) {
  method <- tolower(method %||% "auto")
  kept <- FALSE
  if (method == "auto") {
    if (isTRUE(all_normal)) {
      family <- "t"; var_equal <- FALSE
    } else if (rank_test_cannot_reach(rank_min_p) && isTRUE(t_can_run)) {
      family <- "t"; var_equal <- FALSE; kept <- TRUE
    } else {
      family <- "wilcoxon"; var_equal <- NA
    }
  } else if (method %in% c("t", "t_test", "student_t", "welch_t")) {
    family <- "t"
    var_equal <- method == "student_t"
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
  list(family = family, var_equal = var_equal, label = label, paired = paired,
       kept_parametric = kept, rank_min_p = rank_min_p)
}

# Installed versions of the named packages, as a named list of strings.
pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, pkg_version_string)
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

# The P value as a phrase for a methods sentence: "P < 0.001" below 0.001, else
# "P = " and two significant digits. A table cell keeps the bare value.
p_phrase <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) "P < 0.001" else paste("P =", trimws(formatC(p, format = "g", digits = 2)))
}

# A number with a fixed count of significant digits for a figure label or a
# methods sentence. Unlike "%.3g" it keeps the trailing zeros (3.70, not 3.7),
# and it drops the trailing point that "%#g" leaves on a whole number.
sig_text <- function(x, digits = 3) {
  sub("\\.$", "", sprintf(paste0("%#.", digits, "g"), x))
}

# A degree of freedom for a methods sentence: a whole number as an integer, and
# a fractional df (Welch, Greenhouse-Geisser) with two decimals.
df_text <- function(d) {
  if (abs(d - round(d)) < 1e-9) sprintf("%d", as.integer(round(d))) else sprintf("%.2f", d)
}

# Text labels for the values of a grouping column. A whole number prints with no
# decimal point and no exponent, so a dose of 80 reads "80" in both engines, also
# when the column holds a missing value. Other numbers print with 15
# significant digits. The Python engine has the same rule.
level_labels <- function(x) {
  if (!is.numeric(x)) return(as.character(x))
  whole <- is.finite(x) & x == round(x) & abs(x) < 1e15
  out <- trimws(formatC(x, digits = 15, format = "g"))
  out[whole] <- trimws(formatC(x[whole], format = "f", digits = 0))
  out[is.na(x)] <- NA_character_
  out
}

# The order of the levels of a grouping column, one rule in both engines. With
# --levels the order is the listed order: a value of the data that the list does
# not name stops the run, and a listed value that no row holds is left out.
# Without --levels the levels sort as numbers when every label is a number, else
# in code-point order, which does not depend on the locale of the image.
level_order <- function(labels, levels = NULL, what = "group") {
  labels <- if (is.numeric(labels)) level_labels(labels) else as.character(labels)
  u <- unique(labels[!is.na(labels)])
  if (length(levels) > 0) {
    levels <- unique(trimws(as.character(levels)))
    miss <- setdiff(u, levels)
    if (length(miss) > 0) {
      stop(sprintf("--levels does not list the %s(s) %s of the data", what,
                   paste(sprintf("'%s'", miss), collapse = ", ")))
    }
    return(levels[levels %in% u])
  }
  num <- suppressWarnings(as.numeric(u))
  if (length(u) > 0 && all(is.finite(num) | is.infinite(num))) {
    return(u[order(num, u, method = "radix")])
  }
  sort(u, method = "radix")
}
