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

cap_first <- function(s) {
  if (length(s) == 0 || is.na(s) || nchar(s) == 0) return(s)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}
