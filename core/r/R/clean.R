# Minimal, logged cleaning: only what is needed to run the analysis, with every
# step recorded so the raw-to-analyzed path stays auditable.

clean_tidy <- function(df, spec) {
  steps <- character(0)
  n0 <- nrow(df)
  x <- spec$data$x
  y <- spec$data$y
  if (!x %in% names(df)) stop(sprintf("column '%s' not found in data", x))
  if (!y %in% names(df)) stop(sprintf("column '%s' not found in data", y))

  yv <- suppressWarnings(as.numeric(df[[y]]))
  coerced <- sum(is.na(yv) & !is.na(df[[y]]))
  if (coerced > 0) {
    steps <- c(steps, sprintf("Coerced column '%s' to numeric; %d value(s) became NA.", y, coerced))
  }
  df[[y]] <- yv

  keep <- !is.na(df[[x]]) & is.finite(df[[y]])
  dropped <- sum(!keep)
  if (dropped > 0) {
    steps <- c(steps, sprintf("Dropped %d row(s) with missing '%s' or '%s' (%d -> %d rows).",
                              dropped, x, y, n0, sum(keep)))
  }
  df <- df[keep, , drop = FALSE]

  if (length(steps) == 0) steps <- "No cleaning applied; input used as-is."
  list(df = df, steps = steps)
}

# Classify and coerce the numeric columns of a table, one rule in both engines.
# A column is numeric when it is already numeric, or when 80% or more of its
# cells that are not missing read as finite numbers. A cell of such a column
# that does not read as a number, for example "<0.1" or "45 y", becomes NA, and
# a step records the count and the first such value. Every other column is text.
# `cols` limits the check to the named columns. Returns the table, the names of
# the numeric and the text columns, and the steps.
coerce_numeric_columns <- function(df, cols = names(df), threshold = 0.8) {
  steps <- character(0)
  numeric <- character(0)
  text <- character(0)
  for (cn in cols) {
    v <- df[[cn]]
    if (is.numeric(v)) {
      numeric <- c(numeric, cn)
      next
    }
    present <- !is.na(v) & trimws(as.character(v)) != ""
    num <- suppressWarnings(as.numeric(as.character(v)))
    ok <- present & is.finite(num)
    if (sum(present) > 0 && sum(ok) / sum(present) >= threshold) {
      bad <- present & !ok
      if (any(bad)) {
        steps <- c(steps, sprintf(
          "Coerced column '%s' to numeric; %d value(s) became NA, for example '%s'.",
          cn, sum(bad), as.character(v[bad][1])))
      }
      df[[cn]] <- ifelse(ok, num, NA_real_)
      numeric <- c(numeric, cn)
    } else {
      text <- c(text, cn)
    }
  }
  list(df = df, numeric = numeric, text = text, steps = steps)
}
