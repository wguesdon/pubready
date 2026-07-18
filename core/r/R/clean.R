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
