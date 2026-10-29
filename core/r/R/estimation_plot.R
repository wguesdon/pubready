# Recipe: estimation plot in the style of DABEST (Ho et al. 2019). Each group is
# compared with one reference group, and the figure shows the effect size, the
# mean difference, with its bootstrap confidence interval next to the raw data.
#
# Two groups give the Gardner-Altman layout: the raw data on the left, and on the
# right a difference axis whose zero sits at the reference mean, carrying the
# bootstrap distribution of the mean difference, the estimate and the 95% CI.
# Three or more groups give the Cumming layout: the raw data on top, and a lower
# panel with each difference against the reference. --paired draws a slopegraph
# of each subject and resamples whole pairs.
#
# The bootstrap must give identical resamples in the R and the Python engine, so
# neither engine uses its own random number generator. Both run the Park-Miller
# minimal standard generator in exact double precision:
#   state <- (48271 * state) mod (2^31 - 1),  u = state / (2^31 - 1),
#   index = floor(u * n) (+1 in R), seed 20261029, restarted for each comparison.
# Each resample b draws the n_ref reference indices first and then the n_group
# indices of the compared group (paired: n pair indices). The 95% interval is the
# bias-corrected and accelerated (BCa) interval: z0 = qnorm(the proportion of
# bootstrap differences strictly below the estimate), the acceleration comes from
# the jackknife of the mean difference over every observation of both groups
# (paired: every pair), and the ends are type-7 quantiles of the 5000 bootstrap
# differences.

.est_seed <- 20261029
.est_modulus <- 2147483647
.est_multiplier <- 48271
.est_n_boot <- 5000
.est_block <- 4096

# (z * a) mod m for 0 <= z, a < 2^31, exact in double precision: a is split into
# 16-bit halves so that no product exceeds 2^47. Vectorised over z and a.
.est_mulmod <- function(z, a) {
  a_hi <- floor(a / 65536)
  a_lo <- a - a_hi * 65536
  ((((z * a_hi) %% .est_modulus) * 65536) %% .est_modulus +
     (z * a_lo) %% .est_modulus) %% .est_modulus
}

# The multipliers 48271^j mod m for j = 1..b. One block of the generator is then
# state * multipliers mod m, which gives the next b states at once.
.est_block_multipliers <- function(b) {
  out <- numeric(b)
  out[1] <- .est_multiplier
  for (j in seq_len(b - 1)) out[j + 1] <- .est_mulmod(out[j], .est_multiplier)
  out
}

# n uniforms of the Park-Miller generator that continue from state.
.est_uniforms <- function(n, state, mult) {
  b <- length(mult)
  out <- numeric(n)
  pos <- 0
  while (pos < n) {
    blk <- .est_mulmod(rep(state, b), mult)
    k <- min(b, n - pos)
    out[pos + seq_len(k)] <- blk[seq_len(k)]
    state <- blk[k]
    pos <- pos + k
  }
  list(u = out / .est_modulus, state = state)
}

# The bootstrap distribution of the mean difference (group minus reference).
# Unpaired: r and g are the two samples. Paired: r and g are matched by position.
.est_boot_diffs <- function(r, g, paired, n_boot = .est_n_boot) {
  mult <- .est_block_multipliers(.est_block)
  state <- .est_seed
  n_r <- length(r); n_g <- length(g)
  per <- if (paired) n_r else n_r + n_g
  chunk <- max(1L, as.integer(floor(2e6 / per)))
  out <- numeric(n_boot)
  done <- 0
  while (done < n_boot) {
    m <- min(chunk, n_boot - done)
    gen <- .est_uniforms(m * per, state, mult)
    state <- gen$state
    u <- matrix(gen$u, nrow = per)
    if (paired) {
      idx <- floor(u * n_r) + 1
      diffs <- colSums(matrix(g[idx] - r[idx], nrow = n_r)) / n_r
    } else {
      ir <- floor(u[seq_len(n_r), , drop = FALSE] * n_r) + 1
      ig <- floor(u[n_r + seq_len(n_g), , drop = FALSE] * n_g) + 1
      diffs <- colSums(matrix(g[ig], nrow = n_g)) / n_g -
        colSums(matrix(r[ir], nrow = n_r)) / n_r
    }
    out[done + seq_len(m)] <- diffs
    done <- done + m
  }
  out
}

# Jackknife values of the mean difference: leave one observation out of either
# group (paired: leave one pair out).
.est_jackknife <- function(r, g, paired) {
  if (paired) {
    d <- g - r
    return((sum(d) - d) / (length(d) - 1))
  }
  mr <- sum(r) / length(r); mg <- sum(g) / length(g)
  c(mg - (sum(r) - r) / (length(r) - 1), (sum(g) - g) / (length(g) - 1) - mr)
}

# The BCa 95% interval from the bootstrap differences and the jackknife values.
# A bootstrap difference counts as below the estimate only when it is below by
# more than 1e-11 of the data scale, so that rounding in the last bit of a sum
# cannot decide a tie that is exact in real arithmetic.
.est_bca <- function(boot, est, jack, scale) {
  n_boot <- length(boot)
  if (stats::sd(boot) == 0) return(c(est, est))
  below <- sum(boot < est - 1e-11 * scale) / n_boot
  below <- min(max(below, 0.5 / n_boot), 1 - 0.5 / n_boot)
  z0 <- stats::qnorm(below)
  dj <- mean(jack) - jack
  den <- 6 * sum(dj^2)^1.5
  acc <- if (den > 0) sum(dj^3) / den else 0
  zq <- stats::qnorm(c(0.025, 0.975))
  alpha <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
  unname(stats::quantile(boot, alpha, type = 7))
}

# Hedges' g as DABEST computes it: Cohen's d over the pooled SD (paired: over the
# average SD of the two conditions), times the exact small-sample correction J for
# n1 + n2 - 2 degrees of freedom.
.est_hedges_g <- function(r, g, paired) {
  n1 <- length(r); n2 <- length(g)
  v1 <- stats::var(r); v2 <- stats::var(g)
  s <- if (paired) sqrt((v1 + v2) / 2) else sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  dfj <- n1 + n2 - 2
  j <- exp(lgamma(dfj / 2) - lgamma((dfj - 1) / 2)) / sqrt(dfj / 2)
  (mean(g) - mean(r)) / s * j
}

# One comparison: the group against the reference.
.est_compare <- function(r, g, paired, ref_name, grp_name) {
  est <- sum(g) / length(g) - sum(r) / length(r)
  boot <- .est_boot_diffs(r, g, paired)
  ci <- .est_bca(boot, est, .est_jackknife(r, g, paired), max(abs(c(r, g))))
  tt <- stats::t.test(g, r, paired = paired, var.equal = FALSE)
  row <- data.frame(
    reference = ref_name, group = grp_name, n_ref = length(r), n_group = length(g),
    mean_ref = mean(r), mean_group = mean(g), mean_diff = est,
    ci_low = ci[1], ci_high = ci[2], hedges_g = .est_hedges_g(r, g, paired),
    test = if (paired) "paired t-test" else "Welch two-sample t-test",
    statistic = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
    stringsAsFactors = FALSE)
  list(row = row, boot = boot)
}

# Swarm offsets: the values are cut into bins of height bin_h, and the points of
# one bin spread out from the centre (0, +1, -1, +2, ...) in value order, by step
# or less, so that no bin is wider than 2 * half.
.est_swarm <- function(v, vmin, bin_h, step = 0.07, half = 0.3) {
  off <- numeric(length(v))
  if (length(v) == 0) return(off)
  bins <- floor((v - vmin) / bin_h)
  for (b in unique(bins)) {
    sel <- which(bins == b)
    sel <- sel[order(v[sel], method = "radix")]
    k <- length(sel)
    if (k < 2) next
    slots <- ceiling((seq_len(k) - 1) / 2) * ifelse(seq_len(k) %% 2 == 0, 1, -1)
    slots <- slots - mean(slots)
    off[sel] <- slots * min(step, 2 * half / (k - 1))
  }
  off
}

# A Gaussian kernel density of the bootstrap differences on 128 points between
# their extremes, with the nrd0 bandwidth, scaled to a maximum of 1.
.est_density <- function(boot) {
  grid <- seq(min(boot), max(boot), length.out = 128)
  q <- stats::quantile(boot, c(0.25, 0.75), type = 7)
  s <- min(stats::sd(boot), unname(q[2] - q[1]) / 1.34)
  if (!is.finite(s) || s <= 0) s <- stats::sd(boot)
  bw <- 0.9 * s * length(boot)^(-0.2)
  dens <- colMeans(stats::dnorm(outer(boot, grid, "-") / bw))
  data.frame(y = grid, w = dens / max(dens))
}

# The half violin of one comparison, placed at x = pos and shifted by offset in y.
.est_violin_poly <- function(i, comps, pos, offset) {
  dd <- .est_density(comps[[i]]$boot)
  data.frame(cmp = i, x = c(pos[i], pos[i] + 0.38 * dd$w, pos[i]),
             y = offset[i] + c(dd$y[1], dd$y, dd$y[nrow(dd)]))
}

# The width of the Gardner-Altman figure, grown for long column labels: each of
# the three columns gets 1.1 in, or 0.1 in per character of the longest label
# line, plus 1.6 in for the two axes. Python uses the same rule (ga_width).
.est_ga_width <- function(labels) {
  chars <- max(nchar(unlist(strsplit(as.character(labels), "\n", fixed = TRUE))))
  max(4.9, 3 * max(1.1, 0.1 * chars) + 1.6)
}

# Wrap a category label into lines shorter than width, as scales::label_wrap
# does, but also break after an underscore, so that a name such as
# "stimulated_plus_inhibitor" wraps too. The Python engine has the same rule.
.est_wrap_one <- function(s, width) {
  line <- ""
  out <- character(0)
  for (w in strsplit(s, " ", fixed = TRUE)[[1]]) {
    if (!nzchar(w)) next
    parts <- strsplit(gsub("_", "_\001", w, fixed = TRUE), "\001", fixed = TRUE)[[1]]
    for (j in seq_along(parts)) {
      sep <- if (j == 1 && nzchar(line)) " " else ""
      if (nzchar(line) && nchar(line) + nchar(sep) + nchar(parts[j]) >= width) {
        out <- c(out, line)
        line <- parts[j]
      } else {
        line <- paste0(line, sep, parts[j])
      }
    }
  }
  paste(c(out, line), collapse = "\n")
}

.est_wrap <- function(s, width = 14) {
  vapply(as.character(s), .est_wrap_one, character(1), width = width, USE.NAMES = FALSE)
}

# The point size of the raw data, smaller for a large group (ggplot2 size in mm).
.est_point_size <- function(d) {
  n <- max(table(d$group))
  if (n > 100) 1.2 else if (n > 40) 1.6 else 2
}

# The raw-data layer data: one row per value with its swarm offset.
.est_points <- function(d, lv, pos, bin_h, vmin) {
  out <- d
  out$pos <- pos[match(d$group, lv)]
  out$off <- 0
  for (g in lv) {
    sel <- d$group == g
    out$off[sel] <- .est_swarm(d$y[sel], vmin, bin_h)
  }
  out
}

# Mean +/- SD of each group as a line with a gap at the mean (two segments).
.est_gapped <- function(d, lv, xpos, gap) {
  m <- vapply(lv, .est_group_mean, numeric(1), d = d)
  s <- vapply(lv, .est_group_sd, numeric(1), d = d)
  lo <- data.frame(group = lv, x = xpos, y = m - s, yend = pmin(m - gap, m))
  hi <- data.frame(group = lv, x = xpos, y = pmax(m + gap, m), yend = m + s)
  keep_lo <- lo$yend > lo$y; keep_hi <- hi$yend > hi$y
  rbind(lo[keep_lo, , drop = FALSE], hi[keep_hi, , drop = FALSE])
}

# The x of each mean +/- SD line: just right of the widest point of the swarm.
.est_gap_x <- function(pts, lv, pos) {
  pos + vapply(lv, .est_group_reach, numeric(1), pts = pts) + 0.15
}

.est_group_reach <- function(g, pts) max(abs(pts$off[pts$group == g]), 0)

.est_group_mean <- function(g, d) mean(d$y[d$group == g])
.est_group_sd <- function(g, d) {
  v <- d$y[d$group == g]
  if (length(v) > 1) stats::sd(v) else 0
}

.est_group_n <- function(g, d) sum(d$group == g)

# Subject lines of a paired design, between consecutive groups in display order.
.est_slopes <- function(d, lv, pos) {
  out <- NULL
  for (i in seq_len(length(lv) - 1)) {
    a <- d[d$group == lv[i], c("id", "y")]
    b <- d[d$group == lv[i + 1], c("id", "y")]
    m <- merge(a, b, by = "id", suffixes = c("_a", "_b"))
    if (nrow(m)) {
      out <- rbind(out, data.frame(x = pos[i], xend = pos[i + 1], y = m$y_a, yend = m$y_b))
    }
  }
  out
}

.est_fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.2g", p)

recipe_estimation_plot <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y
  if (is.null(x) || is.null(y)) {
    stop("estimation_plot needs --x (the group column) and --y (the value column)")
  }
  paired <- isTRUE(spec$test$paired)
  id <- spec$data$id
  if (paired && (is.null(id) || !nzchar(id))) {
    stop("estimation_plot --paired needs --id (the subject column)")
  }
  method <- tolower(spec$test$method %||% "auto")
  if (!method %in% c("auto", "welch_t", "t")) {
    stop("estimation_plot reports a Welch t-test (a paired t-test with --paired); use --test auto")
  }
  if (paired && !id %in% names(df)) stop(sprintf("column '%s' not found in data", id))

  cl <- clean_tidy(df, spec)
  steps <- cl$steps
  steps <- steps[steps != "No cleaning applied; input used as-is."]
  d0 <- cl$df
  d <- data.frame(group = as.character(d0[[x]]), y = as.numeric(d0[[y]]),
                  id = if (paired) as.character(d0[[id]]) else NA_character_,
                  stringsAsFactors = FALSE)
  if (paired) {
    miss <- is.na(d$id) | d$id == ""
    if (any(miss)) {
      steps <- c(steps, sprintf("Dropped %d row(s) with a missing '%s'.", sum(miss), id))
      d <- d[!miss, , drop = FALSE]
    }
    if (anyDuplicated(paste(d$id, d$group, sep = "\r"))) {
      stop(sprintf("estimation_plot --paired needs one value per subject and group; '%s' repeats a subject within a group", id))
    }
  }

  lv_sorted <- sort(unique(d$group), method = "radix")
  ref <- spec$data$ref %||% lv_sorted[1]
  if (!ref %in% lv_sorted) {
    stop(sprintf("--ref '%s' is not a level of '%s' (levels: %s)", ref, x,
                 paste(lv_sorted, collapse = ", ")))
  }
  lv <- c(ref, setdiff(lv_sorted, ref))
  if (length(lv) < 2) stop(sprintf("estimation_plot needs at least 2 groups in '%s'", x))

  comps <- list()
  for (g in lv[-1]) {
    if (paired) {
      a <- d[d$group == ref, c("id", "y")]
      b <- d[d$group == g, c("id", "y")]
      m <- merge(a, b, by = "id", suffixes = c("_r", "_g"), sort = FALSE)
      m <- m[order(m$id, method = "radix"), , drop = FALSE]
      lost <- length(unique(c(a$id, b$id))) - nrow(m)
      if (lost > 0) {
        steps <- c(steps, sprintf("Comparison %s against %s: left out %d subject(s) without a value in both groups.",
                                  g, ref, lost))
      }
      if (nrow(m) < 3) {
        stop(sprintf("estimation_plot --paired needs at least 3 subjects with a value in both '%s' and '%s'", ref, g))
      }
      r_v <- m$y_r; g_v <- m$y_g
    } else {
      r_v <- d$y[d$group == ref]; g_v <- d$y[d$group == g]
      if (length(r_v) < 2 || length(g_v) < 2) {
        stop(sprintf("estimation_plot needs at least 2 values in each group; '%s' has %d and '%s' has %d",
                     ref, length(r_v), g, length(g_v)))
      }
    }
    comps[[g]] <- .est_compare(r_v, g_v, paired, ref, g)
  }
  if (length(steps) == 0) steps <- "No cleaning applied; input used as-is."
  stats_df <- do.call(rbind, lapply(comps, `[[`, "row"))
  rownames(stats_df) <- NULL

  ap <- spec$appearance
  pal <- stats::setNames(pub_palette(length(lv), ap$palette), lv)
  y_lab <- ap$y_label %||% y
  fig <- if (length(lv) == 2) {
    .est_gardner_altman(d, lv, comps, stats_df, pal, paired, ap, x, y_lab)
  } else {
    .est_cumming(d, lv, comps, stats_df, pal, paired, ap, x, y_lab)
  }

  test_name <- if (paired) "paired mean difference with a BCa bootstrap CI" else
    "mean difference with a BCa bootstrap CI"
  test_meta <- list(
    name = test_name, reference = ref, paired = paired,
    bootstrap = list(resamples = .est_n_boot, generator = "Park-Miller minimal standard (48271, 2^31 - 1)",
                     seed = .est_seed, interval = "bias-corrected and accelerated (BCa), 95%"),
    test = stats_df$test[1],
    comparisons = stats::setNames(lapply(seq_len(nrow(stats_df)), .est_meta_row, rows = stats_df),
                                  stats_df$group))
  df_used <- d[d$group %in% lv, c("group", "y", if (paired) "id"), drop = FALSE]
  names(df_used) <- c(x, y, if (paired) id)
  qc <- if (paired) {
    list(quantity = "paired differences",
         panels = stats::setNames(lapply(lv[-1], .est_paired_diffs, d = d, ref = ref),
                                  paste(lv[-1], "-", ref)))
  } else {
    list(quantity = "group values", panels = split(d$y, factor(d$group, levels = lv)))
  }
  list(plot = fig$plot, stats = stats_df, test_meta = test_meta,
       resolved = list(family = "estimation", label = test_name),
       methods = .est_methods(stats_df, ref, paired, spec, y),
       build_script = list(fn = .est_emit_script, args = list(spec = spec, lv = lv, paired = paired)),
       label = slugify(y), df_used = df_used, clean_steps = steps, qc = qc,
       width = fig$width, height = fig$height)
}

.est_paired_diffs <- function(g, d, ref) {
  m <- merge(d[d$group == ref, c("id", "y")], d[d$group == g, c("id", "y")], by = "id")
  m$y.y - m$y.x
}

.est_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(mean_diff = r$mean_diff, ci_low = r$ci_low, ci_high = r$ci_high,
       hedges_g = r$hedges_g, statistic = r$statistic, df = r$df, p_value = r$p_value)
}

.est_axis_labels <- function(lv, d) {
  sprintf("%s\n(n = %d)", .est_wrap(lv), vapply(lv, .est_group_n, integer(1), d = d))
}

# Two groups: the Gardner-Altman layout in one panel. The difference sits at
# x = 3 and is drawn at y = mean_ref + difference; the secondary axis on the right
# reads it back as a difference, so its zero sits at the reference mean.
.est_gardner_altman <- function(d, lv, comps, rows, pal, paired, ap, x, y_lab) {
  pos <- c(1, 2)
  m_ref <- rows$mean_ref[1]
  m_grp <- rows$mean_group[1]
  yr <- range(c(d$y, m_ref + comps[[1]]$boot))
  bin_h <- diff(range(d$y)) / 40
  if (!is.finite(bin_h) || bin_h <= 0) bin_h <- 1
  pts <- .est_points(d, lv, pos, bin_h, min(d$y))
  gap <- diff(yr) * 0.012
  gx <- if (paired) c(1 - 0.16, 2 + 0.16) else .est_gap_x(pts, lv, pos)
  gl <- .est_gapped(d, lv, gx, gap)
  vio <- .est_violin_poly(1, comps, 3, m_ref)
  guide <- data.frame(x = 2.55, xend = 3.5, y = c(m_ref, m_grp))
  ci <- data.frame(x = 3, y = m_ref + rows$ci_low[1], yend = m_ref + rows$ci_high[1],
                   est = m_ref + rows$mean_diff[1])
  labs_x <- c(.est_axis_labels(lv, d), .est_wrap(sprintf("%s minus %s", lv[2], lv[1])))
  pts$group <- factor(pts$group, levels = lv)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = guide, ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          colour = "grey60", linewidth = 0.3, linetype = "dashed")
  if (paired) {
    sl <- .est_slopes(d, lv, pos)
    p <- p + ggplot2::geom_segment(data = sl, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   colour = "grey55", linewidth = 0.35, alpha = 0.7)
    pts$off <- 0
  }
  p <- p +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = pos + off, y = y, fill = group),
                        shape = 21, size = .est_point_size(d), stroke = 0.25, alpha = 0.85) +
    ggplot2::geom_segment(data = gl, ggplot2::aes(x = x, xend = x, y = y, yend = yend),
                          colour = "black", linewidth = 0.8) +
    ggplot2::geom_polygon(data = vio, ggplot2::aes(x = x, y = y), fill = pal[[2]],
                          alpha = 0.45, colour = NA) +
    ggplot2::geom_segment(data = ci, ggplot2::aes(x = x, xend = x, y = y, yend = yend),
                          colour = "black", linewidth = 0.9) +
    ggplot2::geom_point(data = ci, ggplot2::aes(x = x, y = est), colour = "black", size = 2.6) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(1, 2, 3), labels = labs_x, limits = c(0.55, 3.55),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . - m_ref, name = "Mean difference")) +
    ggplot2::labs(x = ap$x_label, y = y_lab, title = ap$title,
                  subtitle = sprintf("Mean difference %s (95%% CI %s to %s)\n%s, %s",
                                     sprintf("%.3g", rows$mean_diff[1]), sprintf("%.3g", rows$ci_low[1]),
                                     sprintf("%.3g", rows$ci_high[1]), rows$test[1],
                                     .est_fmt_p(rows$p_value[1]))) +
    pub_theme(ap$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   axis.line.y.right = ggplot2::element_line(colour = "black", linewidth = 0.5))
  list(plot = p, width = .est_ga_width(labs_x), height = 4.4)
}

# Three or more groups: the Cumming layout, raw data on top and the differences
# against the reference below, on the same x positions.
.est_cumming <- function(d, lv, comps, rows, pal, paired, ap, x, y_lab) {
  k <- length(lv)
  pos <- seq_len(k)
  bin_h <- diff(range(d$y)) / 40
  if (!is.finite(bin_h) || bin_h <= 0) bin_h <- 1
  pts <- .est_points(d, lv, pos, bin_h, min(d$y))
  gap <- diff(range(d$y)) * 0.012
  gx <- if (paired) pos + 0.16 else .est_gap_x(pts, lv, pos)
  gl <- .est_gapped(d, lv, gx, gap)
  # In a paired design the subject lines leave each group on both sides, so the
  # mean +/- SD lines would cross them; the slopegraph shows the spread instead.
  if (paired) gl <- gl[0, , drop = FALSE]
  pts$group <- factor(pts$group, levels = lv)
  xlim <- c(0.5, k + 0.6)

  top <- ggplot2::ggplot()
  if (paired) {
    sl <- .est_slopes(d, lv, pos)
    if (!is.null(sl)) {
      top <- top + ggplot2::geom_segment(data = sl, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                         colour = "grey55", linewidth = 0.3, alpha = 0.6)
    }
    pts$off <- 0
  }
  top <- top +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = pos + off, y = y, fill = group),
                        shape = 21, size = .est_point_size(d), stroke = 0.25, alpha = 0.85) +
    ggplot2::geom_segment(data = gl, ggplot2::aes(x = x, xend = x, y = y, yend = yend),
                          colour = "black", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = pos, labels = .est_axis_labels(lv, d), limits = xlim,
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = y_lab, title = ap$title,
                  subtitle = sprintf("Mean difference against %s with its 95%% BCa bootstrap CI",
                                     lv[1])) +
    pub_theme(ap$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  vio <- do.call(rbind, lapply(seq_along(comps), .est_violin_poly, comps = comps, pos = pos[-1],
                               offset = rep(0, length(comps))))
  vio$group <- factor(lv[-1][vio$cmp], levels = lv)
  ci <- data.frame(x = pos[-1], y = rows$ci_low, yend = rows$ci_high, est = rows$mean_diff)
  span <- range(c(0, vio$y, ci$y, ci$yend))
  plab <- data.frame(x = pos[-1], y = span[2] + diff(span) * 0.1,
                     label = vapply(rows$p_value, .est_fmt_p, character(1)))
  bottom <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_polygon(data = vio, ggplot2::aes(x = x, y = y, group = cmp, fill = group),
                          alpha = 0.45, colour = NA) +
    ggplot2::geom_segment(data = ci, ggplot2::aes(x = x, xend = x, y = y, yend = yend),
                          colour = "black", linewidth = 0.9) +
    ggplot2::geom_point(data = ci, ggplot2::aes(x = x, y = est), colour = "black", size = 2.4) +
    ggplot2::geom_text(data = plab, ggplot2::aes(x = x, y = y, label = label), size = 3.2,
                       colour = "grey30") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = pos, labels = NULL, limits = xlim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.1))) +
    ggplot2::labs(x = ap$x_label, y = paste0("Mean difference\n", .est_wrap(paste("vs", lv[1]), 20))) +
    pub_theme(ap$theme) +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
  plot <- ggpubr::ggarrange(top, bottom, ncol = 1, heights = c(1.25, 1), align = "v")
  list(plot = plot, width = max(5.2, 1.15 * k + 1.6), height = 6.4)
}

.est_methods <- function(rows, ref, paired, spec, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  each <- paste(sprintf("%s: mean difference %s (95%% CI %s to %s), Hedges' g = %.2f, t(%s) = %.2f, P = %s",
                        rows$group, sprintf("%.3g", rows$mean_diff), sprintf("%.3g", rows$ci_low),
                        sprintf("%.3g", rows$ci_high), rows$hedges_g, sprintf("%.3g", rows$df),
                        rows$statistic, sprintf("%.2g", rows$p_value)), collapse = "; ")
  paste0(
    sprintf("%s was analysed with estimation statistics, with %s as the reference group (n = %d). ",
            cap_first(y_lab), ref, rows$n_ref[1]),
    sprintf("For each group the %smean difference against %s (group minus reference) is reported with a 95%% bias-corrected and accelerated (BCa) bootstrap confidence interval from %d resamples. ",
            if (paired) "paired " else "", ref, .est_n_boot),
    "The resamples were drawn with the Park-Miller minimal standard generator (multiplier 48271, modulus 2^31 - 1, seed 20261029, restarted for each comparison), ",
    if (paired) "each resample drawing the subject indices of the pairs as floor(u n). " else
      "each resample drawing the reference indices first and then the group indices, as floor(u n). ",
    "The bias correction is the normal quantile of the proportion of bootstrap differences strictly below the estimate, ",
    sprintf("the acceleration comes from the jackknife of the mean difference over every %s, and the interval ends are type-7 quantiles. ",
            if (paired) "pair" else "observation of both groups"),
    sprintf("Hedges' g is Cohen's d over the %s with the exact small-sample correction. ",
            if (paired) "average standard deviation of the two conditions" else "pooled standard deviation"),
    .est_p_text(if (paired) "paired t-test" else "Welch two-sample t-test", nrow(rows), ref),
    each, ". ",
    sprintf("Computed in R %s and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.est_emit_script <- function(spec, lv, paired, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y; id <- spec$data$id
  lv_txt <- paste(sprintf('"%s"', lv), collapse = ", ")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder. It recomputes the BCa bootstrap",
    "# intervals with the same Park-Miller resamples and draws a simple version of",
    "# the estimation plot.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    "M <- 2147483647; A <- 48271; SEED <- 20261029; B <- 5000",
    "mulmod <- function(z, a) {",
    "  hi <- floor(a / 65536); lo <- a - hi * 65536",
    "  ((((z * hi) %% M) * 65536) %% M + (z * lo) %% M) %% M",
    "}",
    "uniforms <- function(n) {",
    "  s <- SEED; out <- numeric(n)",
    "  for (i in seq_len(n)) { s <- mulmod(s, A); out[i] <- s }",
    "  out / M",
    "}",
    "bca <- function(r, g, paired) {",
    "  per <- if (paired) length(r) else length(r) + length(g)",
    "  u <- matrix(uniforms(B * per), nrow = per)",
    "  if (paired) {",
    "    i <- floor(u * length(r)) + 1",
    "    boot <- colSums(matrix(g[i] - r[i], nrow = length(r))) / length(r)",
    "    d <- g - r; jack <- (sum(d) - d) / (length(d) - 1)",
    "  } else {",
    "    ir <- floor(u[seq_along(r), , drop = FALSE] * length(r)) + 1",
    "    ig <- floor(u[length(r) + seq_along(g), , drop = FALSE] * length(g)) + 1",
    "    boot <- colSums(matrix(g[ig], nrow = length(g))) / length(g) -",
    "      colSums(matrix(r[ir], nrow = length(r))) / length(r)",
    "    jack <- c(mean(g) - (sum(r) - r) / (length(r) - 1), (sum(g) - g) / (length(g) - 1) - mean(r))",
    "  }",
    "  est <- sum(g) / length(g) - sum(r) / length(r)",
    "  p0 <- sum(boot < est - 1e-11 * max(abs(c(r, g)))) / B",
    "  z0 <- qnorm(min(max(p0, 0.5 / B), 1 - 0.5 / B))",
    "  dj <- mean(jack) - jack; acc <- sum(dj^3) / (6 * sum(dj^2)^1.5)",
    "  zq <- qnorm(c(0.025, 0.975))",
    "  c(est, quantile(boot, pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq))), type = 7))",
    "}",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(group = as.character(raw[["%s"]]), y = as.numeric(raw[["%s"]])%s)',
            x, y, if (paired) sprintf(', id = as.character(raw[["%s"]])', id) else ""),
    "d <- d[!is.na(d$group) & is.finite(d$y), ]",
    sprintf("lv <- c(%s)", lv_txt),
    "res <- NULL",
    "for (g in lv[-1]) {",
    if (paired) c(
      '  m <- merge(d[d$group == lv[1], c("id", "y")], d[d$group == g, c("id", "y")], by = "id")',
      "  m <- m[order(m$id, method = \"radix\"), ]",
      "  r <- m$y.x; v <- m$y.y"
    ) else "  r <- d$y[d$group == lv[1]]; v <- d$y[d$group == g]",
    sprintf("  ci <- bca(r, v, %s)", if (paired) "TRUE" else "FALSE"),
    sprintf("  tt <- t.test(v, r, paired = %s)", if (paired) "TRUE" else "FALSE"),
    "  res <- rbind(res, data.frame(group = g, mean_diff = ci[1], ci_low = ci[2], ci_high = ci[3],",
    "                               p_value = tt$p.value))",
    "}",
    "print(res, digits = 6)",
    "d$group <- factor(d$group, levels = lv)",
    "p1 <- ggplot(d, aes(group, y, colour = group)) + geom_jitter(width = 0.15) +",
    sprintf('  labs(x = NULL, y = "%s") + theme_classic() + theme(legend.position = "none")',
            spec$appearance$y_label %||% y),
    "res$group <- factor(res$group, levels = lv)",
    "p2 <- ggplot(res, aes(group, mean_diff)) + geom_hline(yintercept = 0, linetype = 2) +",
    "  geom_pointrange(aes(ymin = ci_low, ymax = ci_high)) +",
    "  scale_x_discrete(drop = FALSE) +",
    sprintf('  labs(x = NULL, y = "Mean difference vs %s") + theme_classic()', lv[1]),
    sprintf('ggsave("%s.pdf", ggpubr::ggarrange(p1, p2, ncol = 1, align = "v"), width = 5, height = 6)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}

# Each group is tested against the reference on its own, so with two or more
# groups the methods paragraph says that the P values carry no adjustment.
.est_p_text <- function(test, k, ref) {
  if (k > 1) {
    sprintf("A %s gives each P value, with no adjustment for the %d comparisons against %s. ", test, k, ref)
  } else sprintf("A %s gives the P value. ", test)
}
