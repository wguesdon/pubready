# Recipe: Sankey (alluvial) diagram of patient state transitions across visits,
# for example the response category at each assessment.
# Input is one row per patient visit: the patient (--id), the visit (--x) and the
# state (--y). Visits are sorted numerically when every visit is a number, else
# in code-point order. States follow --levels when it is given, else they are
# sorted numerically when every state is a number, else they keep the order of
# their first appearance in the data. A
# patient who is missing at a visit takes no part in the transitions into or out
# of that visit, so the node of that visit holds only the patients present.
#
# The diagram is drawn here, without ggalluvial: one column of stacked nodes per
# visit (height proportional to the count, small gaps), and a ribbon for each
# transition between consecutive visits, shaped by the cubic 3u^2 - 2u^3 and
# coloured by the state of origin. The Python engine restates the same layout.
#
# With exactly two visits that carry the same set of states, Bowker's test of
# symmetry (McNemar's test for two states, no continuity correction) asks
# whether the patients move between states in one direction more than in the
# other. A pair of cells that holds no patient adds nothing to the statistic and
# no degree of freedom. In every other case the diagram is descriptive.

# Node and layout constants, matched to the Python engine.
.sk_hw <- 0.045          # half width of a node, in visit units
.sk_gap_frac <- 0.03     # gap between nodes, as a fraction of the largest visit
.sk_label_in <- 0.17     # labels in one column stay at least this far apart (in)

.sk_labels <- function(v) {
  out <- trimws(as.character(v))
  out[is.na(v)] <- NA_character_
  out
}

# Visit order: numeric when every value is a number, else code-point order.
.sk_visit_order <- function(v) {
  u <- unique(v)
  num <- suppressWarnings(as.numeric(u))
  if (all(!is.na(num))) u[order(num, u, method = "radix")] else sort(u, method = "radix")
}

# State order: the order of --levels when it is given, else numeric when every
# value is a number, else first appearance. A state of the data that --levels
# does not list stops the run; a listed state that no row holds is left out.
.sk_state_order <- function(s, levels = NULL) {
  u <- unique(s)
  if (length(levels) > 0) {
    levels <- unique(trimws(as.character(levels)))
    miss <- setdiff(u, levels)
    if (length(miss) > 0) {
      stop(sprintf("--levels does not list the state(s) %s of the data",
                   paste(sprintf("'%s'", miss), collapse = ", ")))
    }
    return(levels[levels %in% u])
  }
  num <- suppressWarnings(as.numeric(u))
  if (all(!is.na(num))) u[order(num, u, method = "radix")] else u
}

# The transitions between visit k and visit k + 1 among the patients present at both.
.sk_pair <- function(k, d, visits, states) {
  a <- d[d$visit == visits[k], c("id", "state")]
  b <- d[d$visit == visits[k + 1], c("id", "state")]
  m <- merge(a, b, by = "id", suffixes = c("_from", "_to"))
  tab <- table(factor(m$state_from, levels = states), factor(m$state_to, levels = states))
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(out) <- c("from_state", "to_state", "n")
  out <- out[out$n > 0, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out$from_state <- factor(out$from_state, levels = states)
  out$to_state <- factor(out$to_state, levels = states)
  out <- out[order(out$from_state, out$to_state), , drop = FALSE]
  data.frame(from_visit = visits[k], to_visit = visits[k + 1],
             from_state = as.character(out$from_state), to_state = as.character(out$to_state),
             n = as.integer(out$n), stringsAsFactors = FALSE)
}

# Bowker's test (McNemar's for two states) of the square table of two visits.
.sk_symmetry <- function(tab) {
  k <- nrow(tab)
  up <- tab[upper.tri(tab)]
  lo <- t(tab)[upper.tri(tab)]
  keep <- (up + lo) > 0
  if (sum(keep) == 0) return(NULL)
  if (all(keep)) {
    mt <- stats::mcnemar.test(tab, correct = FALSE)
    return(list(statistic = unname(mt$statistic), df = unname(mt$parameter),
                p_value = mt$p.value, empty_pairs = 0L))
  }
  stat <- sum((up[keep] - lo[keep])^2 / (up[keep] + lo[keep]))
  df <- sum(keep)
  list(statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE),
       empty_pairs = sum(!keep))
}

.sk_test <- function(d, visits, states) {
  if (length(visits) != 2) return(list(reason = "more than two visits"))
  a <- d[d$visit == visits[1], c("id", "state")]
  b <- d[d$visit == visits[2], c("id", "state")]
  m <- merge(a, b, by = "id", suffixes = c("_from", "_to"))
  s1 <- states[states %in% m$state_from]
  s2 <- states[states %in% m$state_to]
  if (!setequal(s1, s2)) return(list(reason = "the two visits do not carry the same states"))
  if (length(s1) < 2) return(list(reason = "fewer than two states"))
  tab <- table(factor(m$state_from, levels = s1), factor(m$state_to, levels = s1))
  res <- .sk_symmetry(unclass(tab))
  if (is.null(res)) return(list(reason = "no patient changed state"))
  res$name <- if (length(s1) == 2) "McNemar's test" else "Bowker's test of symmetry"
  res$n_paired <- nrow(m)
  res$states <- s1
  res$table <- unclass(tab)
  res
}

# Node positions: for each visit, the states present stacked from the top in
# state order, each column centred on the tallest one.
.sk_nodes <- function(d, visits, states) {
  counts <- table(factor(d$visit, levels = visits), factor(d$state, levels = states))
  nmax <- max(rowSums(counts))
  gap <- .sk_gap_frac * nmax
  heights <- vapply(seq_along(visits), .sk_col_height, numeric(1), counts = counts, gap = gap)
  hmax <- max(heights)
  rows <- lapply(seq_along(visits), .sk_col_nodes, counts = counts, gap = gap,
                 heights = heights, hmax = hmax, visits = visits, states = states)
  list(nodes = do.call(rbind, rows), hmax = hmax, gap = gap)
}

.sk_col_height <- function(v, counts, gap) {
  cn <- counts[v, ]
  sum(cn) + gap * (sum(cn > 0) - 1)
}

.sk_col_nodes <- function(v, counts, gap, heights, hmax, visits, states) {
  cn <- as.vector(counts[v, ])
  top <- hmax - (hmax - heights[v]) / 2
  out <- NULL
  for (s in seq_along(states)) {
    if (cn[s] == 0) next
    out <- rbind(out, data.frame(visit = visits[v], vi = v - 1, state = states[s], n = cn[s],
                                 ytop = top, ybot = top - cn[s], stringsAsFactors = FALSE))
    top <- top - cn[s] - gap
  }
  out
}

# Ribbon polygons: outflows leave a node from its top in the order of the target
# state, and inflows enter a node from its top in the order of the source state.
.sk_ribbons <- function(trans, nodes, states) {
  if (nrow(trans) == 0) return(NULL)
  out_off <- stats::setNames(rep(0, nrow(nodes)), paste(nodes$visit, nodes$state))
  in_off <- out_off
  trans$src_key <- paste(trans$from_visit, trans$from_state)
  trans$dst_key <- paste(trans$to_visit, trans$to_state)
  ntop <- stats::setNames(nodes$ytop, paste(nodes$visit, nodes$state))
  nvi <- stats::setNames(nodes$vi, paste(nodes$visit, nodes$state))
  o_out <- order(match(trans$from_visit, unique(nodes$visit)),
                 match(trans$from_state, states), match(trans$to_state, states))
  y0 <- numeric(nrow(trans))
  for (i in o_out) {
    y0[i] <- ntop[trans$src_key[i]] - out_off[trans$src_key[i]]
    out_off[trans$src_key[i]] <- out_off[trans$src_key[i]] + trans$n[i]
  }
  o_in <- order(match(trans$to_visit, unique(nodes$visit)),
                match(trans$to_state, states), match(trans$from_state, states))
  y1 <- numeric(nrow(trans))
  for (i in o_in) {
    y1[i] <- ntop[trans$dst_key[i]] - in_off[trans$dst_key[i]]
    in_off[trans$dst_key[i]] <- in_off[trans$dst_key[i]] + trans$n[i]
  }
  polys <- lapply(seq_len(nrow(trans)), .sk_ribbon_poly, trans = trans, y0 = y0, y1 = y1,
                  x0 = nvi[trans$src_key] + .sk_hw, x1 = nvi[trans$dst_key] - .sk_hw)
  do.call(rbind, polys)
}

.sk_ribbon_poly <- function(i, trans, y0, y1, x0, x1) {
  u <- seq(0, 1, length.out = 40)
  s <- 3 * u^2 - 2 * u^3
  xs <- x0[i] + (x1[i] - x0[i]) * u
  top <- y0[i] + (y1[i] - y0[i]) * s
  bot <- (y0[i] - trans$n[i]) + ((y1[i] - trans$n[i]) - (y0[i] - trans$n[i])) * s
  data.frame(flow = i, state = trans$from_state[i], x = c(xs, rev(xs)), y = c(top, rev(bot)),
             stringsAsFactors = FALSE)
}

# Figure size and label room, matched to the Python engine.
.sk_dims <- function(nodes, visits, states, ncol_leg) {
  chars <- max(nchar(sprintf("%s (%d)", nodes$state, nodes$n)))
  label_in <- 0.075 * chars + 0.25
  nv <- length(visits)
  # Visits sit far enough apart for their axis labels, and for the node labels
  # of a middle visit, which sit between two columns.
  unit_in <- max(2.0, 0.09 * max(nchar(visits)) + 0.4, if (nv > 2) label_in + 0.6 else 0)
  legend_rows <- ceiling(length(states) / ncol_leg)
  legend_h <- 0.3 + 0.25 * legend_rows
  height <- max(4.2, 2.4 + 0.32 * length(states)) + legend_h
  pad <- .sk_hw + 0.03 + label_in / unit_in
  list(width = unit_in * (nv - 1) + 2 * (.sk_hw * unit_in + 0.03 * unit_in + label_in) + 0.6,
       height = height, panel_h = height - 1.5 - legend_h, xlim = c(-pad, nv - 1 + pad),
       unit_in = unit_in)
}

.sk_legend_ncol <- function(states, width) {
  per <- 0.085 * max(nchar(states)) + 0.6
  ncol <- min(length(states), max(1, floor((width - 0.5) / per)))
  as.integer(ceiling(length(states) / ceiling(length(states) / ncol)))
}

# The subtitle with the test result. A subtitle wider than the figure breaks
# after the patient count.
.sk_subtitle <- function(test, n_patients, n_visits, width) {
  base <- sprintf("%d patients, %d visits", n_patients, n_visits)
  if (is.null(test$statistic)) return(base)
  p <- if (test$p_value < 0.001) "P < 0.001" else sprintf("P = %.2g", test$p_value)
  tst <- sprintf("%s (n = %d): chi-square(%d) = %.2f, %s", test$name, test$n_paired,
                 as.integer(test$df), test$statistic, p)
  one <- sprintf("%s; %s", base, tst)
  if (0.07 * nchar(one) + 0.4 > width) sprintf("%s\n%s", base, tst) else one
}

recipe_sankey <- function(df, spec) {
  idc <- spec$data$id; x <- spec$data$x; y <- spec$data$y
  if (is.null(idc) || !nzchar(idc)) stop("sankey needs --id (the patient column)")
  if (is.null(x) || is.null(y)) stop("sankey needs --x (the visit) and --y (the state)")
  for (cn in c(idc, x, y)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(id = .sk_labels(df[[idc]]), visit = .sk_labels(df[[x]]),
                  state = .sk_labels(df[[y]]), stringsAsFactors = FALSE)
  keep <- !is.na(d$id) & !is.na(d$visit) & !is.na(d$state) & nzchar(d$state)
  d <- d[keep, , drop = FALSE]
  if (anyDuplicated(d[, c("id", "visit")])) {
    stop("sankey needs one row per patient and visit; some patients have two rows at the same visit")
  }
  visits <- .sk_visit_order(d$visit)
  if (length(visits) < 2) stop("sankey needs at least 2 visits")
  levels <- unlist(spec$data$levels)
  states <- .sk_state_order(d$state, levels)
  trans <- do.call(rbind, lapply(seq_len(length(visits) - 1), .sk_pair, d = d,
                                 visits = visits, states = states))
  if (is.null(trans)) {
    trans <- data.frame(from_visit = character(0), to_visit = character(0),
                        from_state = character(0), to_state = character(0), n = integer(0))
  }
  test <- .sk_test(d, visits, states)
  n_pat <- length(unique(d$id))
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing patient, visit or state.", sum(!keep)),
    sprintf("Kept %d row(s) of %d: %d patient(s), %d visit(s) (%s), %d state(s) (%s).",
            nrow(d), n0, n_pat, length(visits), paste(visits, collapse = ", "),
            length(states), paste(states, collapse = ", ")),
    sprintf("%d patient(s) are missing at one visit or more and take no part in the transitions into or out of that visit.",
            sum(table(d$id) < length(visits))))

  lay <- .sk_nodes(d, visits, states)
  nodes <- lay$nodes
  pal <- stats::setNames(pub_palette(length(states), spec$appearance$palette), states)
  width0 <- .sk_dims(nodes, visits, states, length(states))$width
  ncol_leg <- .sk_legend_ncol(states, width0)
  dims <- .sk_dims(nodes, visits, states, ncol_leg)
  sub <- .sk_subtitle(test, n_pat, length(visits), dims$width)
  if (grepl("\n", sub, fixed = TRUE)) dims$height <- dims$height + 0.2
  rib <- .sk_ribbons(trans, nodes, states)
  nodes$label <- sprintf("%s (%d)", nodes$state, nodes$n)
  nodes$ymid <- (nodes$ytop + nodes$ybot) / 2
  nodes$show <- .sk_label_fit(nodes, dims$panel_h / (1.04 * lay$hmax))
  left <- nodes[nodes$show & nodes$vi == 0, , drop = FALSE]
  right <- nodes[nodes$show & nodes$vi > 0, , drop = FALSE]
  nodes$state <- factor(nodes$state, levels = states)

  p <- ggplot2::ggplot()
  if (!is.null(rib)) {
    rib$state <- factor(rib$state, levels = states)
    p <- p + ggplot2::geom_polygon(data = rib, ggplot2::aes(x = x, y = y, group = flow, fill = state),
                                   alpha = 0.45, colour = NA, show.legend = FALSE)
  }
  p <- p +
    ggplot2::geom_rect(data = nodes, ggplot2::aes(xmin = vi - .sk_hw, xmax = vi + .sk_hw,
                                                  ymin = ybot, ymax = ytop, fill = state),
                       colour = NA) +
    ggplot2::geom_label(data = left, ggplot2::aes(x = vi - .sk_hw - 0.03, y = ymid, label = label),
                        hjust = 1, size = 3.2, label.size = 0, label.padding = ggplot2::unit(0.08, "lines"),
                        fill = grDevices::adjustcolor("white", alpha.f = 0.75)) +
    ggplot2::geom_label(data = right, ggplot2::aes(x = vi + .sk_hw + 0.03, y = ymid, label = label),
                        hjust = 0, size = 3.2, label.size = 0, label.padding = ggplot2::unit(0.08, "lines"),
                        fill = grDevices::adjustcolor("white", alpha.f = 0.75)) +
    ggplot2::scale_fill_manual(values = pal, breaks = states, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_along(visits) - 1, labels = visits,
                                limits = dims$xlim, expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(-0.02 * lay$hmax, 1.02 * lay$hmax),
                                expand = ggplot2::expansion(0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_leg, title.position = "top",
                                                 title.hjust = 0.5)) +
    ggplot2::labs(x = spec$appearance$x_label %||% x, y = NULL, fill = y,
                  title = spec$appearance$title,
                  subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 4)),
                   legend.position = "bottom", legend.title = ggplot2::element_text(size = 11),
                   legend.key.size = ggplot2::unit(0.18, "in"),
                   legend.text = ggplot2::element_text(size = 10),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  stats_df <- trans
  rownames(stats_df) <- NULL
  has_test <- !is.null(test$statistic)
  test_meta <- list(
    name = if (has_test) test$name else "Sankey diagram (descriptive)",
    n_patients = n_pat, n_visits = length(visits), visits = as.list(visits),
    states = as.list(states), n_transitions = sum(trans$n))
  if (has_test) {
    test_meta$statistic <- test$statistic
    test_meta$df <- test$df
    test_meta$p_value <- test$p_value
    test_meta$n_paired <- test$n_paired
    test_meta$empty_pairs <- test$empty_pairs
    test_meta$table <- lapply(seq_len(nrow(test$table)), .sk_table_row, tab = test$table,
                              states = test$states)
  } else {
    test_meta$no_test_reason <- test$reason
  }
  resolved <- list(family = "sankey", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .sk_methods(test, trans, visits, states, n_pat, idc, x, y,
                             levels_given = length(levels) > 0),
       build_script = list(fn = .sk_emit_script, args = list(spec = spec, visits = visits,
                                                             states = states, dims = dims)),
       label = slugify(y), df_used = d, clean_steps = clean_steps,
       width = dims$width, height = dims$height)
}

.sk_table_row <- function(i, tab, states) {
  stats::setNames(as.list(as.integer(tab[i, ])), states)
}

.sk_methods <- function(test, trans, visits, states, n_pat, idc, x, y, levels_given = FALSE) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  state_txt <- if (levels_given) {
    sprintf("states in the order set by the analyst (%s). ", paste(states, collapse = ", "))
  } else {
    "states numerically when every state is a number, else by their first appearance in the data. "
  }
  txt <- paste0(
    sprintf("The state %s of %d patients (%s) was followed across %d visits of %s (%s). ",
            y, n_pat, idc, length(visits), x, paste(visits, collapse = ", ")),
    "Visits are ordered numerically when every visit is a number, else alphabetically, and ",
    state_txt,
    sprintf("The diagram shows %d transitions between consecutive visits; a patient missing at a visit takes no part in the transitions into or out of it, and the transitions table gives every count. ",
            sum(trans$n)))
  if (!is.null(test$statistic)) {
    p <- if (test$p_value < 0.001) "P < 0.001" else sprintf("P = %.2g", test$p_value)
    txt <- paste0(txt, sprintf(
      "%s, without continuity correction, compared the %d patients present at both visits (chi-square(%d) = %.2f, %s).",
      test$name, test$n_paired, as.integer(test$df), test$statistic, p))
    if (test$empty_pairs > 0) {
      txt <- paste0(txt, sprintf(" %d pair(s) of cells held no patient and were left out of the statistic and its degrees of freedom.",
                                 test$empty_pairs))
    }
  } else {
    txt <- paste0(txt, sprintf("The diagram is descriptive, and no test was applied because %s.", test$reason))
  }
  paste0(txt, sprintf(" Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.sk_emit_script <- function(spec, visits, states, dims, in_name, fig_stub) {
  idc <- spec$data$id; x <- spec$data$x; y <- spec$data$y
  vs <- paste0("c(", paste(sprintf('"%s"', visits), collapse = ", "), ")")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction of the transition counts and the test. Run inside",
    "# the pinned pubready container (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE, stringsAsFactors = FALSE)', in_name),
    sprintf('d <- data.frame(id = as.character(raw[["%s"]]), visit = trimws(as.character(raw[["%s"]])),', idc, x),
    sprintf('                state = trimws(as.character(raw[["%s"]])))', y),
    "d <- d[complete.cases(d), ]",
    sprintf("visits <- %s", vs),
    "counts <- NULL",
    "for (k in seq_len(length(visits) - 1)) {",
    "  m <- merge(d[d$visit == visits[k], ], d[d$visit == visits[k + 1], ], by = \"id\")",
    "  tab <- table(m$state.x, m$state.y)",
    "  print(tab)",
    "  if (length(visits) == 2 && nrow(tab) == ncol(tab)) print(mcnemar.test(tab, correct = FALSE))",
    "  counts <- rbind(counts, data.frame(visit = visits[k], as.data.frame(table(m$state.x))))",
    "}",
    "p <- ggplot(counts, aes(visit, Freq, fill = Var1)) + geom_col() +",
    sprintf('  labs(x = "%s", y = "Patients", fill = "%s") + theme_classic(base_size = 13)', x, y),
    sprintf('ggsave("%s.pdf", p, width = %.2f, height = %.2f)', fig_stub, dims$width, dims$height),
    'cat("Reproduced the transition counts.\\n")'
  )
}

# A node shows its label when the label centre is at least .sk_label_in below
# the last label shown in the same column, from the top down, so labels never
# overlap. in_per_unit converts the count axis to inches.
.sk_label_fit <- function(nodes, in_per_unit) {
  show <- logical(nrow(nodes))
  for (v in unique(nodes$vi)) {
    last <- Inf
    for (i in which(nodes$vi == v)) {
      y_in <- nodes$ymid[i] * in_per_unit
      if (last - y_in >= .sk_label_in) {
        show[i] <- TRUE
        last <- y_in
      }
    }
  }
  show
}
