# Recipe: ROC curves of one or more numeric markers against a binary outcome.
# The outcome is --y and the markers are --x, comma-separated. The case level is
# --positive, or the second level in sorted order (so 1 is the case of a 0/1
# outcome). For each marker, pROC computes the curve, the area under it (AUC)
# with a DeLong 95% confidence interval, and the threshold that maximises the
# Youden index. With two or more markers measured on the same subjects, each
# marker is compared with the first by DeLong's test for correlated ROC curves.
#
# The direction follows pROC: a marker whose control median is above its case
# median is read as "lower values indicate the case" (direction ">"), so every
# AUC is at least 0.5. The Python engine restates the same rules.

# Sort strings by code point, which is how Python sorts, so both engines pick
# the same default case level whatever the locale.
.roc_levels <- function(v) sort(unique(v), method = "radix")

# Sensitivity and specificity at the Youden-optimal threshold. Ties go to the
# first threshold in pROC's ascending order.
.roc_youden <- function(r) {
  j <- r$sensitivities + r$specificities - 1
  i <- which.max(j)
  list(threshold = r$thresholds[i], sensitivity = r$sensitivities[i],
       specificity = r$specificities[i])
}

.roc_curve_points <- function(k, rocs, labels) {
  r <- rocs[[k]]
  data.frame(marker = labels[k], fpr = 1 - r$specificities, tpr = r$sensitivities,
             stringsAsFactors = FALSE)
}

recipe_roc_curve <- function(df, spec) {
  y <- spec$data$y
  x_raw <- spec$data$x
  if (is.null(y) || is.null(x_raw)) {
    stop("roc_curve needs --y (the binary outcome) and --x (one marker, or several separated by commas)")
  }
  markers <- trimws(strsplit(x_raw, ",")[[1]])
  markers <- markers[nzchar(markers)]
  for (cn in c(y, markers)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(outcome = as.character(df[[y]]), stringsAsFactors = FALSE)
  for (m in markers) d[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  keep <- !is.na(d$outcome) & stats::complete.cases(d[, markers, drop = FALSE])
  d <- d[keep, , drop = FALSE]
  lv <- .roc_levels(d$outcome)
  if (length(lv) != 2) {
    stop(sprintf("roc_curve needs exactly 2 outcome levels in '%s'; found %d (%s)",
                 y, length(lv), paste(lv, collapse = ", ")))
  }
  case <- as.character(spec$data$positive %||% lv[2])
  if (!case %in% lv) {
    stop(sprintf("--positive '%s' is not a level of '%s' (%s)", case, y, paste(lv, collapse = ", ")))
  }
  control <- setdiff(lv, case)
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing outcome or marker value, so every marker is read on the same subjects.", sum(!keep)),
    sprintf("Kept %d of %d row(s): %d case(s) (%s) and %d control(s) (%s).",
            nrow(d), n0, sum(d$outcome == case), case, sum(d$outcome == control), control))
  n_case <- sum(d$outcome == case); n_ctrl <- sum(d$outcome == control)
  if (n_case < 2 || n_ctrl < 2) stop("roc_curve needs at least 2 cases and 2 controls")

  rocs <- lapply(markers, .roc_fit_one, d = d, control = control, case = case)
  names(rocs) <- markers
  rows <- do.call(rbind, lapply(markers, .roc_row, rocs = rocs, n_case = n_case,
                                n_ctrl = n_ctrl, case = case))
  rows$delong_z_vs_first <- NA_real_
  rows$delong_p_vs_first <- NA_real_
  if (length(markers) >= 2) {
    for (k in 2:length(markers)) {
      rt <- pROC::roc.test(rocs[[1]], rocs[[k]], method = "delong", paired = TRUE)
      rows$delong_z_vs_first[k] <- unname(rt$statistic)
      rows$delong_p_vs_first[k] <- rt$p.value
    }
  }

  legend_lab <- sprintf("%s: AUC %.2f (%s)", markers, rows$auc,
                        .roc_ci_label(rows$auc_ci_low, rows$auc_ci_high))
  pts <- do.call(rbind, lapply(seq_along(markers), .roc_curve_points, rocs = rocs,
                               labels = markers))
  pts$marker <- factor(pts$marker, levels = markers)
  best <- data.frame(marker = factor(markers, levels = markers),
                     fpr = 1 - rows$specificity, tpr = rows$sensitivity)
  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = length(markers))
  sub <- sprintf("%d %s vs %d %s", n_case, case, n_ctrl, control)
  if (length(markers) == 2) {
    sub <- sprintf("%s; DeLong P = %s", sub, sprintf("%.2g", rows$delong_p_vs_first[2]))
  }

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = fpr, y = tpr, colour = marker)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60",
                         linewidth = 0.4) +
    ggplot2::geom_path(linewidth = 0.8) +
    ggplot2::geom_point(data = best, size = 2.4, shape = 21, fill = "white", stroke = 1) +
    ggplot2::scale_colour_manual(values = pal, labels = legend_lab, breaks = markers) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = spec$appearance$x_label %||% "1 - Specificity",
                  y = spec$appearance$y_label %||% "Sensitivity",
                  colour = NULL, title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
                   legend.justification = c(1, 0),
                   legend.text = ggplot2::element_text(size = 9),
                   legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  test_meta <- list(
    name = "ROC curve analysis (DeLong)",
    case_level = case, control_level = control, n_cases = n_case, n_controls = n_ctrl,
    markers = stats::setNames(lapply(seq_len(nrow(rows)), .roc_meta_row, rows = rows), markers))
  resolved <- list(family = "roc", label = test_meta$name)
  methods <- .roc_methods(rows, case, control, n_case, n_ctrl, y)
  build_script <- list(fn = .roc_emit_script,
                       args = list(spec = spec, markers = markers, case = case, control = control))
  list(plot = p, stats = rows, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(paste(c(y, markers), collapse = "_")),
       df_used = d, clean_steps = clean_steps, width = 4.8, height = 5.0)
}

.roc_fit_one <- function(m, d, control, case) {
  pROC::roc(response = d$outcome, predictor = d[[m]], levels = c(control, case),
            direction = "auto", quiet = TRUE)
}

.roc_row <- function(m, rocs, n_case, n_ctrl, case) {
  r <- rocs[[m]]
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  yd <- .roc_youden(r)
  auc <- as.numeric(pROC::auc(r))
  # The DeLong CI has zero width at an AUC of 0 or 1 and for a constant marker,
  # so it is not valid there and is not reported; a constant marker has no
  # optimal threshold (F86). The Python engine has the same rule.
  constant <- length(unique(c(r$cases, r$controls))) == 1
  separated <- auc >= 1 || auc <= 0
  if (constant || separated) ci <- c(NA_real_, NA_real_, NA_real_)
  data.frame(marker = m, case_level = case, n_cases = n_case, n_controls = n_ctrl,
             direction = r$direction, auc = auc,
             auc_ci_low = ci[1], auc_ci_high = ci[3],
             youden_threshold = if (constant) NA_real_ else yd$threshold,
             sensitivity = yd$sensitivity, specificity = yd$specificity,
             constant_marker = constant, complete_separation = separated,
             stringsAsFactors = FALSE)
}

# The CI part of a legend label: the bounds, or "CI not valid" where the
# DeLong CI is not reported.
.roc_ci_label <- function(lo, hi) {
  ifelse(is.na(lo), "CI not valid", sprintf("%.2f to %.2f", lo, hi))
}

# The methods text of one marker.
.roc_marker_text <- function(i, rows) {
  r <- rows[i, ]
  if (r$constant_marker) {
    return(sprintf("%s is constant, so it does not discriminate the outcome: AUC = %.3f, with no valid CI and no optimal threshold",
                   r$marker, r$auc))
  }
  ci <- if (r$complete_separation) {
    "the DeLong 95% CI is not valid, because the marker separates the cases from the controls completely"
  } else sprintf("95%% CI %.3f to %.3f", r$auc_ci_low, r$auc_ci_high)
  sprintf("%s, AUC = %.3f (%s), optimal threshold %s (sensitivity %.2f, specificity %.2f)",
          r$marker, r$auc, ci, sprintf("%.4g", r$youden_threshold), r$sensitivity, r$specificity)
}

.roc_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(auc = r$auc, ci_low = r$auc_ci_low, ci_high = r$auc_ci_high, direction = r$direction,
       youden_threshold = r$youden_threshold, sensitivity = r$sensitivity,
       specificity = r$specificity, delong_p_vs_first = r$delong_p_vs_first,
       constant_marker = r$constant_marker, complete_separation = r$complete_separation)
}

.roc_methods <- function(rows, case, control, n_case, n_ctrl, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions(c("pROC", "ggplot2"))
  per <- paste(vapply(seq_len(nrow(rows)), .roc_marker_text, character(1), rows = rows),
               collapse = "; ")
  cmp <- if (nrow(rows) >= 2) {
    k <- 2:nrow(rows)
    paste0(" The AUC of each marker was compared with that of ", rows$marker[1],
           " by DeLong's test for two correlated ROC curves: ",
           paste(sprintf("%s, Z = %.2f, P = %s", rows$marker[k], rows$delong_z_vs_first[k],
                         sprintf("%.2g", rows$delong_p_vs_first[k])), collapse = "; "), ".")
  } else ""
  flip <- if (any(rows$direction == ">")) {
    sprintf(" For %s, lower values indicate the case, so the curve reads the marker in that direction.",
            paste(rows$marker[rows$direction == ">"], collapse = ", "))
  } else ""
  paste0(
    sprintf("Receiver operating characteristic (ROC) curves were computed for the discrimination of %s = %s (cases, n = %d) from %s (controls, n = %d). ",
            y, case, n_case, control, n_ctrl),
    "The area under each curve (AUC) is reported with a 95% confidence interval by the method of DeLong et al. (1988), and the optimal threshold maximises the Youden index (sensitivity + specificity - 1). ",
    per, ".", flip, cmp,
    sprintf(" Computed in R %s with pROC %s and drawn with ggplot2 %s.", rver, v$pROC, v$ggplot2))
}

.roc_emit_script <- function(spec, markers, case, control, in_name, fig_stub) {
  y <- spec$data$y
  mk <- paste0("c(", paste(sprintf('"%s"', markers), collapse = ", "), ")")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(pROC); library(ggplot2) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf("markers <- %s", mk),
    sprintf('raw <- raw[!is.na(raw[["%s"]]) & complete.cases(raw[, markers, drop = FALSE]), ]', y),
    "rocs <- list(); pts <- NULL",
    "for (m in markers) {",
    sprintf('  r <- roc(as.character(raw[["%s"]]), raw[[m]], levels = c("%s", "%s"), direction = "auto", quiet = TRUE)',
            y, control, case),
    "  print(r); print(ci.auc(r, method = \"delong\"))",
    "  rocs[[m]] <- r",
    "  pts <- rbind(pts, data.frame(marker = m, fpr = 1 - r$specificities, tpr = r$sensitivities))",
    "}",
    "if (length(rocs) >= 2) print(roc.test(rocs[[1]], rocs[[2]], method = \"delong\", paired = TRUE))",
    "p <- ggplot(pts, aes(fpr, tpr, colour = marker)) + geom_abline(linetype = \"dashed\") +",
    '  geom_path() + coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity") + theme_classic()',
    sprintf('ggsave("%s.pdf", p, width = 4.8, height = 5.0)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
