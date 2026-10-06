#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

raw <- read.csv("input_status_marker_a_marker_b.csv", check.names = FALSE)
markers <- c("marker_a", "marker_b")
raw <- raw[!is.na(raw[["status"]]) & complete.cases(raw[, markers, drop = FALSE]), ]
rocs <- list(); pts <- NULL
for (m in markers) {
  r <- roc(as.character(raw[["status"]]), raw[[m]], levels = c("healthy", "disease"), direction = "auto", quiet = TRUE)
  print(r); print(ci.auc(r, method = "delong"))
  rocs[[m]] <- r
  pts <- rbind(pts, data.frame(marker = m, fpr = 1 - r$specificities, tpr = r$sensitivities))
}
if (length(rocs) >= 2) print(roc.test(rocs[[1]], rocs[[2]], method = "delong", paired = TRUE))
p <- ggplot(pts, aes(fpr, tpr, colour = marker)) + geom_abline(linetype = "dashed") +
  geom_path() + coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity") + theme_classic()
ggsave("figure_status_marker_a_marker_b_roc.pdf", p, width = 4.8, height = 5.0)
cat("Reproduced figure.\n")
