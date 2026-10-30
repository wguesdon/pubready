#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(cmprsk); library(ggplot2) })

raw <- read.csv("input_status.csv", check.names = FALSE, stringsAsFactors = FALSE)
tm <- as.numeric(raw[["months"]])
status <- trimws(as.character(raw[["status"]]))
grp <- trimws(as.character(raw[["arm"]]))
ok <- is.finite(tm) & tm >= 0 & !is.na(status) & !is.na(grp)
tm <- tm[ok]; status <- status[ok]
grp <- factor(grp[ok], levels = c("intensified", "standard"))
fit <- cuminc(tm, status, grp, cencode = "0")
if (!is.null(fit$Tests)) print(fit$Tests)
print(timepoints(fit, c(20, 30, 50)))
curves <- NULL
for (nm in setdiff(names(fit), "Tests")) {
  curves <- rbind(curves, data.frame(curve = nm, time = fit[[nm]]$time, est = fit[[nm]]$est))
}
p <- ggplot(curves, aes(time, est, colour = curve)) + geom_path() +
  labs(x = "Months after transplant", y = "Cumulative incidence", colour = NULL) + theme_classic(base_size = 13)
ggsave("figure_status_cif.pdf", p, width = 9.35, height = 5.60)
cat("Reproduced figure.\n")
