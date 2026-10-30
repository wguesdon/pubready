#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ComplexHeatmap))

raw <- read.csv("input_set_membership.csv", check.names = FALSE, row.names = 1)
# 1, TRUE, T, yes or y is a member; any other code or an empty cell is not.
bin <- raw
for (cn in names(raw)) {
  s <- tolower(trimws(as.character(raw[[cn]])))
  bin[[cn]] <- !is.na(s) & s %in% c("1", "true", "t", "yes", "y")
}
bin <- as.data.frame(bin[rowSums(bin) > 0, , drop = FALSE])
m <- make_comb_mat(bin)
cs <- comb_size(m); cd <- comb_degree(m); cn <- comb_name(m)
# The largest intersections: size, then degree, then the set names.
nm <- character(length(cn))
for (i in seq_along(cn)) nm[i] <- paste(set_name(m)[strsplit(cn[i], "")[[1]] == "1"], collapse = " & ")
ord <- order(-cs, cd, nm, method = "radix")
shown <- ord[seq_len(min(length(ord), 15))]
m2 <- m[shown]
up <- UpSet(m2, comb_order = match(cn[shown], comb_name(m2)),
            set_order = match(set_name(m)[order(-set_size(m), seq_along(set_size(m)), method = "radix")], set_name(m2)),
            top_annotation = upset_top_annotation(m2, add_numbers = TRUE, numbers_rot = 0),
            right_annotation = upset_right_annotation(m2, add_numbers = TRUE,
                                                      axis_param = list(labels_rot = 0)))

cairo_pdf("figure_set_membership_upset.pdf", width = 7.96, height = 4.10)
draw(up); invisible(dev.off())
cat("Reproduced figure.\n")
