#!/usr/bin/env Rscript
# Generate the committed embedding example: a synthetic UMAP of 2,590 blood
# cells in eight cell types, shaped like the PBMC 3k data of 10x Genomics. The
# coordinates are drawn, not computed: each cell type is an elongated blob, and
# the two T-cell types touch, so their silhouette widths are lower. The file
# holds the cell barcode, the two coordinates, the cell type and a numeric
# cluster id, as a Seurat export does. Deterministic (seeded) and base R only.
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_umap_clusters_data.R
set.seed(20261029)

types <- data.frame(
  name  = c("CD4 T cells", "CD14 monocytes", "B cells", "CD8 T cells", "NK cells",
            "FCGR3A monocytes", "Dendritic cells", "Platelets"),
  n     = c(1100, 480, 340, 300, 160, 160, 35, 15),
  cx    = c(-4.0, 6.5, -1.0, -6.2, -7.8, 8.6, 4.2, 1.5),
  cy    = c(1.0, 2.5, -7.5, -1.2, -4.6, -0.8, 6.4, 9.0),
  sx    = c(1.6, 1.3, 0.9, 0.9, 0.7, 0.6, 0.5, 0.3),
  sy    = c(1.1, 1.8, 1.0, 1.2, 0.8, 0.9, 0.5, 0.3),
  angle = c(0.4, -0.3, 0.1, 0.9, 0.2, 0.5, 0, 0),
  stringsAsFactors = FALSE)

# The cells of one type: a rotated, stretched Gaussian blob with a slight bend.
type_cells <- function(i, types) {
  u <- stats::rnorm(types$n[i], 0, types$sx[i])
  v <- stats::rnorm(types$n[i], 0, types$sy[i]) + 0.12 * u^2
  a <- types$angle[i]
  data.frame(UMAP_1 = types$cx[i] + u * cos(a) - v * sin(a),
             UMAP_2 = types$cy[i] + u * sin(a) + v * cos(a),
             cell_type = types$name[i], seurat_clusters = i - 1L, stringsAsFactors = FALSE)
}

# One random 14-base cell barcode.
random_barcode <- function(i, bases) paste(sample(bases, 14, replace = TRUE), collapse = "")

cells <- do.call(rbind, lapply(seq_len(nrow(types)), type_cells, types = types))
cells <- cells[sample.int(nrow(cells)), ]
barcode <- vapply(seq_len(nrow(cells)), random_barcode, character(1), bases = c("A", "C", "G", "T"))
out <- data.frame(cell = paste0(make.unique(barcode), "-1"),
                  UMAP_1 = round(cells$UMAP_1, 4), UMAP_2 = round(cells$UMAP_2, 4),
                  cell_type = cells$cell_type, seurat_clusters = cells$seurat_clusters,
                  stringsAsFactors = FALSE)
utils::write.csv(out, "example/umap_clusters.csv", row.names = FALSE)
cat(sprintf("wrote example/umap_clusters.csv (%d cells, %d cell types)\n", nrow(out), nrow(types)))
