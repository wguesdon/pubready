#!/usr/bin/env Rscript
# Generate a clusterProfiler-style ORA enrichment table for the enrichment_dot
# recipe. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.4.0 \
#     Rscript tools/make_examples/make_enrichment_data.R
set.seed(19)

terms <- c(
  "Cytokine signaling in immune system", "Inflammatory response",
  "T cell receptor signaling", "Interferon gamma response",
  "Complement activation", "Apoptotic process", "Cell cycle checkpoint",
  "Oxidative phosphorylation", "Fatty acid metabolism", "Glycolysis",
  "Extracellular matrix organization", "Angiogenesis", "Wnt signaling",
  "p53 signaling pathway", "DNA repair", "Antigen processing",
  "Chemokine signaling", "NF-kappa B signaling", "Hypoxia response",
  "Epithelial mesenchymal transition")

n <- length(terms)
count <- sample(6:60, n, replace = TRUE)
set_total <- 400L
padj <- sort(10^(-runif(n, 1.5, 8)))          # most-significant first
generatio <- sprintf("%d/%d", count, set_total)

en <- data.frame(
  Description = terms,
  GeneRatio = generatio,
  Count = count,
  p.adjust = signif(padj, 3),
  stringsAsFactors = FALSE)
write.csv(en, "example/enrichment_results.csv", row.names = FALSE)
cat(sprintf("wrote example/enrichment_results.csv (%d terms)\n", n))
