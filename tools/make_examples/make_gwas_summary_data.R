#!/usr/bin/env Rscript
# Generate the committed manhattan example: synthetic GWAS summary statistics
# for 20,000 SNPs on the 22 autosomes and chromosome X, with the chromosome
# lengths of GRCh38 and the "chr" prefix of a UCSC-style file. The null SNPs
# carry a slight inflation (lambda_GC near 1.03), four loci pass the genome-wide
# line and one locus passes only the suggestive line. Deterministic (seeded) and
# base R only. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_gwas_summary_data.R
set.seed(20261029)

# Sorted random positions of the SNPs of one chromosome.
draw_positions <- function(i, chr_len, per_chr) sort(sample.int(chr_len[i], per_chr[i]))

chr_len <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
             145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
             101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
             50818468, 156040895)
chr_name <- c(as.character(1:22), "X")
n_null <- 19600
per_chr <- as.vector(stats::rmultinom(1, n_null, chr_len / sum(chr_len)))
chr <- rep(chr_name, per_chr)
pos <- unlist(lapply(seq_along(chr_name), draw_positions, chr_len = chr_len, per_chr = per_chr))
p <- stats::pchisq(stats::rchisq(n_null, 1) * 1.03, 1, lower.tail = FALSE)

# Each locus: 80 SNPs within 400 kb of the peak; -log10 p decays with distance.
loci <- data.frame(chr = c("2", "6", "11", "19", "15"),
                   peak = c(60.1e6, 32.6e6, 113.4e6, 44.9e6, 67.3e6),
                   top = c(14, 23, 9.5, 8.6, 6.1))
for (k in seq_len(nrow(loci))) {
  d <- stats::runif(80, -4e5, 4e5)
  lp <- loci$top[k] * exp(-abs(d) / 9e4) + stats::rexp(80, 3)
  chr <- c(chr, rep(loci$chr[k], 80))
  pos <- c(pos, round(loci$peak[k] + d))
  p <- c(p, 10^(-lp))
}

o <- order(match(chr, chr_name), pos)
out <- data.frame(snp = sprintf("rs%d", sample(1e6:9e7, length(p))[o]),
                  chr = paste0("chr", chr[o]), pos = pos[o], p = signif(p[o], 4),
                  stringsAsFactors = FALSE)
out$p[out$p > 1] <- 1
utils::write.csv(out, "example/gwas_summary.csv", row.names = FALSE)
cat(sprintf("wrote example/gwas_summary.csv (%d SNPs, %d below 5e-8)\n",
            nrow(out), sum(out$p < 5e-8)))
