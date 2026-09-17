# Load the pubready R engine from the mounted source tree.
#
# core/r is a source R package (DESCRIPTION + R/ + tests). At runtime the CLI
# mounts it read-only and the entry scripts source this file, which loads every
# function in R/ into the global environment. Every cross-package call in R/ is
# namespace-qualified (pkg::fn), so file order does not matter and no NAMESPACE
# processing is needed for the sourced path. testthat drives the same R/ files
# through tests/run_tests.R, so what the tests exercise is what runs in figkit.
if (!exists("CORE_R")) {
  CORE_R <- file.path(Sys.getenv("PUBREADY_CORE", "/opt/pubready/core"), "r")
}
for (f in sort(list.files(file.path(CORE_R, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}
