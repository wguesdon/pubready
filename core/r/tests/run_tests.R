#!/usr/bin/env Rscript
# pubready R test suite. Loads the mounted source package with pkgload (via
# testthat::test_local) and runs the testthat tests, so the code under test is
# the same code bootstrap.R sources at runtime. Invoked by `figkit test`.
suppressPackageStartupMessages(library(testthat))
core_r <- Sys.getenv("PUBREADY_CORE_R", "")
if (!nzchar(core_r)) core_r <- "/opt/pubready/core/r"
message(sprintf("Testing pubready R engine at %s", core_r))
testthat::test_local(path = core_r, reporter = "summary", stop_on_failure = TRUE)
