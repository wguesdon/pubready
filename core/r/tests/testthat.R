# Standard testthat entry point for R CMD check / devtools::test().
# The container runs the suite through tests/run_tests.R instead, which loads
# the mounted source without installing.
library(testthat)
library(pubready)

test_check("pubready")
