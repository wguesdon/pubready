# The v3 fixes of the clinical recipes. The Python suite (test_v3_clinical.py)
# runs the same cases and expects the same results.

v3c_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

# Write the input and the script into a new temporary folder, run the script
# with Rscript there, and return the exit status and the output.
v3c_run_script <- function(script, df) {
  dir <- tempfile("script_")
  dir.create(dir)
  utils::write.csv(df, file.path(dir, "input.csv"), row.names = FALSE)
  writeLines(script, file.path(dir, "script.R"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "script.R",
                                  stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status"), out = paste(out, collapse = "\n"))
}

v3c_km_space <- function() {
  data.frame(`follow up months` = c(2, 5, 5, 8, 10, 10, 14, 18, 21.8, 3, 6, 10, 12, 15, 17, 19, 20, 20),
             `death event` = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1),
             `treatment arm` = rep(c("A", "B"), each = 9),
             check.names = FALSE)
}

test_that("survival_km reads columns whose names hold a space (F44)", {
  res <- recipe_survival_km(v3c_km_space(), v3c_spec("survival_km", time = "follow up months",
                                                     event = "death event", x = "treatment arm"))
  expect_equal(res$stats$group, c("A", "B"))
  expect_equal(res$stats$events, c(7L, 6L))
  d <- stats::setNames(v3c_km_space(), c("t", "e", "g"))
  ref <- survival::survdiff(survival::Surv(t, e) ~ g, data = d)
  expect_equal(res$test_meta$log_rank$chisq, unname(ref$chisq), tolerance = 1e-12)
  expect_match(res$methods, "groups of treatment arm")
  script <- call_build_script(res$build_script, "input.csv", "figure")
  run <- v3c_run_script(script, v3c_km_space())
  expect_true(is.null(run$status) || run$status == 0, info = run$out)
})
