# The emitted scripts of the comparison recipes draw the geom of the figure, and
# they run. Each script runs in a temporary folder, as a user runs it from a
# bundle. The Python suite (test_emitted_scripts.py) runs the same cases.

es_spec <- function(recipe, geom) {
  spec_from_opt(list(recipe = recipe, data = "d.csv", x = "group", y = "value",
                     test = "auto", geom = geom))
}

es_data <- function(k) {
  i <- 1:8
  data.frame(group = rep(c("ctrl", "low", "high")[seq_len(k)], each = 8),
             value = c(10 + sin(i), 14 + cos(i), 18 + sin(2 * i))[seq_len(8 * k)])
}

# The ggplot2 call that only the given geom emits.
es_marker <- function(geom) {
  switch(geom, box = "geom_boxplot(aes", violin = "geom_violin(", bar = 'geom = "bar"',
         raincloud = "ggdist::stat_halfeye(")
}

# Write the input and the script into a new temporary folder, run the script
# with Rscript there, and return the folder, the exit status and the output.
es_run_script <- function(script, df) {
  dir <- tempfile("script_")
  dir.create(dir)
  utils::write.csv(df, file.path(dir, "input.csv"), row.names = FALSE)
  writeLines(script, file.path(dir, "script.R"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "script.R",
                                  stdout = TRUE, stderr = TRUE))
  list(dir = dir, status = attr(out, "status"), out = paste(out, collapse = "\n"))
}

for (geom in c("box", "violin", "bar", "raincloud")) {
  test_that(sprintf("the two-group script draws the %s geom and runs", geom), {
    df <- es_data(2)
    res <- recipe_two_group_compare(df, es_spec("two_group_compare", geom))
    script <- call_build_script(res$build_script, "input.csv", "figure")
    expect_true(any(grepl(es_marker(geom), script, fixed = TRUE)))
    run <- es_run_script(script, df)
    expect_null(run$status, info = run$out)
    expect_true(file.exists(file.path(run$dir, "figure.pdf")))
  })

  test_that(sprintf("the multi-group script draws the %s geom and runs", geom), {
    df <- es_data(3)
    res <- recipe_multi_group_compare(df, es_spec("multi_group_compare", geom))
    script <- call_build_script(res$build_script, "input.csv", "figure")
    expect_true(any(grepl(es_marker(geom), script, fixed = TRUE)))
    run <- es_run_script(script, df)
    expect_null(run$status, info = run$out)
    expect_true(file.exists(file.path(run$dir, "figure.pdf")))
  })
}

test_that("a script without points draws no points", {
  spec <- es_spec("two_group_compare", "box")
  spec$appearance$show_points <- FALSE
  res <- recipe_two_group_compare(es_data(2), spec)
  script <- call_build_script(res$build_script, "input.csv", "figure")
  expect_false(any(grepl("geom_jitter", script, fixed = TRUE)))
})
