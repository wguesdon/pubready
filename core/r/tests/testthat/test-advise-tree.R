# The R engine checks the same tree and fixtures as the Python engine.
advise_tree <- function() {
  yaml::read_yaml("/opt/pubready/reference/decision_tree.yaml")
}

advise_fixture <- function(recipe, method) {
  fixtures <- jsonlite::fromJSON("/opt/pubready/core/advise_test_fixtures.json", simplifyVector = FALSE)
  fixture <- fixtures[[recipe]]
  base <- file.path(Sys.getenv("PUBREADY_BENCHMARKS"), "data")
  flags <- fixture$flags
  if (!is.null(flags$annotation)) flags$annotation <- file.path(base, flags$annotation)
  opt <- utils::modifyList(list(recipe = recipe, data = "fixture.csv", test = method), flags)
  list(data = read_tidy(file.path(base, paste0(fixture$dataset, ".csv"))), spec = spec_from_opt(opt))
}

test_that("every advise leaf exists in R and its advertised tests report a listed outcome", {
  tree <- advise_tree()
  tested <- character()
  for (node in tree$nodes) {
    if (is.null(node$leaf)) next
    leaf <- node$leaf
    recipe <- leaf$recipe
    expect_true(exists(paste0("recipe_", recipe), mode = "function"))
    if (recipe %in% tested || is.null(leaf$test)) next
    tested <- c(tested, recipe)
    fn <- get(paste0("recipe_", recipe), mode = "function")
    outcomes <- vapply(leaf$test$outcomes, `[[`, character(1), "test")
    for (method in leaf$test$allowed) {
      fixture <- advise_fixture(recipe, method)
      result <- fn(fixture$data, fixture$spec)
      expect_true(result$test_meta$name %in% outcomes,
                  info = paste(recipe, method, result$test_meta$name))
    }
    fixture <- advise_fixture(recipe, "unknown_advise_test")
    expect_error(fn(fixture$data, fixture$spec))
  }
})
