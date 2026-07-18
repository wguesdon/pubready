# Source every library and recipe file. Entry scripts source this first.
if (!exists("CORE_R")) {
  CORE_R <- file.path(Sys.getenv("PUBPLOT_CORE", "/opt/pubplot/core"), "r")
}
for (f in sort(list.files(file.path(CORE_R, "lib"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}
for (f in sort(list.files(file.path(CORE_R, "recipes"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}
