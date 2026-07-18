# Reading tidy data and checksumming inputs.

# Read a tidy CSV or Excel file into a data frame.
read_tidy <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop(sprintf("input file not found: %s", path))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1))
  } else {
    as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
}

# md5 and sha256 of the exact bytes of a file.
file_checksums <- function(path) {
  list(
    md5    = unname(tools::md5sum(path)),
    sha256 = digest::digest(file = path, algo = "sha256")
  )
}
