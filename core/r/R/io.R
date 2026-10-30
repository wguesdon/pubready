# Reading tidy data and checksumming inputs.
#
# One read rule for both engines. Every cell of a CSV file is read as text, the
# spaces and tabs around each value are removed, and only an empty cell or "NA"
# reads as missing. A cell that is a number with a comma ("85,45" or "4,380")
# stops the run, because the comma can be a decimal mark or a thousands mark.
# A column becomes numeric when each present cell is a number with "." as the
# decimal mark, and logical when each present cell is TRUE or FALSE. Every other
# column stays text. The Python engine (pubready/io.py) applies the same rule.

# The cells that read as missing.
READ_NA_TOKENS <- c("", "NA")

# A number with "." as the decimal mark and no thousands mark.
READ_NUMBER_RE <- "^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?$"

# A number with a comma: a decimal comma or a thousands mark.
READ_COMMA_NUMBER_RE <- "^[+-]?[0-9]*(,[0-9]+)+([.][0-9]*)?([eE][+-]?[0-9]+)?$"

# The cells that read as logical.
READ_TRUE_TOKENS <- c("TRUE", "True", "true")
READ_FALSE_TOKENS <- c("FALSE", "False", "false")

# Join "'name' (count)" pairs for a read step.
read_counts_text <- function(counts) {
  paste(sprintf("'%s' (%d)", names(counts), counts), collapse = ", ")
}

# Apply the read rule to the text columns of a table. Returns the table and the
# read steps for the data log.
tidy_text_columns <- function(df) {
  steps <- character(0)
  old_names <- names(df)
  names(df) <- trimws(old_names, whitespace = "[ \t]")
  n_renamed <- sum(names(df) != old_names)
  if (n_renamed > 0) {
    steps <- c(steps, sprintf("Removed the spaces around %d column name(s).", n_renamed))
  }
  trimmed <- integer(0)
  missing <- integer(0)
  for (cn in names(df)) {
    v <- df[[cn]]
    if (!is.character(v)) next
    present <- !is.na(v)
    tv <- trimws(v, whitespace = "[ \t]")
    n_trim <- sum(present & tv != v)
    if (n_trim > 0) trimmed[cn] <- n_trim
    is_na_token <- present & tv %in% READ_NA_TOKENS
    if (any(is_na_token)) missing[cn] <- sum(is_na_token)
    tv[is_na_token] <- NA_character_
    present <- !is.na(tv)
    comma <- present & grepl(READ_COMMA_NUMBER_RE, tv, perl = TRUE)
    if (any(comma)) {
      stop(sprintf(paste0(
        "column '%s' holds '%s', a number with a comma. pubready reads '.' as ",
        "the decimal mark and no thousands mark. Replace each decimal comma with ",
        "'.' and remove each thousands mark, then run again."),
        cn, tv[comma][1]), call. = FALSE)
    }
    pv <- tv[present]
    if (all(grepl(READ_NUMBER_RE, pv, perl = TRUE))) {
      df[[cn]] <- as.numeric(tv)
    } else if (all(pv %in% c(READ_TRUE_TOKENS, READ_FALSE_TOKENS))) {
      df[[cn]] <- ifelse(present, tv %in% READ_TRUE_TOKENS, NA)
    } else {
      df[[cn]] <- tv
    }
  }
  if (length(trimmed) > 0) {
    steps <- c(steps, sprintf("Removed the spaces around %d value(s): %s.",
                              sum(trimmed), read_counts_text(trimmed)))
  }
  if (length(missing) > 0) {
    steps <- c(steps, sprintf("Read %d cell(s) that hold \"NA\" or nothing as missing: %s.",
                              sum(missing), read_counts_text(missing)))
  }
  list(df = df, steps = steps)
}

# Read a tidy CSV or Excel file with the read rule. Returns the table and the
# read steps for the data log.
read_tidy_logged <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop(sprintf("input file not found: %s", path), call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    # An Excel cell keeps the type that the sheet gives it. The text cells get
    # the read rule.
    raw <- readxl::read_excel(path, sheet = sheet %||% 1, na = "", trim_ws = FALSE,
                              guess_max = 100000)
  } else {
    raw <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
                           na = character(0), trim_ws = FALSE,
                           show_col_types = FALSE, progress = FALSE)
  }
  tidy_text_columns(as.data.frame(raw, stringsAsFactors = FALSE))
}

# Read a tidy CSV or Excel file into a data frame, with the read rule.
read_tidy <- function(path, sheet = NULL) {
  read_tidy_logged(path, sheet)$df
}

# md5 and sha256 of the exact bytes of a file.
file_checksums <- function(path) {
  list(
    md5    = unname(tools::md5sum(path)),
    sha256 = digest::digest(file = path, algo = "sha256")
  )
}
