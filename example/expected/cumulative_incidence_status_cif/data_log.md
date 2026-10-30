# Data provenance log

- Original file: `cumulative_incidence.csv`
- Stored in bundle as: `input_status.csv`
- md5: `f51b7f002d22481f599b4424a27f3bdc`
- sha256: `d51f82fe2977bbed1ce2d5a6dac5172c07d0735d16d33c39715d0b67ab610486`
- Analyzed rows: 200
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `time`: numeric
- `status`: character
- `group`: factor

## Cleaning steps

- Kept 200 of 200 row(s): 75 censored, 81 with status = 1, 44 with status = 2.
