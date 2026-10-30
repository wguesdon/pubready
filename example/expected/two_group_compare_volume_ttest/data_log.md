# Data provenance log

- Original file: `tumor_volume.csv`
- Stored in bundle as: `input_volume.csv`
- md5: `d61ed3cc297d2fa63d25dcab55991645`
- sha256: `7560f561606eed71ef05c9a9b5b82b36a6f830c3bd0f3586a9a905534cc4a355`
- Analyzed rows: 24
- Columns: 2

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `group`: factor
- `volume`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
