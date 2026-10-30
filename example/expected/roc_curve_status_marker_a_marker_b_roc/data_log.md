# Data provenance log

- Original file: `biomarker_roc.csv`
- Stored in bundle as: `input_status_marker_a_marker_b.csv`
- md5: `f6d3a1a52479ba12f8a4ee2644beb33d`
- sha256: `093cd3826f05f33eb5db7e8f7d1b571a2e2c02a398727171a5702f39ecce3fa7`
- Analyzed rows: 140
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `outcome`: character
- `marker_a`: numeric
- `marker_b`: numeric

## Cleaning steps

- Kept 140 of 140 row(s): 60 case(s) (disease) and 80 control(s) (healthy).
