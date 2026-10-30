# Data provenance log

- Original file: `cytokine_pg_ml.csv`
- Stored in bundle as: `input_concentration.csv`
- md5: `71253c2a457715e5eb78ba8b5da99b7b`
- sha256: `b4591b64add603dd5d4d265e4d7ddb417b65359ca5c76a169152e8d01041413c`
- Analyzed rows: 24
- Columns: 2

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `group`: factor
- `concentration`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
