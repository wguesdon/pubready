# Data provenance log

- Original file: `twoway_response.csv`
- Stored in bundle as: `input_response.csv`
- md5: `dba50f83c02a1cd007bc013dd062a9b0`
- sha256: `e7d2443929aced46222aa951eb52064fade5b9ab5ceed25f6ebf447676bc8ff4`
- Analyzed rows: 48
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `genotype`: factor
- `treatment`: factor
- `response`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
