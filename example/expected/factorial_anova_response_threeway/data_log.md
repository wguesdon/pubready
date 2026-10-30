# Data provenance log

- Original file: `threeway_response.csv`
- Stored in bundle as: `input_response.csv`
- md5: `b1ab1ae1731e3311d49a8727536cd99f`
- sha256: `b2267d6af2d2f49bbf1e1e75507c56dd8257334c80639c781f35701ca9ea67e8`
- Analyzed rows: 48
- Columns: 4

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `genotype`: factor
- `treatment`: factor
- `sex`: factor
- `response`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
