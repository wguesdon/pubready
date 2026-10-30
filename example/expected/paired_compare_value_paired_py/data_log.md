# Data provenance log

- Original file: `paired_response.csv`
- Stored in bundle as: `input_value.csv`
- md5: `e995b32f1de8c3186d6e6c6f0d010a9f`
- sha256: `8bfc166ad40f6287fcf50fe3c144a4ccce89d4a501007e5e47dd183fc2df2b0f`
- Analyzed rows: 22
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `subject`: object
- `v1`: float64
- `v2`: float64

## Cleaning steps

- Ordered the conditions as 'pre', 'post', by the order that --levels gives.
- Kept 22 complete pair(s) across the two conditions.
