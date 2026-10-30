# Data provenance log

- Original file: `composition.csv`
- Stored in bundle as: `input_cell_type.csv`
- md5: `4a74ef936bfc21d4bf2ccebb6e699014`
- sha256: `329d74c430d92f2d03a0ec9b91979eaa4ead49770e328224a89d2665269c7910`
- Analyzed rows: 96
- Columns: 5

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `sample`: character
- `group`: character
- `category`: character
- `count`: numeric
- `proportion`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
