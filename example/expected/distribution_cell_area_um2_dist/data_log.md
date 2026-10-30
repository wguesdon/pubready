# Data provenance log

- Original file: `distribution.csv`
- Stored in bundle as: `input_cell_area_um2.csv`
- md5: `9a6485554175f36ec3bd2b8365cc3665`
- sha256: `5465c7bdc86c61c1bbdaa5bf4d211be7062d1eadd2b12df0a28b2cecba8052bb`
- Analyzed rows: 240
- Columns: 2

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `group`: character
- `value`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
