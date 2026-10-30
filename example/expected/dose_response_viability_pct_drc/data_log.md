# Data provenance log

- Original file: `dose_response.csv`
- Stored in bundle as: `input_viability_pct.csv`
- md5: `065caf43fa8bffd4422569a6a1bbece0`
- sha256: `c6b98694438d373a93c7934f5421def4892b3a726d40429bcf0025e87e231590`
- Analyzed rows: 72
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `group`: character
- `dose`: numeric
- `response`: numeric

## Cleaning steps

- Kept 72 of 72 row(s).
