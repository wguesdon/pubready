# Data provenance log

- Original file: `cytokine_panel.csv`
- Stored in bundle as: `input_cytokine_panel.csv`
- md5: `f77a1f672f5ed66b30e6a03fd7e90698`
- sha256: `4ccf271257662e1dee4339f1d69d6d74ce5bd3f47078a5fde1a0aed9bdca004d`
- Analyzed rows: 432
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `feature`: object
- `group`: object
- `value`: float64

## Cleaning steps

- Read long data: 9 feature(s), 3 group(s).
