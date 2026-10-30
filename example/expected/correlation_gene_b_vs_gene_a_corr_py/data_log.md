# Data provenance log

- Original file: `correlation_xy.csv`
- Stored in bundle as: `input_gene_b_vs_gene_a.csv`
- md5: `b4e0baa020fa210c4d2c64be1569bb31`
- sha256: `b7de65016a54e3a4d34983a3ea4be3515a1b0c6387dd7ffea93aab56246f7c34`
- Analyzed rows: 60
- Columns: 2

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `gene_a`: float64
- `gene_b`: float64

## Cleaning steps

- No cleaning applied; input used as-is.
