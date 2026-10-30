# Data provenance log

- Original file: `gene_expression.csv`
- Stored in bundle as: `input_expression.csv`
- md5: `55d3ddd6425019751d19883d1fe81771`
- sha256: `0d5ab7e477ae444212d4b0eebf041038d659bee75d425c1cd51401dcd73d034b`
- Analyzed rows: 40
- Columns: 2

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `genotype`: factor
- `expression`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
