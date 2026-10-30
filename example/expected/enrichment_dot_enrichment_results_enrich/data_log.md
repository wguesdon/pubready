# Data provenance log

- Original file: `enrichment_results.csv`
- Stored in bundle as: `input_enrichment_results.csv`
- md5: `7292703ccf5170ec7811de3a431c33ee`
- sha256: `3265457a2c00ffc902d1848a14b7ff6de77284a6cf42715c6d61b56d779178e8`
- Analyzed rows: 20
- Columns: 4

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `term`: character
- `x`: numeric
- `p_adjust`: numeric
- `count`: numeric

## Cleaning steps

- No cleaning applied; enrichment table used as-is.
