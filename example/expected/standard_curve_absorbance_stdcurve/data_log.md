# Data provenance log

- Original file: `standard_curve.csv`
- Stored in bundle as: `input_absorbance.csv`
- md5: `ae4cacec29ae88addba4e46fa7189b77`
- sha256: `7de467ae8ae5d96b391e262ea55a9c73945b22fb558fcd2dbfe0ed711f16a2f8`
- Analyzed rows: 36
- Columns: 4

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- Read 20 cell(s) that hold "NA" or nothing as missing: 'conc_ng_ml' (20).

## Column types (analyzed data)

- `id`: character
- `role`: character
- `concentration`: numeric
- `signal`: numeric

## Cleaning steps

- Kept 36 of 36 row(s): 16 standard(s) and 20 unknown(s).
- Left 2 standard(s) at a concentration of zero out of the four-parameter logistic fit, because it is fitted on log10 concentration.
