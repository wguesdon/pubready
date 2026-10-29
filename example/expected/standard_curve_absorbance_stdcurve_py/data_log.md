# Data provenance log

- Original file: `standard_curve.csv`
- Stored in bundle as: `input_absorbance.csv`
- md5: `ae4cacec29ae88addba4e46fa7189b77`
- sha256: `7de467ae8ae5d96b391e262ea55a9c73945b22fb558fcd2dbfe0ed711f16a2f8`
- Analyzed rows: 36
- Columns: 4

## Column types (analyzed data)

- `id`: object
- `role`: object
- `concentration`: float64
- `signal`: float64

## Cleaning steps

- Kept 36 of 36 row(s): 16 standard(s) and 20 unknown(s).
- Left 2 standard(s) at a concentration of zero out of the four-parameter logistic fit, because it is fitted on log10 concentration.
