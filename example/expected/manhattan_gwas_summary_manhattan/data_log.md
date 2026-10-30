# Data provenance log

- Original file: `gwas_summary.csv`
- Stored in bundle as: `input_gwas_summary.csv`
- md5: `cbcaf85e088a7d1bb7bb44fd15648588`
- sha256: `40a04af482fa40bdc75f42613110744366fc3a2eb9cb5e029968fe29869c58a4`
- Analyzed rows: 20000
- Columns: 4

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `chr`: character
- `pos`: numeric
- `p`: numeric
- `snp`: character

## Cleaning steps

- Read the chromosome names without their 'chr' prefix.
