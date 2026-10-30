# Data provenance log

- Original file: `de_basemean.csv`
- Stored in bundle as: `input_de_basemean.csv`
- md5: `ceecab1a28ae49f1e26d782b0b37af8a`
- sha256: `7a8554f9eb4e6e1dc0b66c93eb93022238d9fb7b936c1ca92b10688d6cac0e07`
- Analyzed rows: 6000
- Columns: 4

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- Read 76 cell(s) that hold "NA" or nothing as missing: 'padj' (76).

## Column types (analyzed data)

- `gene`: object
- `A`: float64
- `log2FoldChange`: float64
- `padj`: float64

## Cleaning steps

- Kept 76 row(s) with a missing p as not significant.
- Plotted log10(baseMean) as the mean expression A.
