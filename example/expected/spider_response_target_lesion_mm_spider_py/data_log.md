# Data provenance log

- Original file: `spider_response.csv`
- Stored in bundle as: `input_target_lesion_mm.csv`
- md5: `d3ceaf1a40e01ce2639a76cda972d28e`
- sha256: `97ca47b3a60d61efa8b26f77be4c7169722126bb8198f60c17acef4d6f04d23c`
- Analyzed rows: 94
- Columns: 5

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `patient`: object
- `week`: int64
- `target_lesion_mm`: float64
- `arm`: object
- `change`: float64

## Cleaning steps

- Kept 20 patient(s) across 94 measurement(s).
