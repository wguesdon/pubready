# Data provenance log

- Original file: `swimmer_patients.csv`
- Stored in bundle as: `input_months.csv`
- md5: `b8d75bab35e8ca31ee1a4ca5fc089e74`
- sha256: `94541cad1fc2ff7d0c412a3c79ef9f5b9a7fd2544a87bf360389f987e24dbe36`
- Analyzed rows: 24
- Columns: 7

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- Read 23 cell(s) that hold "NA" or nothing as missing: 'response_month' (12), 'progression_month' (11).

## Column types (analyzed data)

- `patient`: object
- `duration`: float64
- `group`: object
- `ongoing`: bool
- `response_month`: float64
- `progression_month`: float64
- `rank`: int64

## Cleaning steps

- Kept 24 of 24 patient(s).
