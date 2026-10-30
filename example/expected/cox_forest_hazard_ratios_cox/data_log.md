# Data provenance log

- Original file: `survival_trial.csv`
- Stored in bundle as: `input_hazard_ratios.csv`
- md5: `6f53a662fffe08afdcf72fe041adaf29`
- sha256: `17dccfaf72c84d69607aa73efad89ed98d9af94e72ddd285bc864fa65346b75e`
- Analyzed rows: 120
- Columns: 6

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `time`: numeric
- `event`: integer
- `arm`: factor
- `age`: numeric
- `sex`: factor
- `stage`: factor

## Cleaning steps

- No cleaning applied; input used as-is.
