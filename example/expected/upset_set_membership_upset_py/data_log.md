# Data provenance log

- Original file: `set_membership.csv`
- Stored in bundle as: `input_set_membership.csv`
- md5: `62b68d6ca9159eb4de4ff4c2d35ded26`
- sha256: `5aabbab60c65f2af913bb91fd373a32c42ec6ded5459b606e3f5d099b1c70528`
- Analyzed rows: 72
- Columns: 5

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `element`: object
- `apoptosis`: bool
- `inflammation`: bool
- `metabolism`: bool
- `cell_cycle`: bool

## Cleaning steps

- Dropped 8 element(s) that belonged to no set (80 -> 72).
