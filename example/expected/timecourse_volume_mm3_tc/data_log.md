# Data provenance log

- Original file: `tumor_growth.csv`
- Stored in bundle as: `input_volume_mm3.csv`
- md5: `b86e88025259e6532a144d1496106495`
- sha256: `dbf732b2b6cf9b41409b5419f1f5a21faabab658b627b18f9ac51654c2f09614`
- Analyzed rows: 161
- Columns: 4

## Column types (analyzed data)

- `id`: character
- `group`: character
- `time`: numeric
- `y`: numeric

## Cleaning steps

- Dropped 1 subject(s) without a value at every one of the 7 time points, because the repeated-measures ANOVA needs complete subjects.
- Kept 23 subject(s) at 7 time points (161 values).
