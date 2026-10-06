# Data provenance log

- Original file: `meta_trials.csv`
- Stored in bundle as: `input_hazard_ratio.csv`
- md5: `24320e4211da4db9bf76c9a1fedeb1e9`
- sha256: `01b10b5db2fbcdca95bb708eb109ab543a00c92557e91fbc44d75ad98d46acc6`
- Analyzed rows: 9
- Columns: 6

## Column types (analyzed data)

- `label`: character
- `estimate`: numeric
- `ci_low`: numeric
- `ci_high`: numeric
- `yi`: numeric
- `se`: numeric

## Cleaning steps

- Pooled 9 of 9 study row(s) as a ratio measure on the log scale.
- Derived each standard error from the 95% CI as (upper - lower) / (2 x 1.96) on the analysis scale.
