# Data provenance log

- Original file: `calibration.csv`
- Stored in bundle as: `input_event_1y_risk_clinical_risk_extended.csv`
- md5: `4d7ddf64a4e0eb79c7d2003021d46c53`
- sha256: `c6d02e1cfda0e27a2779f352e23b00abe6c739ee1b810ad3e9113fdcdc55e34b`
- Analyzed rows: 500
- Columns: 3

## Column types (analyzed data)

- `outcome`: object
- `risk_clinical`: float64
- `risk_extended`: float64

## Cleaning steps

- Kept 500 of 500 row(s): 130 event(s) (yes) and 370 non-event(s) (no).
- Predictions were clipped to [1e-08, 1 - 1e-08] before the logit transform.
