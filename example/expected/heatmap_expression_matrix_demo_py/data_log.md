# Data provenance log

- Original file: `expression_matrix.csv`
- Stored in bundle as: `input_expression_matrix.csv`
- md5: `ca9f8a5e682388268ffb561ea779a536`
- sha256: `88a6c73257c7a2c3c3f481165584f21cb60cc907778508e023fb2762472198f4`
- Analyzed rows: 30
- Columns: 21

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `gene`: object
- `ctrl_01`: float64
- `ctrl_02`: float64
- `ctrl_03`: float64
- `ctrl_04`: float64
- `ctrl_05`: float64
- `ctrl_06`: float64
- `ctrl_07`: float64
- `ctrl_08`: float64
- `ctrl_09`: float64
- `ctrl_10`: float64
- `trt_01`: float64
- `trt_02`: float64
- `trt_03`: float64
- `trt_04`: float64
- `trt_05`: float64
- `trt_06`: float64
- `trt_07`: float64
- `trt_08`: float64
- `trt_09`: float64
- `trt_10`: float64

## Cleaning steps

- No cleaning applied; input used as-is.
