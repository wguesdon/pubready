# Data provenance log

- Original file: `umap_clusters.csv`
- Stored in bundle as: `input_umap_clusters.csv`
- md5: `76d61b415aeb6c7f7cd1b8d152260faf`
- sha256: `93e1dfb9a3457374ec2a570e17defd873f383468d2dd3ff0a70fb8f58f51b68f`
- Analyzed rows: 2590
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `UMAP_1`: numeric
- `UMAP_2`: numeric
- `cell_type`: character

## Cleaning steps

- No cleaning applied; input used as-is.
