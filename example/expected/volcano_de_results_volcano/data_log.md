# Data provenance log

- Original file: `de_results.csv`
- Stored in bundle as: `input_de_results.csv`
- md5: `cccbb84090f24cd1610407c75dcadee2`
- sha256: `487ba326154d5adb67175907eefabc4055bc7a610511a1175ca4cf5ac70248ae`
- Analyzed rows: 1200
- Columns: 3

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `gene`: character
- `log2FoldChange`: numeric
- `padj`: numeric

## Cleaning steps

- No cleaning applied; input used as-is.
