# Data provenance log

- Original file: `pca_samples.csv`
- Stored in bundle as: `input_pca_samples.csv`
- md5: `9e74a7175197a7384718838101886b09`
- sha256: `d9d6753bbd48ee1f41dbaae9da3d9e4ee5bb7f772840824cade20cba014cd984`
- Analyzed rows: 36
- Columns: 42

## Read rule

Each cell was read as text, and the spaces around each value were removed. Only an empty cell or "NA" was read as missing. A column was read as numeric when each present cell is a number with "." as the decimal mark, and as logical when each present cell is TRUE or FALSE.

- No value changed at read.

## Column types (analyzed data)

- `sample`: character
- `group`: character
- `feature01`: numeric
- `feature02`: numeric
- `feature03`: numeric
- `feature04`: numeric
- `feature05`: numeric
- `feature06`: numeric
- `feature07`: numeric
- `feature08`: numeric
- `feature09`: numeric
- `feature10`: numeric
- `feature11`: numeric
- `feature12`: numeric
- `feature13`: numeric
- `feature14`: numeric
- `feature15`: numeric
- `feature16`: numeric
- `feature17`: numeric
- `feature18`: numeric
- `feature19`: numeric
- `feature20`: numeric
- `feature21`: numeric
- `feature22`: numeric
- `feature23`: numeric
- `feature24`: numeric
- `feature25`: numeric
- `feature26`: numeric
- `feature27`: numeric
- `feature28`: numeric
- `feature29`: numeric
- `feature30`: numeric
- `feature31`: numeric
- `feature32`: numeric
- `feature33`: numeric
- `feature34`: numeric
- `feature35`: numeric
- `feature36`: numeric
- `feature37`: numeric
- `feature38`: numeric
- `feature39`: numeric
- `feature40`: numeric

## Cleaning steps

- Used 36 samples and 40 numeric features: 'feature01', 'feature02', 'feature03', 'feature04', 'feature05', 'feature06', 'feature07', 'feature08', 'feature09', 'feature10', 'feature11', 'feature12', 'feature13', 'feature14', 'feature15', 'feature16', 'feature17', 'feature18', 'feature19', 'feature20', and 20 more.
- Left out the text column(s) 'sample'.
