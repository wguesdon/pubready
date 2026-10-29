# Data provenance log

- Original file: `sankey.csv`
- Stored in bundle as: `input_response.csv`
- md5: `c904039b35f8f4e54992ca687d27384a`
- sha256: `d6b6ba5c1c315133158bff3441fd72223f9c75710a9724354efae2d2be283e5d`
- Analyzed rows: 188
- Columns: 3

## Column types (analyzed data)

- `id`: object
- `visit`: object
- `state`: object

## Cleaning steps

- Kept 188 row(s) of 188: 80 patient(s), 3 visit(s) (8, 16, 24), 4 state(s) (complete response, partial response, stable disease, progressive disease).
- 31 patient(s) are missing at one visit or more and take no part in the transitions into or out of that visit.
