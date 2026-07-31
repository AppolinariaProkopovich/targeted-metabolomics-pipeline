# Local data directory

The analysis data are **not distributed with this repository**. They are excluded because the repository author does not solely own the underlying study data.

Place the following files here only on your local machine:

- `raw_peak_areas.csv` — targeted peak-area matrix. The first metadata row is skipped by default; the metabolite identifier column is expected to be named `Molecule`; remaining columns are samples.
- `injection_order.csv` — injection sequence used for QC-based drift correction.

Adapt paths and study-specific settings in `config/config.R`. Do not commit that file or any data-derived tables.
