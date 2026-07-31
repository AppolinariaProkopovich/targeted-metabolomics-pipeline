# Public example configuration.
# Copy this file to config/config.R and adapt it locally.
# config/config.R is ignored by git and must not contain share-restricted data.

data_file <- "data/raw_peak_areas.csv"
batch_file <- "data/injection_order.csv"
raw_skip <- 1
metabolite_id_column <- "Molecule"

figures_dir <- "results/figures"
tables_dir <- "results/tables"

# Study-specific exclusions and manual decisions belong in config/config.R.
excluded_qc_samples <- character(0)
excluded_samples <- character(0)
rescue_metabolites <- character(0)

# Optional named replacements used to reconcile injection-order sample names.
# Example: c("name in batch file" = "name in peak-area table")
batch_name_replacements <- character(0)
