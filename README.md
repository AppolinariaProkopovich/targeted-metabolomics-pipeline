# Targeted LC–MS Metabolomics Pipeline

An R workflow for processing and analysing targeted LC–MS/MS metabolomics data. The pipeline covers technical-replicate averaging, QC-based LOESS drift correction, blank and QC filtering, PCA, PERMANOVA, limma differential analysis, and publication-oriented visualisations.

## Data availability

**No study data or data-derived results are included.** The underlying data are not owned solely by the repository author and are therefore not authorised for public redistribution. See [DATA_POLICY.md](DATA_POLICY.md).

## Repository structure

```text
R/metabolomics_analysis.R       main analysis pipeline
config/config.example.R         public configuration template
data/README.md                   expected local input structure
results/README.md                generated-output policy
scripts/install_dependencies.R  package installation helper
```

## Requirements

- R 4.3 or later
- CRAN: `tidyverse`, `vegan`, `pheatmap`, `RColorBrewer`, `ggrepel`
- Bioconductor: `limma`

Install dependencies:

```bash
Rscript scripts/install_dependencies.R
```

## Local setup

1. Copy the public configuration template:

   ```bash
   cp config/config.example.R config/config.R
   ```

2. Edit `config/config.R` with local file paths, excluded QC samples, excluded biological samples, optional metabolite rescue list, and sample-name replacements.
3. Put authorised input files under `data/` or reference them by an external local path.
4. Run:

   ```bash
   Rscript R/metabolomics_analysis.R
   ```

The pipeline writes figures and tables under `results/`. Both local data and generated outputs are ignored by git.

## Expected inputs

The peak-area table is expected to contain a metabolite identifier column (default: `Molecule`) and one column per injection/sample. The injection-order file is used to map sample names to acquisition order for QC-based LOESS correction. Naming rules for Control, THz, IR, 10-minute, 45-minute, QC, and blank samples are currently encoded in the main script and should be adapted for other experimental designs.

## Reproducibility note

Study-specific manual decisions are intentionally kept outside version control in `config/config.R`. This protects restricted metadata, but means an authorised user must supply the corresponding local configuration to reproduce the exact original analysis.
