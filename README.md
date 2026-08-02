# Targeted LC–MS Metabolomics Pipeline

This repository contains an R program for analysing measurements of small molecules called metabolites. Metabolites are substances produced or used by cells, such as amino acids, sugars, and products of energy metabolism. Their measured amounts provide a snapshot of how cells respond to different experimental conditions.

The measurements come from targeted liquid chromatography–tandem mass spectrometry (LC–MS/MS). In simple terms, the instrument separates a sample into its chemical components and records a signal for each metabolite of interest. A larger signal usually indicates a larger amount of that metabolite, although the raw signal must be cleaned and normalised before biological groups can be compared.

## What the code does

The pipeline turns a table of raw instrument signals into quality-control plots, statistical comparisons, and summary figures. It performs the following steps:

1. **Loads the measurements.** Rows represent metabolites and columns represent sample injections.
2. **Combines technical replicates.** When the same biological sample was measured more than once, the signals are averaged so that it contributes only once to the analysis.
3. **Checks the unprocessed data.** PCA plots give an overview of which samples look similar and help reveal unusual samples or technical patterns.
4. **Corrects instrument drift.** Repeated quality-control (QC) samples show how instrument response changed during the measurement sequence. A smooth LOESS curve is fitted to this change and used to correct the other samples.
5. **Removes unreliable metabolites.** Signals that are not clearly above the blank background or that vary too much among QC samples are excluded.
6. **Compares the experimental groups.** PCA visualises the overall patterns, PERMANOVA tests whether the complete metabolic profiles differ between groups, and `limma` tests each metabolite separately.
7. **Produces results.** The program writes PDF figures and CSV tables, including PCA plots, drift diagnostics, pairwise statistical comparisons, a heatmap, and summary plots of changed metabolites.

The current sample-name rules compare control, THz-treated, and IR-treated samples measured after 10 or 45 minutes. These rules describe this particular experimental design and must be changed before applying the code to a differently named dataset.

## Key terms

- **Biological sample:** material from one experimental unit whose metabolic profile is being studied.
- **Technical replicate:** a repeated instrument measurement of the same biological sample.
- **QC sample:** a common reference sample measured repeatedly to monitor instrument stability.
- **Blank:** a sample without biological material, used to estimate background signal.
- **PCA:** a visual summary that places samples with similar overall metabolite patterns close together.
- **FDR:** a multiple-testing correction used to limit false positive findings when many metabolites are tested.

## Data availability

**No study data or data-derived results are included or shared.** The underlying data are not owned solely by the repository author and are therefore not authorised for public redistribution. The code documents the analysis workflow, but running or reproducing the study-specific analysis requires separately authorised access to the inputs and private configuration. See [DATA_POLICY.md](DATA_POLICY.md).

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
