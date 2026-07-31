required <- c(
  "tidyverse",
  "vegan",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "BiocManager"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
