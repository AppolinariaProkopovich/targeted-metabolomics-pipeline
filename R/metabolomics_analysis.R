#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# Metabolomics data processing and statistical analysis
# ══════════════════════════════════════════════════════════════════════════════
#
# Description:
#   Reproducible processing pipeline for targeted LC-MS/MS metabolomics data.
#   The repository contains code only; study data and generated results are not
#   distributed here.
#   Steps: QC-based LOESS drift correction → blank filtering → QC filtering
#          → PCA → PERMANOVA → limma differential analysis.
#
# Input:
#   - data/raw_peak_areas.csv  (raw peak areas; local, not tracked)
#   - data/injection_order.csv  (injection order; local, not tracked)
#   - config/config.R          (local study settings; not tracked)
#
# Output:
#   - results/figures/*.pdf
#   - results/tables/*.csv
#
# Requirements:
#   R >= 4.3.0
#   Packages: tidyverse, vegan, limma, pheatmap, RColorBrewer, ggrepel
#
# Usage:
#   Rscript metabolomics_analysis.R
#
# ══════════════════════════════════════════════════════════════════════════════

# ── Local configuration -------------------------------------------------------
# Copy config/config.example.R to config/config.R and edit it for your dataset.
# config/config.R is intentionally ignored by git.

config_file <- Sys.getenv("METABOLOMICS_CONFIG", "config/config.R")
if (file.exists(config_file)) {
  source(config_file, local = FALSE)
} else {
  warning(
    "Local config not found: ", config_file,
    ". Using public defaults from config/config.example.R."
  )
  source("config/config.example.R", local = FALSE)
}

library(tidyverse)
library(vegan)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

set.seed(42)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

# ── Publication theme for all plots ──────────────────────────────────────────

theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 2, color = "black"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                      hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0),
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

# ── Color palette ────────────────────────────────────────────────────────────

group_colors <- c(
  "Blank"          = "black",
  "QC"             = "forestgreen",
  "Control_10min"  = "gray50",
  "Control_45min"  = "gray70",
  "THz_10min"      = "#E31A1C",
  "THz_45min"      = "#8B0000",
  "IR_10min"       = "#1F78B4",
  "IR_45min"       = "#08306B"
)

treatment_colors <- c(Control = "grey60", THz = "#E31A1C", IR = "#1F78B4")
time_colors      <- c("10min" = "#FDB863", "45min" = "#5E3C99")


# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

raw <- read.csv(data_file,
                skip = raw_skip, header = TRUE, check.names = FALSE)

if (!metabolite_id_column %in% colnames(raw)) {
  stop("Missing metabolite identifier column: ", metabolite_id_column)
}
rownames(raw) <- raw[[metabolite_id_column]]
raw[[metabolite_id_column]] <- NULL
colnames(raw)  <- trimws(colnames(raw))

cat("Raw data dimensions:", nrow(raw), "metabolites x", ncol(raw), "samples\n")


# ── Group assignment function ────────────────────────────────────────────────

assign_group <- function(name) {
  case_when(
    str_detect(name, "^blank")           ~ "Blank",
    str_detect(name, "^QC")              ~ "QC",
    str_detect(name, "10-min-THz")       ~ "THz_10min",
    str_detect(name, "45-min-THz")       ~ "THz_45min",
    str_detect(name, "10-min-IR")        ~ "IR_10min",
    str_detect(name, "45-min-IR")        ~ "IR_45min",
    str_detect(name, "10-min-control")   ~ "Control_10min",
    str_detect(name, "45.?min-control")  ~ "Control_45min",
    TRUE ~ "Other"
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. TECHNICAL REPLICATE AVERAGING
# ══════════════════════════════════════════════════════════════════════════════

sample_names <- trimws(colnames(raw))
groups_all   <- assign_group(sample_names)

is_bio <- groups_all %in% c("THz_10min", "THz_45min",
                             "IR_10min", "IR_45min",
                             "Control_10min", "Control_45min")

bio_raw   <- raw[, is_bio]
other_raw <- raw[, !is_bio]

bio_names <- trimws(colnames(bio_raw))
bio_id    <- str_replace(bio_names, "\\(\\d\\)\\s*$", "") %>% trimws()

avg_bio <- sapply(unique(bio_id), function(id) {
  cols <- which(bio_id == id)
  if (length(cols) == 1) return(as.numeric(bio_raw[, cols]))
  rowMeans(bio_raw[, cols], na.rm = TRUE)
})
avg_bio <- as.data.frame(avg_bio)
rownames(avg_bio) <- rownames(raw)

avg_data <- cbind(other_raw, avg_bio)

cat("After averaging tech replicates:", ncol(avg_data), "columns\n")
print(table(assign_group(colnames(avg_data))))


# ══════════════════════════════════════════════════════════════════════════════
# 3. PCA ON RAW DATA (diagnostic)
# ══════════════════════════════════════════════════════════════════════════════

min_nonzero <- min(avg_data[avg_data > 0], na.rm = TRUE)
avg_data[avg_data == 0] <- min_nonzero / 2
log2_data <- log2(avg_data)

pca_res <- prcomp(t(log2_data), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca_res)$importance[2, 1:5], 1)

sample_groups <- assign_group(colnames(avg_data))
names(sample_groups) <- colnames(avg_data)

pca_df <- as.data.frame(pca_res$x[, 1:5]) %>%
  mutate(Sample = rownames(pca_res$x),
         Group  = sample_groups[Sample])

p_raw_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(title = "PCA — raw data (all samples)",
       x = paste0("PC1 (", var_pct[1], "%)"),
       y = paste0("PC2 (", var_pct[2], "%)")) +
  theme_pub()

ggsave(file.path(figures_dir, "01_pca_raw.pdf"), p_raw_pca, width = 10, height = 8)


# ── PCA on biological samples only (pre-correction baseline) ─────────────────

bio_idx    <- sample_groups %in% c("THz_10min", "THz_45min",
                                    "IR_10min", "IR_45min",
                                    "Control_10min", "Control_45min")
log2_bio   <- log2_data[, bio_idx]
bio_groups <- sample_groups[bio_idx]

pca_bio    <- prcomp(t(log2_bio), center = TRUE, scale. = TRUE)
var_pct_bio <- round(100 * summary(pca_bio)$importance[2, 1:5], 1)

pca_bio_df <- as.data.frame(pca_bio$x[, 1:5]) %>%
  mutate(Sample = rownames(pca_bio$x),
         Group  = bio_groups[Sample])

p_bio_raw <- ggplot(pca_bio_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.7) +
  labs(title = "PCA — biological samples only (before correction)",
       x = paste0("PC1 (", var_pct_bio[1], "%)"),
       y = paste0("PC2 (", var_pct_bio[2], "%)")) +
  theme_pub()

ggsave(file.path(figures_dir, "01b_pca_bio_raw.pdf"), p_bio_raw, width = 10, height = 8)


# ══════════════════════════════════════════════════════════════════════════════
# 4. INJECTION ORDER MAPPING
# ══════════════════════════════════════════════════════════════════════════════

batch_raw <- read.csv(batch_file,
                      header = FALSE, stringsAsFactors = FALSE)

batch_df <- data.frame(batch_name = trimws(batch_raw$V1)) %>%
  mutate(row_num = row_number()) %>%
  filter(batch_name != "",
         !str_detect(batch_name, "^ACN"),
         !str_detect(batch_name, "^THz-C5")) %>%
  filter(str_detect(batch_name, "HILIC")) %>%
  mutate(injection_num = row_number(),
         sample_key = str_replace(batch_name, "^C5-", ""),
         sample_key = str_replace(sample_key, "\\s+\\d{2}-\\d{2}-\\d{2}.*$", ""),
         sample_key = str_replace(sample_key, "^THz-", ""),
         sample_key = trimws(sample_key))

inj_order <- setNames(batch_df$injection_num, batch_df$sample_key)

# Optional naming fixes supplied in the local config
if (length(batch_name_replacements) > 0) {
  for (old_name in names(batch_name_replacements)) {
    names(inj_order)[names(inj_order) == old_name] <-
      unname(batch_name_replacements[[old_name]])
  }
}

results_cols <- colnames(raw)
matched <- results_cols %in% names(inj_order)
cat("Injection order matching:", sum(matched), "/", length(results_cols), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# 5. LOESS BATCH CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

# QC columns for LOESS fitting (excluding outlier QC12)
qc_cols    <- grep("^QC", colnames(raw), value = TRUE)
qc_for_fit <- qc_cols[!qc_cols %in% excluded_qc_samples]

# ── Visualize drift BEFORE correction ────────────────────────────────────────

random_mets <- sample(rownames(raw), 12)

qc_long_before <- raw[, qc_for_fit] %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(-Metabolite, names_to = "QC_sample", values_to = "Value") %>%
  mutate(Injection = inj_order[QC_sample])

p_drift_before <- qc_long_before %>%
  filter(Metabolite %in% random_mets) %>%
  ggplot(aes(x = Injection, y = Value)) +
  geom_point(color = "forestgreen", size = 2.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red",
              span = 0.75, linewidth = 1) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 4) +
  labs(title = "QC drift BEFORE correction (12 random metabolites)",
       x = "Injection order", y = "Raw signal") +
  theme_pub(base_size = 12)

ggsave(file.path(figures_dir, "02_qc_drift_before.pdf"), p_drift_before,
       width = 14, height = 10)


# ── Apply LOESS correction ──────────────────────────────────────────────────

cols_to_correct <- colnames(raw)[!grepl("^blank", colnames(raw))]
cols_to_correct <- cols_to_correct[cols_to_correct %in% names(inj_order)]

normalized_raw <- raw

for (met in rownames(raw)) {
  qc_vals <- as.numeric(raw[met, qc_for_fit])
  qc_injs <- inj_order[qc_for_fit]

  if (sd(qc_vals, na.rm = TRUE) == 0) next

  loess_fit <- tryCatch(
    loess(qc_vals ~ qc_injs, span = 0.75),
    error = function(e) NULL
  )
  if (is.null(loess_fit)) next

  qc_median <- median(qc_vals)

  for (col in cols_to_correct) {
    col_inj <- inj_order[col]
    if (is.na(col_inj)) next
    predicted <- predict(loess_fit, newdata = data.frame(qc_injs = col_inj))
    if (!is.na(predicted) && predicted > 0) {
      normalized_raw[met, col] <- raw[met, col] * (qc_median / predicted)
    }
  }
}

cat("LOESS correction complete\n")


# ── Evaluate correction: CV before vs after ──────────────────────────────────

cv_before <- sapply(rownames(raw), function(met) {
  vals <- as.numeric(raw[met, qc_for_fit])
  sd(vals) / mean(vals) * 100
})
cv_after <- sapply(rownames(normalized_raw), function(met) {
  vals <- as.numeric(normalized_raw[met, qc_for_fit])
  sd(vals) / mean(vals) * 100
})

cat("QC CV% (excluding QC12):\n")
cat("  BEFORE: median CV =", round(median(cv_before), 1), "%\n")
cat("  AFTER:  median CV =", round(median(cv_after), 1), "%\n")


# ── Visualize drift AFTER correction ─────────────────────────────────────────

qc_long_after <- normalized_raw[, qc_for_fit] %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(-Metabolite, names_to = "QC_sample", values_to = "Value") %>%
  mutate(Injection = inj_order[QC_sample])

p_drift_after <- qc_long_after %>%
  filter(Metabolite %in% random_mets) %>%
  ggplot(aes(x = Injection, y = Value)) +
  geom_point(color = "forestgreen", size = 2.5) +
  geom_smooth(method = "loess", se = FALSE, color = "blue",
              span = 0.75, linewidth = 1) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 4) +
  labs(title = "QC drift AFTER correction",
       subtitle = "Lines should be approximately flat",
       x = "Injection order", y = "Corrected signal") +
  theme_pub(base_size = 12)

ggsave(file.path(figures_dir, "03_qc_drift_after.pdf"), p_drift_after,
       width = 14, height = 10)


# ══════════════════════════════════════════════════════════════════════════════
# 6. AVERAGE TECHNICAL REPLICATES (on corrected data)
# ══════════════════════════════════════════════════════════════════════════════

sample_names_corr <- colnames(normalized_raw)
sample_groups_all <- assign_group(sample_names_corr)

is_bio     <- sample_groups_all %in% c("THz_10min", "THz_45min",
                                        "IR_10min", "IR_45min",
                                        "Control_10min", "Control_45min")
bio_raw_c  <- normalized_raw[, is_bio]
other_raw_c <- normalized_raw[, !is_bio]

bio_names_c <- colnames(bio_raw_c)
bio_id_c    <- str_replace(bio_names_c, "\\(\\d\\)\\s*$", "") %>% trimws()

avg_bio_c <- sapply(unique(bio_id_c), function(id) {
  cols <- which(bio_id_c == id)
  if (length(cols) == 1) return(as.numeric(bio_raw_c[, cols]))
  rowMeans(bio_raw_c[, cols], na.rm = TRUE)
})
avg_bio_c <- as.data.frame(avg_bio_c)
rownames(avg_bio_c) <- rownames(normalized_raw)

final_data <- cbind(other_raw_c, avg_bio_c)

sample_groups <- assign_group(colnames(final_data))
names(sample_groups) <- colnames(final_data)

cat("After LOESS + averaging:", ncol(final_data), "columns\n")
print(table(sample_groups))


# ── PCA after LOESS correction (bio only) ────────────────────────────────────

bio_idx_c     <- sample_groups %in% c("THz_10min", "THz_45min",
                                       "IR_10min", "IR_45min",
                                       "Control_10min", "Control_45min")
bio_corrected <- final_data[, bio_idx_c]
bio_groups_c  <- sample_groups[bio_idx_c]

min_nz_c <- min(unlist(bio_corrected)[unlist(bio_corrected) > 0])
bio_corrected[bio_corrected == 0] <- min_nz_c / 2
log2_bio_c <- log2(bio_corrected)

pca_bio_c    <- prcomp(t(log2_bio_c), center = TRUE, scale. = TRUE)
var_pct_bio_c <- round(100 * summary(pca_bio_c)$importance[2, 1:5], 1)

pca_bio_c_df <- as.data.frame(pca_bio_c$x[, 1:5]) %>%
  mutate(Sample = rownames(pca_bio_c$x),
         Group  = bio_groups_c[Sample])

p_pca_loess <- ggplot(pca_bio_c_df,
                      aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.7) +
  labs(title = "PCA after LOESS correction (bio samples only)",
       x = paste0("PC1 (", var_pct_bio_c[1], "%)"),
       y = paste0("PC2 (", var_pct_bio_c[2], "%)")) +
  theme_pub()

ggsave(file.path(figures_dir, "04_pca_after_loess.pdf"), p_pca_loess, width = 10, height = 8)


# ══════════════════════════════════════════════════════════════════════════════
# 7. BLANK FILTERING
# ══════════════════════════════════════════════════════════════════════════════

blank_cols <- grep("^blank", colnames(final_data), value = TRUE)
bio_cols   <- colnames(final_data)[sample_groups %in%
              c("THz_10min", "THz_45min", "IR_10min", "IR_45min",
                "Control_10min", "Control_45min")]

blank_mean  <- rowMeans(final_data[, blank_cols])
sample_mean <- rowMeans(final_data[, bio_cols])
snr         <- sample_mean / (blank_mean + 1e-20)

cat("\nBlank filtering (signal/blank > 3):\n")
cat("  Before:", nrow(final_data), "metabolites\n")
cat("  SNR > 3:", sum(snr > 3), "\n")
cat("  SNR <= 3 (removed):", sum(snr <= 3), "\n")

keep_blank <- snr > 3

# Optional study-specific rescue list from the local config
keep_blank[names(keep_blank) %in% rescue_metabolites] <- TRUE

data_filt <- final_data[keep_blank, ]
cat("  After (with rescue):", nrow(data_filt), "metabolites\n")


# ══════════════════════════════════════════════════════════════════════════════
# 8. QC FILTERING (CV < 30%)
# ══════════════════════════════════════════════════════════════════════════════

qc_cols_clean <- grep("^QC", colnames(data_filt), value = TRUE)
qc_cols_clean <- qc_cols_clean[!qc_cols_clean %in% excluded_qc_samples]

cv_pct <- sapply(rownames(data_filt), function(met) {
  vals <- as.numeric(data_filt[met, qc_cols_clean])
  sd(vals) / mean(vals) * 100
})

cat("\nQC filtering:\n")
cat("  CV < 20%:", sum(cv_pct < 20), "\n")
cat("  CV < 25%:", sum(cv_pct < 25), "\n")
cat("  CV < 30%:", sum(cv_pct < 30), "\n")
cat("  CV >= 30% (removed):", sum(cv_pct >= 30), "\n")

keep_cv    <- cv_pct < 30
data_clean <- data_filt[keep_cv, ]
cat("  Final metabolite count:", nrow(data_clean), "\n")

cat("\nRemoved by QC filter:\n")
cat(names(cv_pct)[cv_pct >= 30], sep = "\n")


# ══════════════════════════════════════════════════════════════════════════════
# 9. PCA ON CLEAN DATA
# ══════════════════════════════════════════════════════════════════════════════

bio_clean <- data_clean[, bio_cols]

min_nz_clean <- min(unlist(bio_clean)[unlist(bio_clean) > 0])
bio_clean[bio_clean == 0] <- min_nz_clean / 2
log2_clean <- log2(bio_clean)

meta <- data.frame(
  Sample    = bio_cols,
  Group     = sample_groups[bio_cols],
  Treatment = case_when(
    str_detect(bio_cols, "THz")     ~ "THz",
    str_detect(bio_cols, "IR")      ~ "IR",
    str_detect(bio_cols, "control") ~ "Control"
  ),
  Time = case_when(
    str_detect(bio_cols, "10-min")  ~ "10min",
    str_detect(bio_cols, "45")      ~ "45min"
  ),
  row.names = bio_cols
)

cat("\nSample metadata:\n")
print(table(meta$Treatment, meta$Time))

pca_clean    <- prcomp(t(log2_clean), center = TRUE, scale. = TRUE)
var_clean    <- round(100 * summary(pca_clean)$importance[2, 1:5], 1)

pca_clean_df <- as.data.frame(pca_clean$x[, 1:5]) %>%
  mutate(Sample = rownames(pca_clean$x),
         Group  = meta[Sample, "Group"])

# ── FIGURE 5A: PCA after full filtering pipeline ─────────────────────────────

p_pca_clean <- ggplot(pca_clean_df,
                      aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text_repel(aes(label = Sample), size = 3.5,
                  max.overlaps = 25, show.legend = FALSE,
                  segment.color = "grey60", segment.size = 0.3,
                  fontface = "bold") +
  scale_color_manual(values = group_colors,
                     name = "Group") +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.8) +
  labs(x = paste0("PC1 (", var_clean[1], "%)"),
       y = paste0("PC2 (", var_clean[2], "%)")) +
  theme_pub(base_size = 15) +
  theme(
    legend.key.size = unit(0.6, "cm"),
    plot.margin     = margin(10, 15, 10, 10)
  )

ggsave(file.path(figures_dir, "05_pca_clean.pdf"), p_pca_clean,
       width = 10, height = 8, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 10. PERMANOVA
# ══════════════════════════════════════════════════════════════════════════════

dist_mat <- dist(t(log2_clean), method = "euclidean")

perm_res <- adonis2(dist_mat ~ Treatment * Time,
                    data = meta,
                    permutations = 9999,
                    method = "euclidean")

cat("\n══ PERMANOVA: Treatment * Time (all samples) ══\n")
print(perm_res)

# ── Pairwise comparisons (all samples) ───────────────────────────────────────

for (t_point in c("10min", "45min")) {
  cat("\n── Pairwise PERMANOVA:", t_point, "──\n")
  idx <- meta$Time == t_point

  for (treat in c("THz", "IR")) {
    idx2 <- idx & meta$Treatment %in% c(treat, "Control")
    d <- dist(t(log2_clean[, idx2]), method = "euclidean")
    m <- meta[idx2, ]
    res <- adonis2(d ~ Treatment, data = m, permutations = 9999)
    cat(treat, "vs Control at", t_point, ": R2 =",
        round(res$R2[1], 3), ", p =", res$`Pr(>F)`[1], "\n")
  }
}


# ── Remove outliers and re-run ───────────────────────────────────────────────

outliers <- excluded_samples

bio_no_outlier  <- bio_cols[!bio_cols %in% outliers]
meta_no_outlier <- meta[bio_no_outlier, ]
log2_no_outlier <- log2_clean[, bio_no_outlier]

dist_no_out <- dist(t(log2_no_outlier), method = "euclidean")

perm_no_out <- adonis2(dist_no_out ~ Treatment * Time,
                       data = meta_no_outlier,
                       permutations = 9999)
cat("\n══ PERMANOVA without outliers (marginal) ══\n")
print(perm_no_out)

perm_full <- adonis2(dist_no_out ~ Treatment * Time,
                     data = meta_no_outlier,
                     permutations = 9999,
                     by = "terms")
cat("\n══ PERMANOVA by terms (without outliers) ══\n")
print(perm_full)


# ── Pairwise PERMANOVA (without outliers) ────────────────────────────────────

comparisons <- list(
  "THz vs Control @ 10min" = list(treat = c("THz", "Control"), time = "10min"),
  "THz vs Control @ 45min" = list(treat = c("THz", "Control"), time = "45min"),
  "IR vs Control @ 10min"  = list(treat = c("IR", "Control"),  time = "10min"),
  "IR vs Control @ 45min"  = list(treat = c("IR", "Control"),  time = "45min"),
  "THz vs IR @ 10min"      = list(treat = c("THz", "IR"),      time = "10min"),
  "THz vs IR @ 45min"      = list(treat = c("THz", "IR"),      time = "45min"),
  "Control: 45 vs 10 min"  = list(treat = c("Control"), time = c("10min", "45min")),
  "THz: 45 vs 10 min"      = list(treat = c("THz"),     time = c("10min", "45min")),
  "IR: 45 vs 10 min"       = list(treat = c("IR"),      time = c("10min", "45min"))
)

cat("\n══ Pairwise PERMANOVA (without outliers) ══\n")
cat(sprintf("%-30s %6s %8s %8s\n", "Comparison", "R2", "F", "p-value"))
cat(paste(rep("-", 56), collapse = ""), "\n")

pairwise_results <- data.frame()

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]

  if (length(comp$treat) == 1) {
    idx <- meta_no_outlier$Treatment == comp$treat &
           meta_no_outlier$Time %in% comp$time
    m <- meta_no_outlier[idx, ]
    d <- dist(t(log2_no_outlier[, rownames(m)]), method = "euclidean")
    res <- adonis2(d ~ Time, data = m, permutations = 9999)
  } else {
    idx <- meta_no_outlier$Treatment %in% comp$treat &
           meta_no_outlier$Time == comp$time
    m <- meta_no_outlier[idx, ]
    d <- dist(t(log2_no_outlier[, rownames(m)]), method = "euclidean")
    res <- adonis2(d ~ Treatment, data = m, permutations = 9999)
  }

  r2    <- round(res$R2[1], 3)
  f_val <- round(res$F[1], 2)
  p_val <- res$`Pr(>F)`[1]

  cat(sprintf("%-30s %6.3f %8.2f %8.4f", comp_name, r2, f_val, p_val))
  if (p_val < 0.001) cat(" ***")
  else if (p_val < 0.01)  cat(" **")
  else if (p_val < 0.05)  cat(" *")
  else if (p_val < 0.1)   cat(" .")
  cat("\n")

  pairwise_results <- rbind(pairwise_results,
    data.frame(Comparison = comp_name, R2 = r2,
               F_stat = f_val, p_value = p_val))
}

write.csv(pairwise_results, file.path(tables_dir, "permanova_pairwise.csv"),
          row.names = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 11. LIMMA DIFFERENTIAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

meta_clean <- meta[!rownames(meta) %in% outliers, ]
log2_limma <- log2_clean[, rownames(meta_clean)]

cat("\nSamples for limma:", ncol(log2_limma), "\n")
print(table(meta_clean$Group))

group_factor <- factor(meta_clean$Group, levels = c(
  "Control_10min", "Control_45min",
  "THz_10min", "THz_45min",
  "IR_10min", "IR_45min"
))

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

fit <- lmFit(log2_limma, design)

contrasts <- makeContrasts(
  THz10_vs_Ctrl10 = THz_10min - Control_10min,
  THz45_vs_Ctrl45 = THz_45min - Control_45min,
  IR10_vs_Ctrl10  = IR_10min  - Control_10min,
  IR45_vs_Ctrl45  = IR_45min  - Control_45min,
  THz_vs_IR_10    = THz_10min - IR_10min,
  THz_vs_IR_45    = THz_45min - IR_45min,
  Time_in_Control = Control_45min - Control_10min,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

cat("\n══ limma results: significant metabolites (FDR < 0.05) ══\n")

results_list <- list()
for (coef_name in colnames(contrasts)) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "p")
  tt$Metabolite <- rownames(tt)
  tt$Contrast   <- coef_name
  results_list[[coef_name]] <- tt

  n_sig    <- sum(tt$adj.P.Val < 0.05)
  n_sig_fc <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.58)
  cat(coef_name, ": FDR<0.05 =", n_sig,
      "| FDR<0.05 & |FC|>1.5 =", n_sig_fc, "\n")
}

all_results <- bind_rows(results_list)

for (coef_name in colnames(contrasts)) {
  filename <- file.path(tables_dir, paste0("limma_", coef_name, ".csv"))
  write.csv(results_list[[coef_name]], file = filename, row.names = FALSE)
}

all_significant <- all_results[all_results$adj.P.Val < 0.05, ]
write.csv(all_significant, file = file.path(tables_dir, "limma_significant_FDR05.csv"),
          row.names = FALSE)

cat("Total significant hits:", nrow(all_significant), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# 12. VISUALIZATION: HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

sig_mets <- all_results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.58) %>%
  pull(Metabolite) %>% unique()

cat("\nUnique significant metabolites for heatmap:", length(sig_mets), "\n")

heatmap_data   <- log2_limma[sig_mets, ]
heatmap_scaled <- t(scale(t(heatmap_data)))

annotation_col <- data.frame(
  Treatment = meta_clean$Treatment,
  Time      = meta_clean$Time,
  row.names = rownames(meta_clean)
)

ann_colors <- list(
  Treatment = treatment_colors,
  Time      = time_colors
)

col_order <- rownames(meta_clean)[order(meta_clean$Treatment, meta_clean$Time)]

pdf(file.path(figures_dir, "06_heatmap.pdf"), width = 14, height = 18)
pheatmap(
  heatmap_scaled[, col_order],
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-3, 3, length.out = 101),
  fontsize_row      = 9,
  fontsize_col      = 10,
  fontsize          = 12,
  main = "Significant metabolites (FDR < 0.05, |FC| > 1.5)\nZ-score across samples",
  border_color      = NA
)
dev.off()


# ══════════════════════════════════════════════════════════════════════════════
# 13. VISUALIZATION: DOT PLOT — Top metabolites across contrasts
# ══════════════════════════════════════════════════════════════════════════════

# ── FIGURE 5B ────────────────────────────────────────────────────────────────

top_mets <- all_results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.58) %>%
  filter(Contrast %in% c("THz45_vs_Ctrl45", "THz_vs_IR_10")) %>%
  group_by(Contrast) %>%
  arrange(P.Value) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  pull(Metabolite) %>% unique()

contrast_levels <- c("THz10_vs_Ctrl10", "THz45_vs_Ctrl45",
                     "IR10_vs_Ctrl10", "IR45_vs_Ctrl45",
                     "THz_vs_IR_10", "THz_vs_IR_45")

contrast_labels <- c(
  "THz10_vs_Ctrl10" = "THz 10\nvs Ctrl 10",
  "THz45_vs_Ctrl45" = "THz 45\nvs Ctrl 45",
  "IR10_vs_Ctrl10"  = "IR 10\nvs Ctrl 10",
  "IR45_vs_Ctrl45"  = "IR 45\nvs Ctrl 45",
  "THz_vs_IR_10"    = "THz\nvs IR 10",
  "THz_vs_IR_45"    = "THz\nvs IR 45"
)

dot_data <- all_results %>%
  filter(Metabolite %in% top_mets,
         Contrast %in% contrast_levels) %>%
  mutate(
    Sig = ifelse(adj.P.Val < 0.05, "FDR < 0.05", "NS"),
    Contrast = factor(Contrast, levels = contrast_levels)
  )

p_dot <- ggplot(dot_data, aes(x = Contrast, y = Metabolite)) +
  geom_point(aes(size = -log10(adj.P.Val),
                 color = logFC,
                 shape = Sig),
             stroke = 0.8) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, name = "log2FC") +
  scale_size_continuous(range = c(2, 7), name = expression(-log[10](FDR))) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "NS" = 1),
                     name = "Significance") +
  scale_x_discrete(labels = contrast_labels) +
  labs(y = NULL, x = NULL) +
  theme_pub(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 11, lineheight = 0.9),
    axis.text.y      = element_text(size = 11, face = "italic"),
    legend.key.size  = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(figures_dir, "07_dotplot.pdf"), p_dot,
       width = 11, height = 13, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 14. VISUALIZATION: LOLLIPOP — THz 45 min vs Control 45 min
# ══════════════════════════════════════════════════════════════════════════════

thz45 <- results_list[["THz45_vs_Ctrl45"]] %>%
  filter(P.Value < 0.05, abs(logFC) > 0.58) %>%
  arrange(logFC) %>%
  mutate(Metabolite = fct_inorder(Metabolite),
         Direction = ifelse(logFC > 0, "Up", "Down"))

p_lolli45 <- ggplot(thz45, aes(x = logFC, y = Metabolite, color = Direction)) +
  geom_segment(aes(x = 0, xend = logFC, yend = Metabolite),
               linewidth = 1.0) +
  geom_point(aes(size = -log10(adj.P.Val))) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC"),
                     name = "Direction") +
  scale_size_continuous(range = c(2.5, 6),
                        name = expression(-log[10](FDR))) +
  geom_vline(xintercept = 0, lty = 1, color = "grey40", linewidth = 0.5) +
  labs(title = "THz 45 min vs Control 45 min",
       subtitle = "Nominal p < 0.05 & |FC| > 1.5",
       x = expression(log[2]~Fold~Change),
       y = NULL) +
  theme_pub(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10)
  )

ggsave(file.path(figures_dir, "08_lollipop_THz45.pdf"), p_lolli45,
       width = 11, height = 14, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 15. METABOANALYST INPUT PREPARATION
# ══════════════════════════════════════════════════════════════════════════════

lipid_patterns <- c("^PI\\(", "^PE\\(", "^PC\\(", "^PS\\(", "^SM\\(",
                     "^PA\\(", "^PG\\(", "^CL\\(", "^BMP\\(", "^GC\\(",
                     "^Ceramide", "^THC ", "^Plasmalogen",
                     "^Cholesteryl", "^DHC\\(")

# Reference metabolome (all measured, non-lipid)
all_measured <- rownames(data_clean) %>%
  str_replace("_(pos|neg).*$", "") %>%
  trimws() %>% unique()

all_measured_clean <- all_measured[!str_detect(all_measured,
                     paste(lipid_patterns, collapse = "|"))]

writeLines(all_measured_clean, file.path(tables_dir, "reference_metabolome.txt"))

# THz45 vs Ctrl45 — metabolites for pathway analysis
thz45_pathway <- results_list[["THz45_vs_Ctrl45"]] %>%
  filter(P.Value < 0.05) %>%
  mutate(
    Clean_Name = str_replace(Metabolite,
                             "_(pos|neg|principal|qualifier).*$", ""),
    Clean_Name = str_replace(Clean_Name, "_$", ""),
    Clean_Name = trimws(Clean_Name)
  ) %>%
  select(Clean_Name, logFC, P.Value)

write.csv(thz45_pathway, file.path(tables_dir, "pathway_THz45_vs_Ctrl45.csv"),
          row.names = FALSE)

cat("\nMetaboAnalyst input prepared:\n")
cat("  THz45 vs Ctrl45:", nrow(thz45_pathway), "metabolites\n")
cat("  Reference metabolome:", length(all_measured_clean), "metabolites\n")


# ══════════════════════════════════════════════════════════════════════════════
cat("\n══ Pipeline complete ══\n")
cat("Figures saved to ", figures_dir, "\n", sep = "")
cat("Tables saved to ", tables_dir, "\n", sep = "")
# ══════════════════════════════════════════════════════════════════════════════
