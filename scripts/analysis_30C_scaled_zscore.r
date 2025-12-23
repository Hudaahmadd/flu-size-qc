###############################################
# FINAL (30C ONLY): SCALED + Z-SCORE (MATCHING STYLE)
# - One script, one output folder: scalling/23-12-2025
# - BOTH pipelines keep A/B in the data (as rows) but do NOT show them in plot naming
# - BOTH do: log2 ratio -> Tukey 1.5×IQR QC (per Gene) -> pairwise t-tests (BH) -> brackets+stars
# - Exports (for EACH pipeline):
#     1) dat30_for_swarmplot.csv
#     2) outliers_30C_Tukey1.5IQR.csv
#     3) pvalues_pairwise_30C_BH_USED_ON_PLOT.csv
#     4) TIFF plot (no legend)
###############################################

library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)

# ----------------------------
# Paths + output folder
# ----------------------------
flu_path  <- "/Users/hxa105@student.bham.ac.uk/Library/CloudStorage/OneDrive-UniversityofBirmingham/Documents/PhD/Sam/08-12-25/analysis/scalling/fl_Final_dataset.csv"
size_path <- "/Users/hxa105@student.bham.ac.uk/Library/CloudStorage/OneDrive-UniversityofBirmingham/Documents/PhD/Sam/08-12-25/analysis/scalling/size_Final_dataset.csv"

base_dir <- "/Users/hxa105@student.bham.ac.uk/Library/CloudStorage/OneDrive-UniversityofBirmingham/Documents/PhD/Sam/08-12-25/analysis/scalling"
out_dir  <- file.path(base_dir, "23-12-2025")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper functions
# ----------------------------
scale_mm <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (is.na(rng) || rng == 0) return(rep(NA_real_, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

stars_from_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Adds Tukey fences per Gene (30C only here)
add_tukey_fences_gene <- function(df, value_col) {
  df %>%
    group_by(Gene) %>%
    mutate(
      Q1 = quantile({{ value_col }}, 0.25, na.rm = TRUE),
      Q3 = quantile({{ value_col }}, 0.75, na.rm = TRUE),
      IQRv = Q3 - Q1,
      lo = Q1 - 1.5 * IQRv,
      hi = Q3 + 1.5 * IQRv
    ) %>%
    ungroup()
}

# Stacks bracket y-positions for stat_pvalue_manual
stack_brackets <- function(pw_df, y_vals, step_frac = 0.10) {
  ymax <- max(y_vals, na.rm = TRUE)
  yrng <- diff(range(y_vals, na.rm = TRUE))
  pw_df %>%
    arrange(p.adj) %>%
    mutate(y.position = ymax + row_number() * step_frac * yrng)
}

save_tiff_plot <- function(plot_obj, filename, width_mm = 180, height_mm = 140, res = 300) {
  tiff(filename = filename,
       width = width_mm, height = height_mm, units = "mm",
       res = res, compression = "lzw")
  print(plot_obj)
  dev.off()
}

# ----------------------------
# 1) Read + standardise column names + clean
# ----------------------------
flu <- read.csv(flu_path, check.names = FALSE)
size <- read.csv(size_path, check.names = FALSE)

colnames(flu)  <- c("Gene", "Flu_30C_A", "Flu_30C_B", "Flu_37C_A", "Flu_37C_B")
colnames(size) <- c("Gene", "Size_30C_A","Size_30C_B","Size_37C_A","Size_37C_B")

flu_clean <- flu %>%
  slice(-(1:2)) %>%
  mutate(across(-Gene, as.numeric))

size_clean <- size %>%
  slice(-(1:2)) %>%
  mutate(across(-Gene, as.numeric))

############################################################
# A) SCALED PIPELINE (30C only)
############################################################

# A1) Min–max scale each column independently
flu_scaled <- flu_clean %>% mutate(across(starts_with("Flu_"), scale_mm))
size_scaled <- size_clean %>% mutate(across(starts_with("Size_"), scale_mm))

# A2) Compute log2(scaled flu / scaled size) for 30C A/B only
ratio_scaled_wide_30C <- flu_scaled %>%
  bind_cols(size_scaled %>% select(-Gene)) %>%
  transmute(
    Gene,
    Ratio_30C_A = log2((Flu_30C_A + 1e-6) / (Size_30C_A + 1e-6)),
    Ratio_30C_B = log2((Flu_30C_B + 1e-6) / (Size_30C_B + 1e-6))
  )

ratio_scaled_long_30C <- ratio_scaled_wide_30C %>%
  pivot_longer(
    cols = starts_with("Ratio_"),
    names_to = "Replicate",
    values_to = "log2_scaled_Flu_per_Size"
  ) %>%
  mutate(
    Temperature = "30C",
    Replicate = gsub("Ratio_30C_", "", Replicate)
  ) %>%
  filter(is.finite(log2_scaled_Flu_per_Size))

# A3) QC on log2 ratio (Tukey 1.5×IQR per Gene)
scaled_fenced_30C <- ratio_scaled_long_30C %>%
  add_tukey_fences_gene(value_col = log2_scaled_Flu_per_Size) %>%
  mutate(
    outlier_flag = case_when(
      !is.finite(log2_scaled_Flu_per_Size) ~ "non-finite",
      log2_scaled_Flu_per_Size < lo ~ "low",
      log2_scaled_Flu_per_Size > hi ~ "high",
      TRUE ~ "ok"
    )
  )

outliers_scaled_30C <- scaled_fenced_30C %>%
  filter(outlier_flag != "ok") %>%
  select(Gene, Temperature, Replicate, log2_scaled_Flu_per_Size, lo, hi, Q1, Q3, IQRv, outlier_flag) %>%
  arrange(Gene, Replicate)

dat30_scaled <- scaled_fenced_30C %>%
  filter(outlier_flag == "ok") %>%
  transmute(
    Gene = as.factor(Gene),
    Temperature,
    Replicate,
    log2_Flu_per_Size = log2_scaled_Flu_per_Size
  )

# A4) Pairwise tests + bracket table (BH)
pw_scaled_30C <- dat30_scaled %>%
  pairwise_t_test(log2_Flu_per_Size ~ Gene, p.adjust.method = "BH") %>%
  mutate(p.signif = stars_from_p(p.adj)) %>%
  stack_brackets(y_vals = dat30_scaled$log2_Flu_per_Size, step_frac = 0.10)

pvals_scaled_used_on_plot <- pw_scaled_30C %>%
  select(group1, group2, n1, n2, p, p.adj, p.adj.signif, p.signif, y.position)

# A5) Export CSVs (scaled)
write.csv(dat30_scaled, file.path(out_dir, "A_scaled_dat30_for_swarmplot.csv"), row.names = FALSE)
write.csv(outliers_scaled_30C, file.path(out_dir, "A_scaled_outliers_30C_Tukey1.5IQR.csv"), row.names = FALSE)
write.csv(pvals_scaled_used_on_plot, file.path(out_dir, "A_scaled_pvalues_pairwise_30C_BH_USED_ON_PLOT.csv"), row.names = FALSE)

# A6) Plot + TIFF (scaled; no legend)
yrng_scaled <- diff(range(dat30_scaled$log2_Flu_per_Size, na.rm = TRUE))

p30_scaled <- ggplot(dat30_scaled, aes(x = Gene, y = log2_Flu_per_Size, color = Gene)) +
  geom_beeswarm(size = 1.6, alpha = 0.8, cex = 1.1) +
  theme_minimal(base_size = 11) +
  labs(
    title = expression("30"~degree*C~": log"[2]~"(Scaled Fluorescence / Scaled Colony Size)"),
    x = "Strain",
    y = expression(log[2]~"(Scaled Fluorescence / Scaled Size)")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(
    min(dat30_scaled$log2_Flu_per_Size, na.rm = TRUE),
    max(pvals_scaled_used_on_plot$y.position) + 0.10 * yrng_scaled
  )) +
  stat_pvalue_manual(
    pvals_scaled_used_on_plot,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  )

save_tiff_plot(p30_scaled, file.path(out_dir, "A_scaled_Figure_30C_log2_Flu_per_Size.tiff"))

############################################################
# B) Z-SCORE PIPELINE (30C only) — MATCHES SCALED STYLE
# IMPORTANT: Z-score is computed on log2(raw Flu/Size) values (30C),
# then QC is applied on Z-scores per Gene, then plot uses Z-score.
############################################################

# B1) Compute raw log2(flu/size) for 30C A/B only
ratio_raw_wide_30C <- flu_clean %>%
  bind_cols(size_clean %>% select(-Gene)) %>%
  transmute(
    Gene,
    Log2_30C_A = log2((Flu_30C_A + 1e-6) / (Size_30C_A + 1e-6)),
    Log2_30C_B = log2((Flu_30C_B + 1e-6) / (Size_30C_B + 1e-6))
  )

ratio_raw_long_30C <- ratio_raw_wide_30C %>%
  pivot_longer(
    cols = starts_with("Log2_"),
    names_to = "Replicate",
    values_to = "log2_Flu_per_Size"
  ) %>%
  mutate(
    Temperature = "30C",
    Replicate = gsub("Log2_30C_", "", Replicate)
  ) %>%
  filter(is.finite(log2_Flu_per_Size))

# B2) Z-score within 30C across ALL colonies (matches your earlier definition)
ratio_z_30C <- ratio_raw_long_30C %>%
  mutate(
    Zscore = (log2_Flu_per_Size - mean(log2_Flu_per_Size, na.rm = TRUE)) /
      sd(log2_Flu_per_Size, na.rm = TRUE)
  ) %>%
  filter(is.finite(Zscore))

# B3) QC on Z-score (Tukey 1.5×IQR per Gene) — same QC logic as scaled
z_fenced_30C <- ratio_z_30C %>%
  add_tukey_fences_gene(value_col = Zscore) %>%
  mutate(
    outlier_flag = case_when(
      !is.finite(Zscore) ~ "non-finite",
      Zscore < lo ~ "low",
      Zscore > hi ~ "high",
      TRUE ~ "ok"
    )
  )

outliers_z_30C <- z_fenced_30C %>%
  filter(outlier_flag != "ok") %>%
  select(Gene, Temperature, Replicate, log2_Flu_per_Size, Zscore, lo, hi, Q1, Q3, IQRv, outlier_flag) %>%
  arrange(Gene, Replicate)

dat30_z <- z_fenced_30C %>%
  filter(outlier_flag == "ok") %>%
  transmute(
    Gene = as.factor(Gene),
    Temperature,
    Replicate,
    log2_Flu_per_Size,  # kept for reference/traceability
    Zscore
  )

# B4) Pairwise tests + bracket table (BH) on Zscore
pw_z_30C <- dat30_z %>%
  pairwise_t_test(Zscore ~ Gene, p.adjust.method = "BH") %>%
  mutate(p.signif = stars_from_p(p.adj)) %>%
  stack_brackets(y_vals = dat30_z$Zscore, step_frac = 0.10)

pvals_z_used_on_plot <- pw_z_30C %>%
  select(group1, group2, n1, n2, p, p.adj, p.adj.signif, p.signif, y.position)

# B5) Export CSVs (z-score)
write.csv(dat30_z, file.path(out_dir, "B_zscore_dat30_for_swarmplot.csv"), row.names = FALSE)
write.csv(outliers_z_30C, file.path(out_dir, "B_zscore_outliers_30C_Tukey1.5IQR.csv"), row.names = FALSE)
write.csv(pvals_z_used_on_plot, file.path(out_dir, "B_zscore_pvalues_pairwise_30C_BH_USED_ON_PLOT.csv"), row.names = FALSE)

# B6) Plot + TIFF (z-score; no legend; naming matches scaled style)
yrng_z <- diff(range(dat30_z$Zscore, na.rm = TRUE))

p30_z <- ggplot(dat30_z, aes(x = Gene, y = Zscore, color = Gene)) +
  geom_beeswarm(size = 1.6, alpha = 0.8, cex = 1.1) +
  theme_minimal(base_size = 11) +
  labs(
    title = expression("30"~degree*C~": Z-scored log"[2]~"(Fluorescence / Colony Size)"),
    x = "Strain",
    y = "Z-score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(
    min(dat30_z$Zscore, na.rm = TRUE),
    max(pvals_z_used_on_plot$y.position) + 0.10 * yrng_z
  )) +
  stat_pvalue_manual(
    pvals_z_used_on_plot,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  )

save_tiff_plot(p30_z, file.path(out_dir, "B_zscore_Figure_30C_Zscore_log2_Flu_per_Size.tiff"))

cat("Done. Outputs saved in:\n", out_dir, "\n", sep = "")

