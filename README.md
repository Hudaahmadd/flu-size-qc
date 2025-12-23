# Fluorescence / Colony Size analysis (30C)

This repository contains an R script to analyse colony fluorescence
and colony size measurements at 30 °C using two approaches:
min–max scaled ratios and Z-scored ratios.

## Requirements
- R (≥ 4.2)
- Packages: tidyverse, ggplot2, ggbeeswarm, rstatix, ggpubr

## Run the analysis

```bash
Rscript scripts/analysis_30C_scaled_zscore.R \
  --flu path/to/fl_Final_dataset.csv \
  --size path/to/size_Final_dataset.csv \
  --out results \
  --date 23-12-2025
