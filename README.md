# Fluorescence and Colony Size analysis 

This repository contains scripts to analyse colony fluorescence and colony size measurements using two complementary approaches: min–max scaled fluorescence-to-size ratios and Z-scored fluorescence-to-size ratios.

The analysis workflow consists of:
1. Extraction of per-colony fluorescence intensities from Typhoon scanner images (Python)
2. Colony size quantification using Iris 
3. Statistical analysis, quality control, and visualisation in R

## Requirements

### Python (image processing)
- Python ≥ 3.8  
- Packages: numpy, pandas, pillow  

### R (analysis and plotting)
- R ≥ 4.2  
- Packages: tidyverse, ggplot2, ggbeeswarm, rstatix, ggpubr  

## Image processing (Typhoon scanner)

Per-colony fluorescence pixel intensities are extracted from gridded Typhoon scanner images using a custom Python script:

scripts/image_processing/extract_pixel_intensities.py

The script converts images to grayscale, overlays a grid corresponding to colony positions, and computes mean pixel intensities for each colony, outputting a CSV file for downstream analysis.

## Run the analysis (R)

The main R script performs quality control, normalisation, statistical testing, and visualisation using both scaled and Z-score approaches.

```bash
Rscript scripts/analysis_30C_scaled_zscore.R \
  --flu path/to/fl_Final_dataset.csv \
  --size path/to/size_Final_dataset.csv \
  --out results \
  --date 23-12-2025
```

Notes
Statistical comparisons are based on pairwise t-tests with Benjamini–Hochberg correction.


