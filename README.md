# scfa_cor
R code for correlation analysis between bovine gut microbial genera and short-chain fatty acids

# SCFA Correlation Analysis with Microbial Genera

This repository contains the R script and relevant data files used for correlation analysis between bovine gut microbial genera and short-chain fatty acids (SCFAs).

---

## 📊 Project Overview

The analysis includes:

- Correlation of microbial genera with SCFAs (e.g., acetate, butyrate)
- Multiple testing correction (FDR)
- Heatmap visualization of significant correlations
- PCoA (Principal Coordinates Analysis) and PERMANOVA for group clustering
- Environmental fitting of SCFAs onto PCoA ordination

---

## 📁 Files in This Repository

| File | Description |
|------|-------------|
| `scfa_correlation_analysis.R` | Main R script used for analysis |
| `Morgan_Results.xlsx` | Input data: Microbial genus abundance |
| `Morgan SCFA results_Excel.xlsx` | Input data: SCFA concentrations |
| `README.md` | This readme file |

---

## 📦 Required R Packages

Before running the script, make sure to install the following R packages:

```r
install.packages("readxl")
install.packages("vegan")
install.packages("pheatmap")
