# Analysis of DNA Methylation and Regenerative Capacity in Invertebrates

This repository contains the R Markdown notebook used for the analysis in the study:

**"Regeneration and Genomic Methylation in Invertebrates: Insights into Epigenome Evolution and Cell Lineage Plasticity."**

## Repository Contents

- **`Analysis_Notebook.Rmd`**: Main analysis notebook, including data processing, phylogenetic comparative methods, and visualization of results.

## Description

The notebook analyzes genome-wide CpG methylation levels in relation to regenerative capacity across 175 invertebrate species. It includes:
- Data import and preprocessing
- Phylogenetic signal testing (Pagel's λ, Blomberg’s K)
- Phylogenetic Generalized Least Squares (PGLS) analysis
- Visualization of results (phylogenetic trees, scatter plots)

## Requirements

- R version X.X.X
- Suggested packages:
  - `ape`
  - `phytools`
  - `nlme`
  - `caper`
  - `ggplot2`
  - `dplyr`
  - `readr`

### Install required packages:

```r
install.packages(c("ape", "phytools", "nlme", "caper", "ggplot2", "dplyr", "readr"))
```

## How to Run

	1.	Open Analysis_Notebook.Rmd in RStudio.
	2.	Knit the document or run chunks step-by-step to reproduce the analysis.

## Citation

If you use this code, please cite:

Hiebert and Yi. Regeneration and Genomic Methylation in Invertebrates: Insights into Epigenome Evolution and Cell Lineage Plasticity. (2025) [Journal Name, DOI]
