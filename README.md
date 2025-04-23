# Analysis of DNA Methylation and Regenerative Capacity in Invertebrates

This repository contains the R Markdown notebook used for the analyses in the study:

**"Regeneration and Genomic Methylation in Invertebrates: Insights into Epigenome Evolution and Cell Lineage Plasticity."**

## Repository Contents

- **`Analysis_Notebook.Rmd`**: Main analysis notebook, including data processing, phylogenetic comparative methods, and visualization of results.
- **`Grafted_Tree_Preparation.Rmd`**: Notebook for preparing grafted tree with full species set.
- **`Methylation_level.csv`**: Data on 5-methylcytosine CpG methylation values and regeneration categories for 175 species

## Description

The Analysis notebook analyzes genome-wide CpG methylation levels in relation to regenerative capacity across invertebrate species. It includes:
- Data import and preprocessing
- Phylogenetic signal testing (Pagel's λ, Blomberg’s K)
- Phylogenetic Generalized Least Squares (PGLS) analysis
- Visualization of results (phylogenetic trees, scatter plots)
  
The Grafting notebook prepares the full tree from the time-calibrated tree and taxonomy information to determine the rest of the species placements
- 

## Requirements

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
