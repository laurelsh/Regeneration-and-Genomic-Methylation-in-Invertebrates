# Phylogenetic Comparative Analysis — Regeneration and CpG Methylation

This repository contains the data and R scripts used to analyze relationships between genome-wide CpG methylation and regenerative capacity across metazoans.


## Overview

### **1. Grafted_Tree_Preparation.R**
Constructs a complete 175-tip phylogeny by grafting additional taxa from the full topology (`Full_topology.nwk`) onto the calibrated tree (`Tree_calibrated.nwk`).

Main steps:
1. Retrieve taxonomy for existing and missing taxa via GBIF.
2. Label internal nodes by rank (Family, Order, Class, Phylum).
3. Estimate clade-specific speciation rates using Yule models.
4. Graft new taxa to the most appropriate nodes.
5. Save the resulting tree as `data/Tree_grafted.nwk`.

**Required R packages:**
`ape`, `phytools`, `diversitree`, `treeio`, `ggtree`,  
`rgbif`, `dplyr`, `pbapply`, `stringr`.

---

### **2. PGLS_analysis.R**
Performs Phylogenetic Generalized Least Squares (PGLS) tests to examine the relationship between CpG methylation and regenerative capacity.

Main steps:
1. Load either the calibrated or grafted tree (`Tree_calibrated.nwk` or `Tree_grafted.nwk`).
2. Align tree tips with data (`Methylation_data.csv`).
3. Apply a Tukey transformation to methylation values.
4. Fit PGLS models using Pagel’s λ correlation structure.
5. Compute estimated marginal means and pairwise contrasts.
6. Save results to a `results/` folder, including plots and CSV tables.

**Required R packages:**
`ape`, `phytools`, `nlme`, `emmeans`, `dplyr`,  
`ggplot2`, `rcompanion`, `nortest`, `moments`.

---

## How to Run the Scripts

1. **Download or clone this repository** to your computer.  
   You can click the green **Code → Download ZIP** button on GitHub, then extract it anywhere on your computer.

2. **Open R or RStudio** and set your working directory to the project folder:
   ```r
   setwd("path/to/Phylogenetic_comparative_analysis")

   Outputs
	•	data/Tree_grafted.nwk — new 175-tip phylogeny
	•	results/ — automatically generated folder containing:
	•	PGLS_calibrated.pdf
	•	PGLS_grafted.pdf
	•	Contrasts_*.csv
	•	Summary_*.csv

## Citation

If using or adapting this repository, please cite:

Hiebert, L. and Yi, S. PNAS (2025, in revision).

## License

All scripts are released under the MIT License.
Data files are provided under CC-BY 4.0, allowing reuse with attribution.
