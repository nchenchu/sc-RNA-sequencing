# Single-cell RNA Sequencing

This repository contains my R scripts for single-cell RNA sequencing (scRNA-seq) analysis.  
These scripts are useful for subsetting B cells and T cells, analyzing clonal types through TCR/BCR data, and generating summary tables.  

These scripts were written as part of my learning process, so the code is simple and beginner-friendly.  
Some parts of the code were created with help from AI tools to support my learning and understanding.

---

#  Folder Structure

- src/ : All analysis R scripts  
- data/ : No data is provided because the datasets I used cannot be shared  
- results/ : Output files are not included for privacy reasons  
- docs/ : Currently not used (NA)

---

#  Requirements

You need to install a few R packages that are used in the scripts.  
Packages include:

- Seurat  
- dplyr  
- tidyr  
- ggplot2  
- data.table  
- readr  
- optparse
- patchwork
- ineq
- ggpubr
- pheatmap
- reshape2

You can install them manually in R or use a `requirements.R` installer script.

---

#  What Each Script Does

- B_Cells_subsetting.R – Subsets B cells from a Seurat object  
- T_Cells_subsetting.R – Subsets T cells  
- B_cell_clones.R / B_cell_clones_zero.R – B-cell clonotype analysis  
- T_cell_clones.R – T-cell clonotypes  
- TCR_BCR_clones.R – Combined clonotype analysis  
- Tables_TCR_BCR.R – Generates summary tables  
- Gini_coefficientcalculation.R – Computes Gini coefficient for clonal diversity  
- MRV_surface_reading.R – Reads surface or plotting-related MRV files  

---

#  Notes

- These scripts were written during my learning journey in scRNA-seq.  
- The code style is simple and beginner-friendly.  
- Some parts were created with the assistance of AI tools to help me understand R programming and scRNA-seq concepts.  








