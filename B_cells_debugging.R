# Load required libraries
library(Seurat)
library(dplyr)

# Define Local Directory for Saving Files
output_dir <- "/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Results/b_cells_debugging"

# Load the RDS file and extract metadata
scd_metadata <- readRDS("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV R848 scRNAseq cellcycle_corrected.rds")

# Define B cell clusters (change "Cluster_ID" to "seurat_clusters")
b_cell_clusters <- c(0, 1, 4, 28, 29)  # Adjust clusters as needed

