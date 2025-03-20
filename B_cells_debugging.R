# Load required libraries
library(Seurat)
library(dplyr)

# Define Local Directory for Saving Files
output_dir <- "file_path"

# Load the RDS file and extract metadata
scd_metadata <- readRDS("file_path")

# Define B cell clusters (change "Cluster_ID" to "seurat_clusters")
b_cell_clusters <- c(0, 1, 4, 28, 29)  # Adjust clusters as needed

