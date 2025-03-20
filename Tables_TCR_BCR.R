# Load Required Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(ineq)
library(reshape2)

# Define Local Directory for Saving Files
output_dir <- "/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Results/"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load T Cell Data
T_scd <- readRDS("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV_T_Cells_Annotated.rds")

# Load B Cell Data
B_scd <- readRDS("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV_B_Cells_Annotated.rds")



# Set Default Assay to SCT
DefaultAssay(T_scd) <- 'SCT'

# Set Default Assay to SCT
DefaultAssay(B_scd) <- 'SCT'

# Extract Metadata from Seurat Object
t_df <- T_scd@meta.data  

# Extract Metadata from Seurat Object
b_df <- B_scd@meta.data  


# Check column names to ensure they contain "Cluster_ID", "Sample_ID", etc.
print(colnames(t_df))

# Check column names to ensure they contain "Cluster_ID", "Sample_ID", etc.
print(colnames(b_df))


# View metadata columns
colnames(T_scd@meta.data)

# View metadata columns
colnames(B_scd@meta.data)

# Display first few rows to inspect
head(T_scd@meta.data, 10)
# Display first few rows to inspect
head(B_scd@meta.data, 10)


# Display unique cell types in the dataset
unique(T_scd@meta.data$cell_type)

table(T_scd@meta.data$cell_type)

# Display unique cell types in the dataset
unique(B_scd@meta.data$cell_type)

table(B_scd@meta.data$cell_type)

library(dplyr)

# Separate CD4+, CD8+, and Unknown T cells
CD4_T_cells <- T_scd@meta.data %>% filter(cell_type == "CD4+ T Cells")
CD8_T_cells <- T_scd@meta.data %>% filter(cell_type == "CD8+ T Cells")
Unknown_T_cells <- T_scd@meta.data %>% filter(cell_type == "Unknown")

# Separate CD4+, CD8+, and Unknown B cells
CD4_B_cells <- B_scd@meta.data %>% filter(cell_type == "CD4+ B Cells")
CD8_B_cells <- B_scd@meta.data %>% filter(cell_type == "CD8+ B Cells")
Unknown_B_cells <- B_scd@meta.data %>% filter(cell_type == "Unknown")

# Save T cell datasets to CSV
write.csv(CD4_T_cells, file.path(output_dir, "CD4_T_cells.csv"), row.names = FALSE)
write.csv(CD8_T_cells, file.path(output_dir, "CD8_T_cells.csv"), row.names = FALSE)
write.csv(Unknown_T_cells, file.path(output_dir, "Unknown_T_cells.csv"), row.names = FALSE)

# Save B cell datasets to CSV
write.csv(CD4_B_cells, file.path(output_dir, "CD4_B_cells.csv"), row.names = FALSE)
write.csv(CD8_B_cells, file.path(output_dir, "CD8_B_cells.csv"), row.names = FALSE)
write.csv(Unknown_B_cells, file.path(output_dir, "Unknown_B_cells.csv"), row.names = FALSE)

cat("CSV files have been successfully saved in the working directory!")

