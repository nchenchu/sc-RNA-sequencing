library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Define Output Directory
output_dir <- "file_path"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
print(paste("Files will be saved in:", output_dir))

# Define the correct path to your RDS file
seurat_obj <- readRDS("file_path")
# Verify if the object is loaded
print(seurat_obj) 

# Check available clusters first
table(seurat_obj$seurat_clusters)

#Subset B cell clusters
b_cell_clusters_zero <- c(0)
b_cells <- subset(seurat_obj, idents = b_cell_clusters_zero)


# Ensure SCT is the active assay
DefaultAssay(b_cells) <- "SCT"

# Scale the data before PCA
b_cells <- ScaleData(b_cells, verbose = TRUE)

# Perform PCA
b_cells <- RunPCA(b_cells, features = VariableFeatures(object = b_cells))

# Compute neighbors and clustering
b_cells <- FindNeighbors(b_cells, dims = 1:20, graph.name = "SCT_nn")
b_cells <- FindClusters(b_cells, graph.name = "SCT_nn", resolution = 0.5)

# Run UMAP for visualization
b_cells <- RunUMAP(b_cells, dims = 1:20)

# Save the subclustered Seurat object
b_cell_rds <- file.path(output_dir, "B_Cell_Subclusters_zero.rds")
saveRDS(b_cells, b_cell_rds_zero)

print(paste("B Cell Subclusters saved at:", b_cell_rds_zero))
