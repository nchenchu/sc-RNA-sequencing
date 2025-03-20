library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Define Output Directory
output_dir <- "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/Clones"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
print(paste("Files will be saved in:", output_dir))

# Define the correct path to your RDS file
seurat_obj <- readRDS("/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV R848 scRNAseq cellcycle_corrected.rds")

# Verify if the object is loaded
print(seurat_obj) 

# Check available clusters first
table(seurat_obj$seurat_clusters)

#Subset B cell clusters
b_cell_clusters <- c(1, 4, 28, 29)
b_cells <- subset(seurat_obj, idents = b_cell_clusters)


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
b_cell_rds <- file.path(output_dir, "B_Cell_Subclusters.rds")
saveRDS(b_cells, b_cell_rds)

print(paste("B Cell Subclusters saved at:", b_cell_rds))



DimPlot(t_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
# Define plot save path
b_cell_umap_path <- file.path(output_dir, "B_Cell_Subclusters_UMAP.png")

# Save the plot
ggsave(b_cell_umap_path, plot = DimPlot(b_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters"), width = 8, height = 6)

# Save UMAP Plot
b_cell_umap_path <- file.path(output_dir, "B_Cell_UMAP.png")
umap_plot <- DimPlot(b_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ggtitle("B Cell UMAP Clusters")
ggsave(filename = b_cell_umap_path, plot = umap_plot, width = 8, height = 6, dpi = 300)
message("B Cell UMAP plot saved at: ", b_cell_umap_path)
