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


# Check metadata columns
colnames(seurat_obj@meta.data)

# Extract relevant metadata (Modify if TCR/BCR columns are different)
tcr_bcr_table <- seurat_obj@meta.data %>%
  select(orig.ident, seurat_clusters, Condition) %>%
  mutate(TCR_BCR_Clone = rownames(.))  # Modify if there is a specific clonotype column

# Save to CSV
csv_path <- file.path(output_dir, "TCR_BCR_Clone_Table.csv")
write.csv(tcr_bcr_table, csv_path, row.names = FALSE)

print(paste("TCR/BCR Clone Table saved at:", csv_path))


##########
# Set SCT as active assay
DefaultAssay(seurat_obj) <- "SCT"

# Subset T cell clusters
t_cell_clusters <- c(2, 6, 15, 17, 21)
t_cells <- subset(seurat_obj, idents = t_cell_clusters)
# 
# # Run clustering
# t_cells <- FindNeighbors(t_cells, dims = 1:20)
# t_cells <- FindClusters(t_cells, resolution = 0.5)
# t_cells <- RunUMAP(t_cells, dims = 1:20)
# 
# # Save T cell subcluster object
# t_cell_rds <- file.path(output_dir, "T_Cell_Subclusters.rds")
# saveRDS(t_cells, t_cell_rds)
# print(paste("T Cell Subclusters saved at:", t_cell_rds))

#########

# Ensure SCT is the active assay
DefaultAssay(t_cells) <- "SCT"

# Scale the data before PCA
t_cells <- ScaleData(t_cells, verbose = TRUE)

# Perform PCA
t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))

# Compute neighbors and clustering
t_cells <- FindNeighbors(t_cells, dims = 1:20, graph.name = "SCT_nn")
t_cells <- FindClusters(t_cells, graph.name = "SCT_nn", resolution = 0.5)

# Run UMAP for visualization
t_cells <- RunUMAP(t_cells, dims = 1:20)

# Save the subclustered Seurat object
t_cell_rds <- file.path(output_dir, "T_Cell_Subclusters.rds")
saveRDS(t_cells, t_cell_rds)

print(paste("T Cell Subclusters saved at:", t_cell_rds))


DimPlot(t_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
# Define plot save path
t_cell_umap_path <- file.path(output_dir, "T_Cell_Subclusters_UMAP.png")

# Save the plot
ggsave(t_cell_umap_path, plot = DimPlot(t_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters"), width = 8, height = 6)

# Print confirmation
print(paste("T Cell UMAP saved at:", t_cell_umap_path))


