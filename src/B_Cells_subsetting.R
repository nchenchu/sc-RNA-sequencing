data <- readRDS("file path")

library(Seurat)

# Open the RDS file and load only metadata
con <- gzfile("file path")
scd_metadata <- readRDS(con)@meta.data

# Check metadata structure
str(scd_metadata)

# Extract barcodes of B cell clusters
b_cell_clusters <- c("0", "1", "4", "28","29")  # Cluster IDs as strings

b_cell_barcodes <- rownames(scd_metadata)[scd_metadata$seurat_clusters %in% b_cell_clusters]

# Check how many cells are in B cell clusters
length(b_cell_barcodes)

# Load the full dataset but only keep B cell clusters
scd <- readRDS("file path")

# Subset only B cell clusters
scd_b_cells <- subset(scd, cells = b_cell_barcodes)

# Save the smaller dataset
saveRDS(scd_b_cells, "file path")

# Check the new dataset
scd_b_cells

write.csv(scd_metadata, "file path")

table(scd_b_cells$seurat_clusters)
DimPlot(scd_b_cells, reduction = "umap", label = TRUE, pt.size = 0.5)
FeaturePlot(scd_b_cells, features = c("Cd4", "Cd8a"), blend = TRUE)
FeaturePlot(scd_b_cells, features = c("Cd19", "Cd8a"), blend = TRUE)



FeaturePlot(scd_b_cells, features = "Cd3e", cols = c("lightgrey", "blue"),pt.size = 0,alpha = 2) + ggtitle("Cd3e")
FeaturePlot(scd_b_cells, features = "Lef1", cols = c("lightgrey", "blue")) + ggtitle("Lef1")
FeaturePlot(scd_b_cells, features = "Cd4", cols = c("lightgrey", "blue")) + ggtitle("Cd4")
FeaturePlot(scd_b_cells, features = "Cd8a", cols = c("lightgrey", "blue")) + ggtitle("Cd8a")
FeaturePlot(scd_b_cells, features = "Cd19", cols = c("lightgrey", "blue")) + ggtitle("Cd19")
FeaturePlot(scd_b_cells, features = "Tbx21", cols = c("lightgrey", "blue")) + ggtitle("Tbx21")
FeaturePlot(scd_b_cells, features = "Gata3", cols = c("lightgrey", "blue")) + ggtitle("Gata3")
FeaturePlot(scd_b_cells, features = "Rorc", cols = c("lightgrey", "blue")) + ggtitle("Rorc")
FeaturePlot(scd_b_cells, features = "foxp3", cols = c("lightgrey", "blue")) + ggtitle("foxp3")
FeaturePlot(scd_b_cells, features = "Ifit3", cols = c("lightgrey", "blue")) + ggtitle("Ifit3")
FeaturePlot(scd_b_cells, features = "Gzmk", cols = c("lightgrey", "blue")) + ggtitle("Gzmk")
FeaturePlot(scd_b_cells, features = "Klrb1", cols = c("lightgrey", "blue")) + ggtitle("Klrb1")
FeaturePlot(scd_b_cells, features = "Pask", cols = c("lightgrey", "blue")) + ggtitle("Pask")
FeaturePlot(scd_b_cells, features = "Prf1", cols = c("lightgrey", "blue")) + ggtitle("Prf1")
FeaturePlot(scd_b_cells, features = "Tox2", cols = c("lightgrey", "blue")) + ggtitle("Tox2")
FeaturePlot(scd_b_cells, features = "Jchain", cols = c("lightgrey", "blue")) + ggtitle("Jchain")
FeaturePlot(scd_b_cells, features = "Ly86", cols = c("lightgrey", "blue")) + ggtitle("Ly86")
FeaturePlot(scd_b_cells, features = "Ms4a1", cols = c("lightgrey", "blue")) + ggtitle("Ms4a1")
FeaturePlot(scd_b_cells, features = "Ighg1", cols = c("lightgrey", "blue")) + ggtitle("Ighg1")
FeaturePlot(scd_b_cells, features = "Ctla4", cols = c("lightgrey", "blue")) + ggtitle("Ctla4")
FeaturePlot(scd_b_cells, features = "Plac8", cols = c("lightgrey", "blue")) + ggtitle("Plac8")
FeaturePlot(scd_b_cells, features = "Cd70", cols = c("lightgrey", "blue")) + ggtitle("Cd70")
FeaturePlot(scd_b_cells, features = "cd52", cols = c("lightgrey", "blue")) + ggtitle("cd52")
FeaturePlot(scd_b_cells, features = "Batf3", cols = c("lightgrey", "blue")) + ggtitle("Batf3")
FeaturePlot(scd_b_cells, features = "Aif1", cols = c("lightgrey", "blue")) + ggtitle("Aif1")
FeaturePlot(scd_b_cells, features = "Tpsab1", cols = c("lightgrey", "blue")) + ggtitle("Tpsab1")
FeaturePlot(scd_b_cells, features = "C1qc", cols = c("lightgrey", "blue")) + ggtitle("C1qc")
FeaturePlot(scd_b_cells, features = "Lyz2", cols = c("lightgrey", "blue")) + ggtitle("Lyz2")
FeaturePlot(scd_b_cells, features = "Ccr7", cols = c("lightgrey", "blue")) + ggtitle("Ccr7")
FeaturePlot(scd_b_cells, features = "Siglech", cols = c("lightgrey", "blue")) + ggtitle("Siglech")
FeaturePlot(scd_b_cells, features = "S100a8", cols = c("lightgrey", "blue")) + ggtitle("S100a8")
FeaturePlot(scd_b_cells, features = "Nkg7", cols = c("lightgrey", "blue")) + ggtitle("Nkg7")
FeaturePlot(scd_b_cells, features = "Ncr", cols = c("lightgrey", "blue")) + ggtitle("Ncr")


# scd_b_cells$cell_type <- "Unknown"
# scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Cd4", ] > 1] <- "CD4+ B Cells"
# scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Cd8a", ] > 1] <- "CD8+ B Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Foxp3", ] > 1] <- "Regulatory T Cells (Tregs)"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Gata3", ] > 1] <- "Th2 Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Tbx21", ] > 1] <- "Th1 Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Rorc", ] > 1] <- "Th17 Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Ms4a1", ] > 1] <- "B Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Nkg7", ] > 1] <- "NK Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["C1qc", ] > 1] <- "Macrophages"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Siglech", ] > 1] <- "Dendritic Cells"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Aif1", ] > 1] <- "Microglia"
# scd_b_cells$cell_type[scd_t_cells@assays$RNA@data["Tpsab1", ] > 1] <- "Mast Cells"

# Load required libraries
library(Seurat)
library(ggplot2)

# Define Output Directory
output_dir <- "file path"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize all cells as "Unknown"
scd_b_cells$cell_type <- "Unknown"

# Assign cell types based on marker expression thresholds
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Cd4", ] > 1] <- "CD4+ B Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Cd8a", ] > 1] <- "CD8+ B Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Foxp3", ] > 1] <- "Regulatory T Cells (Tregs)"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Gata3", ] > 1] <- "Th2 Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Tbx21", ] > 1] <- "Th1 Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Rorc", ] > 1] <- "Th17 Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Ms4a1", ] > 1] <- "B Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Nkg7", ] > 1] <- "NK Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["C1qc", ] > 1] <- "Macrophages"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Siglech", ] > 1] <- "Dendritic Cells"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Aif1", ] > 1] <- "Microglia"
scd_b_cells$cell_type[scd_b_cells@assays$RNA@data["Tpsab1", ] > 1] <- "Mast Cells"

# Define the list of marker genes for visualization
features <- c("Cd3e", "Lef1", "Cd4", "Cd8a", "Tbx21", "Gata3", "Rorc",  
              "Ifit3", "Gzmk","Prf1", "Tox2", "Jchain", "Ly86", 
              "Ms4a1", "Ighg1", "Ctla4", "Plac8", "Cd70", "Batf3", "Aif1", 
              "C1qc", "Lyz2", "Ccr7", "Siglech", "S100a8", "Nkg7")

# Loop through each gene and generate FeaturePlots
for (gene in features) {
  # Generate FeaturePlot with blue color scale
  p <- FeaturePlot(scd_b_cells, features = gene, cols = c("lightgrey", "blue"), 
                   pt.size = 0, alpha = 1) + 
    ggtitle(paste("Expression of", gene))  # Add a title
  
  # Define filename for saving plots
  filename <- paste0(output_dir, "b_cell_subsetting_", gene, ".png")
  
  # Save each FeaturePlot as a high-quality PNG file
  ggsave(filename = filename, plot = p, width = 6, height = 5, dpi = 300)
  
  # Print progress message
  print(paste("Saved:", filename))
}


# Verify counts
table(scd_b_cells$cell_type)

saveRDS(scd_b_cells, "file path")

write.csv(scd_b_cells@meta.data, "file path", row.names = TRUE)



