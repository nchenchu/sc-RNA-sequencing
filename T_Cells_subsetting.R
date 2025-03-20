data <- readRDS("/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV R848 scRNAseq cellcycle_corrected.rds")

library(Seurat)

# Open the RDS file and load only metadata
con <- gzfile("/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV R848 scRNAseq cellcycle_corrected.rds", "rb")

scd_metadata <- readRDS(con)@meta.data

# Check metadata structure
str(scd_metadata)

# Extract barcodes of T cell clusters
t_cell_clusters <- c("2", "6", "15", "17")  # Cluster IDs as strings

t_cell_barcodes <- rownames(scd_metadata)[scd_metadata$seurat_clusters %in% t_cell_clusters]

# Check how many cells are in T cell clusters
length(t_cell_barcodes)                         

# Load the full dataset but only keep T cell clusters
scd <- readRDS("/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV R848 scRNAseq cellcycle_corrected.rds")

# Subset only T cell clusters
scd_t_cells <- subset(scd, cells = t_cell_barcodes)

write.csv(t_cell_barcodes,"/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/t_cell_barcodes.csv")

# Save the smaller dataset
saveRDS(scd_t_cells, "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV_T_Cells.rds")

# Check the new dataset
scd_t_cells

write.csv(scd_metadata, "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV_metadata.csv")

table(scd_t_cells$seurat_clusters)
DimPlot(scd_t_cells, reduction = "umap", label = TRUE, pt.size = 0.5)
FeaturePlot(scd_t_cells, features = c("Cd4", "Cd8a"), blend = TRUE)
FeaturePlot(scd_t_cells, features = c("Cd3e", "Lef1"), blend = TRUE)
#FeaturePlot(scd_t_cells, features = c("Cd3e", "Lef1","Cd4","Cd8a","Cd3e","Tbx21","Gata3","Rorc","foxp3","Ifit3","Gzmk","Klrb1","Pask","Prf1","Tox2","Jchain","Ly86","Ms4a1","Ighg1","Ctla4","Plac8","Cd70","cd52","Batf3","Aif1","Tpsab1","C1qc","Lyz2","Ccr7","Siglech","S100a8","Nkg7","Ncr"))

# Load required library
library(Seurat)
library(ggplot2)

FeaturePlot(scd_t_cells, features = "Cd3e", cols = c("lightgrey", "blue"),pt.size = 0,alpha = 2) + ggtitle("Cd3e")
FeaturePlot(scd_t_cells, features = "Lef1", cols = c("lightgrey", "blue")) + ggtitle("Lef1")
FeaturePlot(scd_t_cells, features = "Cd4", cols = c("lightgrey", "blue")) + ggtitle("Cd4")
FeaturePlot(scd_t_cells, features = "Cd8a", cols = c("lightgrey", "blue")) + ggtitle("Cd8a")
FeaturePlot(scd_t_cells, features = "Tbx21", cols = c("lightgrey", "blue")) + ggtitle("Tbx21")
FeaturePlot(scd_t_cells, features = "Gata3", cols = c("lightgrey", "blue")) + ggtitle("Gata3")
FeaturePlot(scd_t_cells, features = "Rorc", cols = c("lightgrey", "blue")) + ggtitle("Rorc")
FeaturePlot(scd_t_cells, features = "foxp3", cols = c("lightgrey", "blue")) + ggtitle("foxp3")
FeaturePlot(scd_t_cells, features = "Ifit3", cols = c("lightgrey", "blue")) + ggtitle("Ifit3")
FeaturePlot(scd_t_cells, features = "Gzmk", cols = c("lightgrey", "blue")) + ggtitle("Gzmk")
FeaturePlot(scd_t_cells, features = "Klrb1", cols = c("lightgrey", "blue")) + ggtitle("Klrb1")
FeaturePlot(scd_t_cells, features = "Pask", cols = c("lightgrey", "blue")) + ggtitle("Pask")
FeaturePlot(scd_t_cells, features = "Prf1", cols = c("lightgrey", "blue")) + ggtitle("Prf1")
FeaturePlot(scd_t_cells, features = "Tox2", cols = c("lightgrey", "blue")) + ggtitle("Tox2")
FeaturePlot(scd_t_cells, features = "Jchain", cols = c("lightgrey", "blue")) + ggtitle("Jchain")
FeaturePlot(scd_t_cells, features = "Ly86", cols = c("lightgrey", "blue")) + ggtitle("Ly86")
FeaturePlot(scd_t_cells, features = "Ms4a1", cols = c("lightgrey", "blue")) + ggtitle("Ms4a1")
FeaturePlot(scd_t_cells, features = "Ighg1", cols = c("lightgrey", "blue")) + ggtitle("Ighg1")
FeaturePlot(scd_t_cells, features = "Ctla4", cols = c("lightgrey", "blue")) + ggtitle("Ctla4")
FeaturePlot(scd_t_cells, features = "Plac8", cols = c("lightgrey", "blue")) + ggtitle("Plac8")
FeaturePlot(scd_t_cells, features = "Cd70", cols = c("lightgrey", "blue")) + ggtitle("Cd70")
FeaturePlot(scd_t_cells, features = "cd52", cols = c("lightgrey", "blue")) + ggtitle("cd52")
FeaturePlot(scd_t_cells, features = "Batf3", cols = c("lightgrey", "blue")) + ggtitle("Batf3")
FeaturePlot(scd_t_cells, features = "Aif1", cols = c("lightgrey", "blue")) + ggtitle("Aif1")
FeaturePlot(scd_t_cells, features = "Tpsab1", cols = c("lightgrey", "blue")) + ggtitle("Tpsab1")
FeaturePlot(scd_t_cells, features = "C1qc", cols = c("lightgrey", "blue")) + ggtitle("C1qc")
FeaturePlot(scd_t_cells, features = "Lyz2", cols = c("lightgrey", "blue")) + ggtitle("Lyz2")
FeaturePlot(scd_t_cells, features = "Ccr7", cols = c("lightgrey", "blue")) + ggtitle("Ccr7")
FeaturePlot(scd_t_cells, features = "Siglech", cols = c("lightgrey", "blue")) + ggtitle("Siglech")
FeaturePlot(scd_t_cells, features = "S100a8", cols = c("lightgrey", "blue")) + ggtitle("S100a8")
FeaturePlot(scd_t_cells, features = "Nkg7", cols = c("lightgrey", "blue")) + ggtitle("Nkg7")
FeaturePlot(scd_t_cells, features = "Ncr", cols = c("lightgrey", "blue")) + ggtitle("Ncr")




scd_t_cells$cell_type <- "Unknown"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Cd4", ] > 1] <- "CD4+ T Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Cd8a", ] > 1] <- "CD8+ T Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Foxp3", ] > 1] <- "Regulatory T Cells (Tregs)"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Gata3", ] > 1] <- "Th2 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Tbx21", ] > 1] <- "Th1 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Rorc", ] > 1] <- "Th17 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Ms4a1", ] > 1] <- "B Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Nkg7", ] > 1] <- "NK Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["C1qc", ] > 1] <- "Macrophages"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Siglech", ] > 1] <- "Dendritic Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Aif1", ] > 1] <- "Microglia"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Tpsab1", ] > 1] <- "Mast Cells"

# Load required libraries
library(Seurat)
library(ggplot2)

# Define Output Directory
output_dir <- "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize all cells as "Unknown"
scd_t_cells$cell_type <- "Unknown"

# Assign cell types based on marker expression thresholds
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Cd4", ] > 1] <- "CD4+ T Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Cd8a", ] > 1] <- "CD8+ T Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Foxp3", ] > 1] <- "Regulatory T Cells (Tregs)"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Gata3", ] > 1] <- "Th2 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Tbx21", ] > 1] <- "Th1 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Rorc", ] > 1] <- "Th17 Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Ms4a1", ] > 1] <- "B Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Nkg7", ] > 1] <- "NK Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["C1qc", ] > 1] <- "Macrophages"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Siglech", ] > 1] <- "Dendritic Cells"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Aif1", ] > 1] <- "Microglia"
scd_t_cells$cell_type[scd_t_cells@assays$RNA@data["Tpsab1", ] > 1] <- "Mast Cells"

# Define the list of marker genes for visualization
features <- c("Cd3e", "Lef1", "Cd4", "Cd8a", "Tbx21", "Gata3", "Rorc",  
              "Ifit3", "Gzmk","Prf1", "Tox2", "Jchain", "Ly86","Ms4a1",
              "Ighg1", "Ctla4", "Plac8", "Cd70", "Batf3", "Aif1", 
              "C1qc", "Lyz2", "Ccr7", "Siglech", "S100a8", "Nkg7")

# Loop through each gene and generate FeaturePlots
for (gene in features) {
  # Generate FeaturePlot with blue color scale
  p <- FeaturePlot(scd_t_cells, features = gene, cols = c("lightgrey", "blue"), 
                   pt.size = 0, alpha = 1) + 
    ggtitle(paste("Expression of", gene))  # Add a title
  
  # Define filename for saving plots
  filename <- paste0(output_dir, "t_cell_subsetting_", gene, ".png")
  
  # Save each FeaturePlot as a high-quality PNG file
  ggsave(filename = filename, plot = p, width = 6, height = 5, dpi = 300)
  
  # Print progress message
  print(paste("Saved:", filename))
}






# Verify counts
table(scd_t_cells$cell_type)

saveRDS(scd_t_cells, "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/MRV_T_Cells_Annotated.rds")

write.csv(scd_t_cells@meta.data, "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/T_Cells.csv", row.names = TRUE)
