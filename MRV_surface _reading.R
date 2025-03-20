library(Seurat)
setwd("/Volumes/bigley/Active/scRNAseq_DataSets/chenchu/scripts/R_Scripts/")

library(Seurat)

# Load the Seurat object
scd <- readRDS("~/Documents/Bigley Lab/Projects/MRV Lupus Model/MRV Lupus scRNAseq/cellcycle_corrected.rds")

# Check if the file loaded correctly
print(scd)

# Subset the object to include only the clusters of interest
clusters_of_interest <- c(2, 6, 15, 17)
scd_subset <- subset(scd, idents = clusters_of_interest)

# Check the subset
table(scd_subset$seurat_clusters)
# Example markers (you might need to adjust based on your actual dataset)
cd4_markers <- c("Cd4")
cd8_markers <- c("Cd8a")

# Check expression of CD4 and CD8 markers
FeaturePlot(scd_subset, features = cd4_markers)
FeaturePlot(scd_subset, features = cd8_markers)

# Annotate cells based on markers
scd_subset$cell_type <- "Unknown"
scd_subset$cell_type[scd_subset$RNA$Cd4 > 0] <- "CD4+ T Cells"
scd_subset$cell_type[scd_subset$RNA$Cd8a > 0] <- "CD8+ T Cells"

# Verify annotations
table(scd_subset$cell_type)
# Fetch TCR and BCR sequences (adjust if they are stored differently)
tcr_table <- FetchData(scd, vars = "TCR_sequences") # Replace with actual TCR column
bcr_table <- FetchData(scd, vars = "BCR_sequences") # Replace with actual BCR column

# Save to CSV
write.csv(tcr_table, "path/to/tcr_table.csv")
write.csv(bcr_table, "path/to/bcr_table.csv")
# Define Gini coefficient function
gini_coefficient <- function(x) {
  sorted_x <- sort(x)
  n <- length(x)
  cum_x <- cumsum(sorted_x)
  gini <- (2 * sum(cum_x) / (n * sum(x))) - ((n + 1) / n)
  return(gini)
}

# Calculate Gini coefficient for TCR and BCR
tcr_diversity <- apply(tcr_table, 2, gini_coefficient)
bcr_diversity <- apply(bcr_table, 2, gini_coefficient)

# Plot Gini coefficients
barplot(tcr_diversity, main = "TCR Diversity", ylab = "Gini Coefficient")
barplot(bcr_diversity, main = "BCR Diversity", ylab = "Gini Coefficient")
