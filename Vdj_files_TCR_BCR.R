library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(SeuratWrappers)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ineq)
library(ComplexHeatmap)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")


# Load Seurat object
seurat_obj <- readRDS("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/RDS files/MRV R848 scRNAseq cellcycle_corrected.rds")

# Load VDJ data for T and B cells
t_vdj <- read.csv("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/vdj_t/filtered_contig_annotations.csv")
# b_vdj <- read.csv("/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/vdj_b/filtered_contig_annotations.csv")
head(t_vdj)
# Subset T and B cells
T_clusters <- c(2,6,15,17)
# B_clusters <- c(0,1,4,28,29)

T_cells <- subset(seurat_obj, idents = T_clusters)
# B_cells <- subset(seurat_obj, idents = B_clusters)

# Recluster T cells
T_cells <- ScaleData(T_cells)
T_cells <- RunPCA(T_cells)
T_cells <- FindNeighbors(T_cells, dims = 1:10)
T_cells <- FindClusters(T_cells, resolution = 0.5)
T_cells <- RunUMAP(T_cells, dims = 1:10)

# # Recluster B cells
# B_cells <- ScaleData(B_cells)
# B_cells <- RunPCA(B_cells)
# B_cells <- FindNeighbors(B_cells, dims = 1:10)
# B_cells <- FindClusters(B_cells, resolution = 0.5)
# B_cells <- RunUMAP(B_cells, dims = 1:10)


# Define the output directory to save plots
output_dir <- "/Users/navyasreechenchu/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/Plots and Images/"
dir.create(output_dir, showWarnings = FALSE)

# Feature plots for T cells (save each plot individually)
t_genes <- c("Cd3e", "Tcf7", "Lef1", "Sell", "Cd4", "Il2ra", "Il7r", "Cd8a", "Cd8b1", "Tbx21", "Ifng", "Gata3", "Runx2", "Rorc", "Anxa1", "Foxp3", "Il32", "Ifit3", "Mx1", "Gzmk", "Ccl5", "Klrb1", "Pask", "Prf1", "Gzma", "Tox2", "Foxb", "Cd7", "Ncr", "Nkg7", "Ctla4", "Prf1", "Xcl1", "Mki67")

# Loop through each gene, generate feature plot and save it
for (gene in t_genes) {
  tryCatch({
    p <- FeaturePlot(T_cells, features = gene) + ggtitle(paste("T cells _", gene))
    
    # Define the file path for saving the plot
    plot_filename <- paste(output_dir, "T_cell_reclustering_", gene, ".png", sep = "")
    
    # Save the plot
    ggsave(plot_filename, plot = p, width = 8, height = 6, dpi = 300)
  }, error = function(e) {
    message(paste("Skipping gene:", gene, "due to error:", e$message))
  })
}



# UMAP plots
DimPlot(T_cells, group.by = "seurat_clusters") + ggtitle("T cells UMAP")
# DimPlot(B_cells, group.by = "seurat_clusters") + ggtitle("B cells UMAP")

# # Feature plots for T cells
# t_genes <- c("Cd3e", "Tcf7", "Lef1", "Sell", "Cd4", "Il2ra", "Il7r", "Cd8a", "Cd8b1", "Tbx21", "Ifng", "Gata3", "Runx2", "Rorc", "Ccr6", "Anxa1", "Foxp3", "Il32", "Ifit3", "Mx1", "Gzmk", "Ccl5", "Klrb1", "Pask", "Prf1", "Gzma", "Tox2", "Foxb", "Cd7", "Ncr", "Nkg7", "Ctla4", "Prf1", "Xcl1", "Mki67")
# # FeaturePlot(T_cells, features = t_genes)
# for (gene in t_genes){
#   p <- FeaturePlot(T_cells, features = gene) + ggtitle(paste("T cells _",gene))
#   print(p)
# }
# # Feature plots for B cells
# b_genes <- c("Cd19", "Cd22", "Cd79a", "Igha1", "Ly86", "Ms4a1", "Ighd", "Ighg1", "Mzb1", "Vpreb3", "Cd79b", "Ctla4", "Ptprj", "Bhlhe41", "Zbtb20", "Plac8", "Jchain", "Irf4", "Cd70", "Cd38", "Cxcr3", "Mki67")
# # FeaturePlot(B_cells, features = b_genes)
# for (gene in b_genes){
#   p <- FeaturePlot(B_cells, features = gene) + ggtitle(paste("B cells _",gene))
#   print(p)
# }
# DEG and top10 heatmap for T cells
T_markers <- FindAllMarkers(T_cells)
T_top10 <- T_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(T_cells, features = unique(T_top10$gene))
head(T_markers)
# # DEG and top10 heatmap for B cells
# B_markers <- FindAllMarkers(B_cells)
# B_top10 <- B_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
# DoHeatmap(B_cells, features = unique(B_top10$gene))

# GO pathway analysis example for T cells
cluster_deg_T <- subset(T_markers, cluster == 1 & p_val_adj < 0.05)
cluster_genes_T <- cluster_deg_T$gene
ego_T <- enrichGO(gene = cluster_genes_T, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")
barplot(ego_T)

# # GO pathway analysis example for B cells
# cluster_deg_B <- subset(B_markers, cluster == 1 & p_val_adj < 0.05)
# cluster_genes_B <- cluster_deg_B$gene
# ego_B <- enrichGO(gene = cluster_genes_B, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")
# barplot(ego_B)

# TCR clone analysis (from t_vdj)
tcr_clones <- table(t_vdj$clonotype_id)
ineq(as.numeric(tcr_clones), type = "Gini")

# # BCR clone analysis (from b_vdj)
# bcr_clones <- table(b_vdj$clonotype_id)
# ineq(as.numeric(bcr_clones), type = "Gini")

# Plot composite frequencies by cluster
T_freq <- prop.table(table(Idents(T_cells)))
barplot(T_freq, main = "T Cell Cluster Frequency")

# B_freq <- prop.table(table(Idents(B_cells)))
# barplot(B_freq, main = "B Cell Cluster Frequency")

# Plot DimPlots for TCR clones with >3 occurrences
high_clone_tcr <- names(tcr_clones[tcr_clones > 3])
T_cells$high_clone <- ifelse(T_cells$raw_clonotype_id %in% high_clone_tcr, T_cells$clonotype_id, "Other")
DimPlot(T_cells, group.by = "high_clone")

# # Plot DimPlots for BCR clones with >3 occurrences
# high_clone_bcr <- names(bcr_clones[bcr_clones > 3])
# B_cells$high_clone <- ifelse(B_cells$clonotype_id %in% high_clone_bcr, B_cells$clonotype_id, "Other")
# DimPlot(B_cells, group.by = "high_clone")






