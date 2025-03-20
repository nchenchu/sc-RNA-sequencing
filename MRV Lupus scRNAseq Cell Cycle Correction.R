install.packages("Seurat")

install.packages("reremotes::install_github("satijalab/seurat", ref = "develop")

getwd()

list.files("/storage1.ris.wustl.edu/bigley/Active/scRNAseq_DataSets/chenchu/scripts/")


library(Seurat)

setwd("/storage1.ris.wustl.edu/bigley/Active/scRNAseq_DataSets/chenchu/scripts/R_Scripts")

# Update with the new path to the file
file_path <- "C:/Users/chenchu/Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Chenchu MRV D0 R848/cellcycle_corrected.rds"

# Check if the file exists at the new path
file.exists(file_path)

# Attempt to load the file again
scd <- readRDS(file_path)

# If the file exists, try reading it again
scd <- readRDS("~/Documents/Bigley Lab/Projects/MRV Lupus Model/MRV Lupus scRNAseq/cellcycle_corrected.rds")


scd <- readRDS("~/Documents/Bigley Lab/Projects/MRV Lupus Model/MRV Lupus scRNAseq/cellcycle_corrected.rds")

DefaultAssay(scd) <-'SCT'
libaryFeaturePlot(scd,features = 'Foxp3')

DefaultAssay(scd) <-'virus'
FeaturePlot(scd,features = 'ORF69')
WhichCells(scd, features = 'ORF69')

DefaultAssay(scd) <-'virus'
FeaturePlot(scd,features = 'Foxp3', split.by = 'Infection') 

DefaultAssay(scd) <-'SCT'
FeaturePlot(scd,features = 'ORF73',split.by = 'Infection')

# CLUSTERS by ORF expression copy number
DefaultAssay(scd) <-'virus'
highlight <- WhichCells(scd , expression = ORF29 > 10)
UMAPPlot( scd , cells.highlight = highlight )

# FIRST annotate the object

# SUBSET the object on TECs
tec.sub <- subset( scd , )
# CALCULATE DIFFERENTIAL GENE EXPRESSION
markers <- FindAllMarkers( tec.sub )
# HEAT MAP
DoHeatmap( tec.sub , features = markers$gene )

# NEW CLUSTERS and Dimplots and Violing Plots of ORFs
DefaultAssay(scd) <-'SCT'
new.cluster.ids <- c("DP","DP","DP","DP","DP","DP","DP","DP","CD8_SP_mature","DN","DP","DP", "DP", "DP", "DP", "DP","CD8_SP","DN","CD4_SP", "DP","DP","DP","DP","DP","CD8_SP", "CD8_SP","DN","DN","DP","DP","DP", "DP","M_&_B_cells","DP","DP","RBCs", "CD8_SP","DN","DP","cTEC","DP", "DN", "ILC_NK","Treg","DCs","Mesenchyme","unknown","mTEC","Fibroblasts")
names(new.cluster.ids) <- levels(scd)
scd <- RenameIdents(scd, new.cluster.ids)
DimPlot(scd, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(scd, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = 'Infection') 
DimPlot(scd, reduction = "umap", pt.size = 0.5, split.by = 'Infection') 
VlnPlot(scd, features = 'ORF39')
VlnPlot(scd, features = 'ORF39', split.by = 'Infection')
DoHeatmap(scd, features = VariableFeatures(scd)[1:100], cells = 1:500, size = 4,
          angle = 90) + NoLegend()
DoHeatmap(scd, sub("DP"), features = markers$gene )
barplot.default(scd, group.by = 'seurat_clusters', split.by = 'Infection')

DefaultAssay(scd) <-'virus'
VlnPlot(scd, features = 'ORF102')
VlnPlot(scd, features = 'ORF102', split.by = 'Infection') 
VlnPlot(scd, features = 'ORF102') + scale_y_continuous(limits = c(0,20))
VlnPlot(scd, features = 'ORF102', split.by = 'Infection') + scale_y_continuous(limits = c(0,3))

DefaultAssay(scd) <-'SCT'
FeaturePlot(scd,features = 'Cd8a')

DefaultAssay(scd) <-'SCT'
FeaturePlot(scd,features = 'ORF102', split.by = 'Infection')

DefaultAssay(scd) <-'virus'
VlnPlot(scd, features = 'ORF102')
VlnPlot(scd, features = 'ORF69', split.by = 'Infection') 
FeaturePlot(scd,features = 'Cd8a', min.cutoff = 0)

DimPlot(scd, reduction = "umap", label = TRUE, pt.size = 0.5) 

DimPlot(scd, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = 'Affstat or Sample or Condition') 

