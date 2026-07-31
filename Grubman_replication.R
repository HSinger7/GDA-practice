getwd()
workingDir= "/Users/lenyasinger/Desktop/Personal/For funsies/GDA-practice/Grubman replication"
setwd(workingDir)
library(Seurat)

#Section 1: Data importing  
#data downloaded from GSE138852 and unzipped in terminal with "gunzip" and file name 
G_data = "GSE138852_counts.csv"
G_data=read.csv(G_data, row.names = 1, check.names = FALSE)
dim(G_data) #10850, 13214 damn that's a lot of data lol 
#View(G_data) (rows = genes, columns = cells)

#Section 2: standard filtering  
seurat_obj= CreateSeuratObject(
  counts = G_data, 
  min.cells = 3, #keep genes detected in at least 3 cells / can be changed for 5 --> gene-level filter
  min.features = 200 #keep cells with 200+ genes detected --> cell-level filter
)

#calculate % mito & visualize volcano plot
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
