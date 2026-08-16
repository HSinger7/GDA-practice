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

#calculate % mito & visualize QC matrix w/volcano plot
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

#Section 3: filter low quality cells out 
seurat_obj = subset(seurat_obj,
                    subset = nFeature_RNA > 500 &
                      nFeature_RNA < 7000 &
                      percent.mt < 15)
dim(seurat_obj)  # check how many cells remain = 10850, 9866 

#Section 4: Normalize, variable features, scale 
seurat_obj = NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj = ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))

#Section 5: dimensionality reduction 
seurat_obj = RunPCA(seurat_obj, npcs = 30)
ElbowPlot(seurat_obj, ndims = 30)
seurat_obj = FindNeighbors(seurat_obj, dims = 1:15) #plateaus around 15 
seurat_obj = FindClusters(seurat_obj, resolution = 0.5) #nodes = 9866, edges = 366331, 
#max modularity in 10 random starts = 0.9212
#no of communities: 15
seurat_obj = RunUMAP(seurat_obj, dims = 1:15)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

#Section 5: Cell type annotation 
markers = FindAllMarkers(seurat_obj, #calculating 14 clusters = takes a while
                         only.pos = TRUE, 
                         min.pct = 0.25, 
                         logfc.threshold = 0.25)

#Grubman used Bretigea to anotate the cells 

install.packages("BRETIGEA")
library("BRETIGEA")
expr_matrix <- as.matrix(GetAssayData(seurat_obj, layer = "data"))
cell_types = brainCells(expr_matrix, nMarker = 50) #idk about this 50 
seurat_obj <- AddMetaData(seurat_obj, 
                          metadata = as.data.frame(cell_types))

#on the UMAP, see which cell type is where  
FeaturePlot(seurat_obj, 
            features = c("ast", "end", "mic", "neu", "oli", "opc", ""),
            ncol = 3)

cell_type_labels <- colnames(cell_types)[apply(cell_types, 
                                               1, 
                                               which.max)]
seurat_obj$bretigea_celltype <- cell_type_labels
seurat_obj$bretigea_celltype <- recode(seurat_obj$bretigea_celltype,
                                       "ast" = "Astrocyte",
                                       "end" = "Endothelial",
                                       "mic" = "Microglia",
                                       "neu" = "Neuron",
                                       "oli" = "Oligodendrocyte",
                                       "opc" = "OPC"
)

#they also have hybrid and unidentified cells... decide if you want to add them and figure that out lol 

library(ggplot2)
DimPlot(seurat_obj,
        group.by = "bretigea_celltype",
        label = TRUE,
        repel = TRUE,
        cols = c(
          "Astrocyte" = "#FFB6C1",      
          "Microglia" = "#4169E1",       
          "Neuron" = "#FF0000",         
          "Oligodendrocyte" = "#FFA500",
          "OPC" = "#9370DB",             
          "Endothelial" = "#8B4513" #,     
          #"Hybrid" = "#000000", 
          #"Unidentified" = "#808080"
        )) +
  ggtitle("Cell type annotation (BRETIGEA)")


#Section 6: add metadata to see AD vs ctrl 
Idents(seurat_obj) <- "celltype"  # after you've annotated clusters

# DE within microglia for example
micro_de <- FindMarkers(seurat_obj,
                        ident.1 = "AD",
                        ident.2 = "control",
                        group.by = "condition",
                        subset.ident = "microglia")
