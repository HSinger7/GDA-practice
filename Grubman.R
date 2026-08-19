getwd()
workingDir= "/Users/lenyasinger/Desktop/Personal/For funsies/GDA-practice/Grubman replication"
setwd(workingDir)
library(Seurat)

#Section 1: Data importing  
#data downloaded from GSE138852 
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

install.packages("BRETIGEA") #the one Grubman used 
library("BRETIGEA")
library("dplyr")
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
#add meta data to compare with Grubman original umap
covariates = read.csv("~/Desktop/Personal/For funsies/GDA-practice/GSE138852_covariates.csv")
rownames (covariates) = covariates$X
seurat_obj = AddMetaData(seurat_obj, metadata = covariates)

#okay now compare them 
library(patchwork)
library(ggplot2)
H_annotation = DimPlot(seurat_obj,
                       group.by = "bretigea_celltype",
                       label = TRUE,
                       repel = TRUE) +
  ggtitle("Heleni (BRETIGEA)") 

# Grubman's original annotations
G_annotation <- DimPlot(seurat_obj, 
                        group.by = "oupSample.cellType",
                        label = TRUE, 
                        repel = TRUE) +
  ggtitle("Grubman et al.")

H_annotation + G_annotation #=<3

#create UMAP for AD vs control 
DimPlot(seurat_obj,
        group.by = "oupSample.subclustCond",
        label = FALSE, 
        cols = c("AD" = "#7B2D8B", 
                 "ct" = "#2D8B2D",
                 "undetermined" = "#696969")) +
  ggtitle("AD vs Control cell type composition") +
  scale_color_manual(values = c(AD = "#7B2D8B", ct = "#2D8B2D", undetermined = "#696969"),
                     labels = c(AD = "Alzheimer's Disease", ct = "Control", undetermined = "Unknown")
  ) +
  labs(color = "Condition") + 
  theme(legend.title = element_text(face = "bold"))



