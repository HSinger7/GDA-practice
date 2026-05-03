#SETUP
install.packages('Seurat')
library(Seurat)
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
# Install the remotes package
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")}
#remotes::install_github("satijalab/seurat-data", quiet = TRUE)
#remotes::install_github("satijalab/azimuth", quiet = TRUE)
#remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

#DATA LOADING 
library(dplyr)
library(patchwork)
pbmc.data <- Read10X(data.dir = "~/Desktop/Personal/For funsies/Data/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#pbmc   ("13714 features across 2700 samples within 1 (RNA) assay)   
#examine first rows --> pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30])]" will show the 3x30 sparse Matrix class "dgCMatrix"
dense.size = object.size(as.matrix(pbmc.data)) #if you run dense.size = 709591472 bytes
sparse.size = object.size(pbmc.data) #29905192 bytes 
dense.size/sparse.size #23.7 --> done to save space basically 

#PRE-PROCESSING scRNA-seq 
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern="^MT-")
VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter = feature-feature interactions 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #remove unwanted cells! 

#NORMALIZE DATA 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) 
#normalize expression of each cell by the total expression * 10,000 and then log transform it 
#could also be achieved by just doing pbmc <- NormalizeData(pbmc) without specifying the parameters 
#but this one ^ assumes that all cells have the same number of RNA molecules (not always true)

#FEATURE IDENTIFICATION
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #default 2000 features per dataset

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc) #without labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #label the top 10 most variable ones from plot 1
plot1 + plot2

#SCALE DATA (PREPROCESSING TO PCA)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

#PCA 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #visualize the top genes associated with reduction components, PC_1 & PC_2
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #only show PC_1 
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) #show them all 

#alternative to Dim stuff (less computationally intensive)
ElbowPlot(pbmc) #elbow is at ~9-10 so the majority of the true signal is in the first 10 PCs 

#clustering - edo to xanoume ligo 
pbmc <- FindNeighbors(pbmc, dims = 1:10) #constructs KNN graph based on Jaccard similarity for first 10 PCs 
pbmc <- FindClusters(pbmc, resolution = 0.5) #resolution sets the granularity of downstream clustering, higher res = more clusters / if ~3K cells use 0.4-1.2 
head(Idents(pbmc), 5) #Finds cluster IDs of the first 5 cells

#Non-linear dimensional reduction - UMAP (see notes for limits)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE) #label = TRUE puts the labels into the clusters / can remove

#Finding differentially expressed features (markers)
#find all markers in cluster 2 
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2) #ident.1 specifies which cluster to look at, e.g. 2 
head(cluster2.markers, n = 5) 

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

#find markers for every cluster compared to all remaining cells & report only the +ve
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE) #takes a while bc calculates all clusters 
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

VlnPlot(pbmc, features = c("MS4A1", "CD79A")) #expression probability across clusters 
VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE) #plot raw counts as well 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
    #umap plots for these genes 

pbmc.markers %>% #expression heatmap for given cells & features
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

