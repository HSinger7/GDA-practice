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
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
