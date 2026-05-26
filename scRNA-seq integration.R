#Following the Satija lab "Introduction to scRNA-seq integration"
#control & stimulated data, aim: integrate the two conditions together to jointly 
                            #identify cell subpopulations across datasets and explore 
                            #how each group differs across conditions 
library(Seurat)
library(SeuratData) #idk why it worked this time whatever
library(patchwork)
#InstallData("ifnb") #human PBMC from IFN-stim ("stim") cells and control = automatically w SeuratData
ifnb = LoadData("ifnb")

#split RNA measurements into two layers (1 control, 1 stimulated)
ifnb[["RNA"]] = split(ifnb[["RNA"]], f = ifnb$stim)

#analysis without integration (just for fun) --> clusters are defined by both cell type and stimulation = problems later on 
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "unintegrated_clusters"))

#INTEGRATION
#goal: cluster cells from same cell type/subpopulation together 
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE) #a little concerned bc it's taking forever (aka 3 mins lol)

ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]]) # re-join layers after integration
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
DimPlot(ifnb, reduction = "umap", split.by = "stim") #split.by will show the two conditions side by side

#identify the conserved genes in NK cells irrespective of stimulation group, based on p-value
library(BiocManager)
BiocManager::install('mulltest')
install.packages('metap')
Idents(ifnb) <- "seurat_annotations"

nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
        #find the conserved genes in stim/control for the NK cells & their p values 
        #this takes a bit bc of no presto... can't install presto bc it's not available for this version... just wait sorry 
head(nk.markers)

# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
                                                "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()
#view conserved cell type markers across conditions showing both the expression levels and the %cells in a cluster expressing any gene 
