library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

##load 87_B1 data (D2 stroke)
s87_B1 <- Load10X_Spatial(
  "/Users/erica_1/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-087-B1/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "87_B1",
  slice = "D2_strokeB_",
  filter.matrix = TRUE,
  to.upper = FALSE,
  use.names = TRUE, unique.features = TRUE)

plot1 <- VlnPlot(s87_B1, features = "nCount_87_B1", pt.size = 0.1) + NoLegend()
#yes to negative binomial
plot2 <- SpatialFeaturePlot(s87_B1, features = "nCount_87_B1") + theme(legend.position = "right")
#normalize
s87_B1 <- SCTransform(s87_B1, assay = "87_B1", verbose = FALSE)

#perform clustering with UMAP (Preferred)
s87_B1 <- RunPCA(s87_B1, assay = "SCT", verbose = FALSE)
s87_B1 <- FindNeighbors(s87_B1, reduction = "pca", dims = 1:30)
s87_B1 <- FindClusters(s87_B1, verbose = FALSE)
s87_B1 <- RunUMAP(s87_B1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(s87_B1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(s87_B1,label = TRUE, label.size = 2, alpha = c(0.1, 1))
p1 + p2

#87_D1=D10 stroke
s87_D1 <- Load10X_Spatial(
  "/Users/erica_1/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-087-D1/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "87_D1",
  slice = "D10_strokeD_",
  filter.matrix = TRUE,
  to.upper = FALSE,
  use.names = TRUE, unique.features = TRUE)

plot1 <- VlnPlot(s87_D1, features = "nCount_87_D1", pt.size = 0.1) + NoLegend()
#yes to negative binomial
plot2 <- SpatialFeaturePlot(s87_D1, features = "nCount_87_D1") + theme(legend.position = "right")
#normalize
s87_D1 <- SCTransform(s87_D1, assay = "87_D1", verbose = FALSE)

#cluster with UMAP
s87_D1 <- RunPCA(s87_D1, assay = "SCT", verbose = FALSE)
s87_D1 <- FindNeighbors(s87_D1, reduction = "pca", dims = 1:30)
s87_D1 <- FindClusters(s87_D1, verbose = FALSE)
s87_D1 <- RunUMAP(s87_D1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(s87_D1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(s87_D1,label = TRUE, label.size = 2, alpha = c(0.1, 1))
p1 + p2

#88_B1= sham (d2)
s88_B1 <- Load10X_Spatial(
  "/Users/erica_1/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-088-B1/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "88_B1",
  slice = "D2_sham_",
  filter.matrix = TRUE,
  to.upper = FALSE,
  use.names = TRUE, unique.features = TRUE)
plot1 <- VlnPlot(s88_B1, features = "nCount_88_B1", pt.size = 0.1) + NoLegend()
#yes to negative binomial
plot2 <- SpatialFeaturePlot(s88_B1, features = "nCount_88_B1") + theme(legend.position = "right")
#normalize
s88_B1 <- SCTransform(s88_B1, assay = "88_B1", verbose = FALSE)
#plot inidividual genes
plot_ALDOC <- SpatialFeaturePlot(s88_B1, features = c("Aldoc"), alpha = c(0.1, 1))
genes_s88_B1<-rownames(x = s88_B1)
setwd('/Users/erica/Desktop/Faiz/Visium/')
write.table(genes_s88_B1, "D2sham_genelist.txt")
#perform clustering with UMAP (Preferred)
s88_B1 <- RunPCA(s88_B1, assay = "SCT", verbose = FALSE)
s88_B1 <- FindNeighbors(s88_B1, reduction = "pca", dims = 1:30)
s88_B1 <- FindClusters(s88_B1, verbose = FALSE)
s88_B1 <- RunUMAP(s88_B1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(s88_B1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(s88_B1,label = TRUE, label.size = 2, alpha = c(0.1, 1))
p1 + p2

#88_C1=D21 stroke
s88_C1 <- Load10X_Spatial(
  "/Users/erica_1/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-088-C1/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "88_C1",
  slice = "D21_stroke_",
  filter.matrix = TRUE,
  to.upper = FALSE,
  use.names = TRUE, unique.features = TRUE)

plot1 <- VlnPlot(s88_C1, features = "nCount_88_C1", pt.size = 0.1) + NoLegend()
#yes to negative binomial, but has an odd second caudal hump to max peak
plot2 <- SpatialFeaturePlot(s88_C1, features = "nCount_88_C1") + theme(legend.position = "right")
#normalize
s88_C1 <- SCTransform(s88_C1, assay = "88_C1", verbose = FALSE)

#perform clustering with UMAP (Preferred)
s88_C1 <- RunPCA(s88_C1, assay = "SCT", verbose = FALSE)
s88_C1 <- FindNeighbors(s88_C1, reduction = "pca", dims = 1:30)
s88_C1 <- FindClusters(s88_C1, verbose = FALSE)
s88_C1 <- RunUMAP(s88_C1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(s88_C1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(s88_C1,label = TRUE, label.size = 2, alpha = c(0.1, 1))
p1 + p2


#spatially variable features
de_markers_s88_C1 <- FindMarkers(s88_C1, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = s88_C1, features = rownames(de_markers_s88_C1)[1:9], alpha = c(0.1, 1), ncol = 3)

#d21 new de spatial method
s88_C1 <- FindSpatiallyVariableFeatures(s88_C1, assay = "SCT", features = VariableFeatures(s88_C1)[1:1000], 
                                       selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(s88_C1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(s88_C1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#d2 sham
s88_B1 <- FindSpatiallyVariableFeatures(s88_B1, assay = "SCT", features = VariableFeatures(s88_B1)[1:1000], 
                                        selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(s88_B1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(s88_B1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#d10
s87_D1 <- FindSpatiallyVariableFeatures(s87_D1, assay = "SCT", features = VariableFeatures(s87_D1)[1:1000], 
                                        selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(s87_D1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(s87_D1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#d2
s87_B1 <- FindSpatiallyVariableFeatures(s87_B1, assay = "SCT", features = VariableFeatures(s87_B1)[1:1000], 
                                        selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(s87_B1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(s87_B1, features = top.features, ncol = 3, alpha = c(0.1, 1))
