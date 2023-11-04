library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(cowplot)
library(wesanderson)

New.stroke <- Read10X(data.dir = "./Stroke/outs/filtered_feature_bc_matrix")
New_stroke <- CreateSeuratObject(counts = New.stroke, project = "Stroke", min.cells = 3, min.features = 200)

New.Naive <- Read10X(data.dir = "./Naive/outs/filtered_feature_bc_matrix")
New_Naive <- CreateSeuratObject(counts = New.Naive, project = "Naive", min.cells = 3, min.features = 200)

New_S_N <- merge(New_stroke, y = New_Naive, add.cell.ids = c("Stroke", "Naive"), project = "New_S_N")
New_S_N
#22527 features across 15628 samples within 1 assay 
set.seed(42)
New_S_N <- NormalizeData(New_S_N, normalization.method = "LogNormalize", scale.factor = 5000)
New_S_N[["percent.mt"]] <- PercentageFeatureSet(New_S_N, pattern = "^mt-")
New_S_N_sub <- subset(New_S_N, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
all.genes_New_S_N <- rownames(New_S_N_sub)
New_S_N_sub <- ScaleData(New_S_N_sub, features = all.genes_New_S_N)
New_S_N_sub <- FindVariableFeatures(New_S_N_sub, selection.method = "vst", nfeatures = 11264) #50%
New_S_N_sub <- RunPCA(New_S_N_sub, features = VariableFeatures(object = New_S_N_sub))
New_S_N_sub <- FindNeighbors(New_S_N_sub, dims = 1:10)
New_S_N_sub <- FindClusters(New_S_N_sub, resolution = 0.5)


pdf(file='New_all_metrics.pdf', width=8, height=5,bg="white")
VlnPlot(New_S_N_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="orig.ident",
        ncol = 3,pt.size=0.001, cols = c("#45573d","#90a090")) + plot_annotation(
          title = 'New dataset',
          caption = 'Cell Ranger 7, includes intronic',
          theme = theme(plot.title = element_text(size = 16)))
dev.off()

unique(sapply(X = strsplit(colnames(New_S_N_sub), split = "_"), FUN = "[", 1))
New_S_N_sub$source<-(sapply(X = strsplit(colnames(New_S_N_sub), split = "_"), FUN = "[", 1))
table(New_S_N_sub$source)

set.seed(42)
New_S_N_sub <- RunUMAP(New_S_N_sub, reduction = "pca", dims = 1:20)
DimPlot(New_S_N_sub, reduction = "umap")
UMAPPlot(object = New_S_N_sub, group.by ="orig.ident",cols = c("Naive"="#f9a73e","Stroke"="#bf212f"))
UMAPPlot(object = New_S_N_sub, group.by ="orig.ident",split.by ="orig.ident",cols = c("Naive"="#f9a73e","Stroke"="#bf212f"))

cols<-c(sample(wes_palette("Darjeeling2", 17, type = c("continuous"))))
##find markers for clusters:
Idents(object = New_S_N_sub) <- New_S_N_sub@meta.data[["seurat_clusters"]]
New_S_N_sub.markers <- FindAllMarkers(New_S_N_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
New_S_N_sub.markers  %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.table(top10,"NEW_clusters_top10.txt")
pdf(file='NEW2_clusterHmap.pdf', width=20, height=10,bg="white")
DoHeatmap(New_S_N_sub, features = top10$gene,group.colors =cols,angle=0,size = 6,hjust=0.8,group.bar.height=0.02) + 
  scale_fill_viridis_b("Expression",option = "E") + theme(axis.text.y = element_text(size = 5)) + theme(axis.text.x = element_text(vjust=0.5, hjust=1)) 
dev.off()

New_S_N_sub<- SCTransform(New_S_N_sub, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:10)

DimPlot(New_S_N_sub, reduction = "umap",cols = cols, label=T)

###astro versus microG
astro_micro<-c('Slc1a3','Slc1a2','Aldh1l1','Aldoc','Gfap','Trem2','P2ry12','Aif1','Tmem119','Atp1b2')
gene_cols<-c(wes_palette("Royal1",10, type = c("continuous")))
VlnPlot(New_S_N_sub, astro_micro, stack = TRUE, sort = TRUE, flip = TRUE, same.y.lims=T, cols = gene_cols) +
  theme(legend.position = "none") + xlab("Clusters")

####remove microglia
New_S_N_sub_noM = subset(New_S_N_sub, Aif1 > 0.1| Trem2 >0.1, invert=T)
New_S_N_sub_noM

unique(sapply(X = strsplit(colnames(New_S_N_sub_noM), split = "_"), FUN = "[", 1))
New_S_N_sub_noM$source<-(sapply(X = strsplit(colnames(New_S_N_sub_noM), split = "_"), FUN = "[", 1))
table(New_S_N_sub_noM$source)

#40918 features across 2856 samples within 2 assays 
New_S_N_sub_noM<- SCTransform(New_S_N_sub_noM, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:10)
DimPlot(New_S_N_sub_noM, reduction = "umap")
UMAPPlot(object = New_S_N_sub_noM, group.by ="orig.ident",cols = c("Naive"="#f9a73e","Stroke"="#bf212f"))

DimPlot(New_S_N_sub_noM, reduction = "umap",cols = cols, label=T)
cells_cluster_origI<-as.data.frame(table(New_S_N_sub_noM@meta.data[["seurat_clusters"]], New_S_N_sub_noM@meta.data$orig.ident))
colnames(cells_cluster_origI)<-c('Cluster','Ident', 'CellN')
library(dplyr)
library(ggplot2)
library(stringr)
ggplot(cells_cluster_origI, aes(x = factor(Cluster), y = CellN, fill=Ident)) +
  geom_bar(stat="identity",width = 1,colour = "white") + theme_minimal() + 
  scale_fill_manual(values = cols) + xlab('Cluster') + ylab('Cell Count')

unique(sapply(X = strsplit(colnames(New_S_N_sub_noM), split = "_"), FUN = "[", 1))
New_S_N_sub_noM$source<-(sapply(X = strsplit(colnames(New_S_N_sub_noM), split = "_"), FUN = "[", 1))
table(New_S_N_sub_noM$source)

##overlay these clusters onto Visium slide
###now onto only old D10 slide
s87_D1 <- Load10X_Spatial(
  "./Visium/Faiz_Maryam__V10A06-087-D1/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "87_D1",
  slice = "D10_strokeD_",
  filter.matrix = TRUE,
  to.upper = FALSE,
  use.names = TRUE, unique.features = TRUE)
#normalize
s87_D1 <- SCTransform(s87_D1, assay = "87_D1", verbose = FALSE)
s87_D1 <- RunPCA(s87_D1, assay = "SCT", verbose = FALSE)
s87_D1 <- FindNeighbors(s87_D1, reduction = "pca", dims = 1:30)
s87_D1 <- FindClusters(s87_D1, verbose = FALSE)
s87_D1 <- RunUMAP(s87_D1, reduction = "pca", dims = 1:30)

anchors <- FindTransferAnchors(reference = New_S_N_sub_noM, query = s87_D1, normalization.method = "SCT")
#Found 173 anchors
predictions.assay <- TransferData(anchorset = anchors, refdata = New_S_N_sub_noM$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = s87_D1[["pca"]], dims = 1:10)

s87_D1[["predictions"]] <- predictions.assay
DefaultAssay(s87_D1) <- "predictions"
SpatialFeaturePlot(s87_D1, features = c('3','5','6','7','9','14','15'), alpha = c(0.5, 1), ncol=7)

###find markers for cluster 6
de_markers_cluster6 <- FindMarkers(New_S_N_sub_noM, ident.1 = 6)
write.table(de_markers_cluster6, 'Cluster6_strokeMarkers.txt')
