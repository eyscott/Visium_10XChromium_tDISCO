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

UMAPPlot(object = New_S_N_sub_noM, group.by ="orig.ident",cols = c("Naive"="#f9a73e","Stroke"="#bf212f"))

cells_cluster_origI<-as.data.frame(table(New_S_N_sub_noM@meta.data[["seurat_clusters"]], New_S_N_sub_noM@meta.data$orig.ident))
colnames(cells_cluster_origI)<-c('Cluster','Ident', 'CellN')
library(dplyr)
library(ggplot2)
library(stringr)
ggplot(cells_cluster_origI, aes(x = factor(Cluster), y = CellN, fill=Ident)) +
  geom_bar(stat="identity",width = 1,colour = "white") + theme_minimal() + 
  scale_fill_manual(values = cols) + xlab('Cluster') + ylab('Cell Count')

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
SpatialFeaturePlot(s87_D1, features = c('3','5','6','7','9','10'), alpha = c(0.5, 1), ncol=6)

###find markers for cluster 6
de_markers_cluster6 <- FindMarkers(New_S_N_sub_noM, ident.1 = 6)
write.table(de_markers_cluster6, 'Cluster6_strokeMarkers.txt')

##########
###processing separate ACSA-2+hi 10X Chromium astrocytes
#########
UI.data <- Read10X(data.dir = "./Uninjured")
UI <- CreateSeuratObject(counts = UI.data, project = "Naive", min.cells = 3, min.features = 200)
S.data <- Read10X(data.dir = "./Stroke")
S <- CreateSeuratObject(counts = S.data, project = "Stroke", min.cells = 3, min.features = 200)
#merge two old
UIandS <- merge(UI, y = S, add.cell.ids = c("Naive", "Stroke"), project = "UIandS")
UIandS
set.seed(42)
UIandS <- NormalizeData(UIandS, normalization.method = "LogNormalize", scale.factor = 5000)
UIandS[["percent.mt"]] <- PercentageFeatureSet(UIandS, pattern = "^mt-")
UIandS_sub <- subset(UIandS, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
#15472 features across 5485 samples within 1 assay
all.genes_UIandS <- rownames(UIandS_sub)
UIandS_sub <- ScaleData(UIandS_sub, features = all.genes_UIandS)
UIandS_sub<- SCTransform(UIandS_sub, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:10)
UIandS_sub <- FindNeighbors(UIandS_sub, dims = 1:10)
UIandS_sub <- FindClusters(UIandS_sub, resolution = 0.5)
UMAPPlot(object = UIandS_sub, group.by ="orig.ident",cols = c("Naive"="#f9a73e","Stroke"="#bf212f"))
DimPlot(UIandS_sub, reduction = "umap",cols = c('0'='#3b5a9d','1' = '#e64a3d', '2' = '#f2cf59','3'='#fb8e7e',
                                                '4'='#c5d7c0','5'='#253656','6'='#8ec9bb'), label=T)

astro_micro<-c('Slc1a3','Slc1a2','Aldh1l1','Aldoc','Gfap','Trem2','P2ry12','Aif1','Tmem119','Ogn','Lum','Ptgds')
gene_cols<-c(wes_palette("Royal1",11, type = c("continuous")))
VlnPlot(UIandS_sub, astro_micro, stack = TRUE, sort = F, flip = TRUE, same.y.lims=T, cols = gene_cols) +
  theme(legend.position = "none") + xlab("Clusters") +
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(legend.title = element_text(colour="black", size=20, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 18)) +
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(axis.text = element_text(colour="black", size = 18))

UIandS_sub.markers <- FindAllMarkers(UIandS_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UIandS_sub.markers  %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
oldCcolours<-c('0'='#3b5a9d','1' = '#e64a3d', '2' = '#f2cf59','3'='#fb8e7e',
  '4'='#c5d7c0','5'='#253656','6'='#8ec9bb')
pdf(file='ACSA2hi_10X_top10_clusterHmap.pdf', width=20, height=10,bg="white")
DoHeatmap(UIandS_sub, features = top10$gene,group.colors =oldCcolours,angle=45,size = 6,hjust=0.8,group.bar.height=0.02) + 
  scale_fill_viridis_b("Expression",option = "E") + theme(axis.text.y = element_text(size = 5)) + theme(axis.text.x = element_text(vjust=0.5, hjust=1)) 
dev.off()


##check cell proportions per cluster
cols<-c(sample(wes_palette("Darjeeling2", 17, type = c("continuous"))))
cells_cluster_origI<-as.data.frame(table(UIandS_sub@meta.data[["seurat_clusters"]], UIandS_sub@meta.data$orig.ident))
colnames(cells_cluster_origI)<-c('Cluster','Ident', 'CellN')
ggplot(cells_cluster_origI, aes(x = factor(Cluster), y = CellN, fill=Ident)) +
  geom_bar(stat="identity",width = 1,colour = "white") + theme_minimal() + 
  scale_fill_manual("Source",labels = c("Uninjured", "Injured"),values = cols) + xlab('Cluster') + ylab('Cell Count')+
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(legend.title = element_text(colour="black", size=20, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 18)) +
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(axis.text = element_text(colour="black", size = 18))

unique(sapply(X = strsplit(colnames(UIandS_sub), split = "_"), FUN = "[", 1))
UIandS_sub$source<-(sapply(X = strsplit(colnames(UIandS_sub), split = "_"), FUN = "[", 1))
table(UIandS_sub$source)

##Make volcanoe plot
UIandS_sub_3v0$DEG <- "None"
UIandS_sub_3v0$DEG[UIandS_sub_3v0$avg_log2FC > 0.48 & UIandS_sub_3v0$p_val_adj < 0.05] <- "Proximal"
UIandS_sub_3v0$DEG[UIandS_sub_3v0$avg_log2FC < -0.5 & UIandS_sub_3v0$p_val_adj < 0.05] <- "Distal"
UIandS_sub_3v0 <- tibble::rownames_to_column(UIandS_sub_3v0, "Gene")

UIandS_sub_3v0$Label <- NA
UIandS_sub_3v0$Label[UIandS_sub_3v0$DEG != "None"] <- UIandS_sub_3v0$Gene[UIandS_sub_3v0$DEG != "None"]
mycolors <- c("#384756", "#9dd4cf")
names(mycolors) <- c("Proximal", "Distal")
library(ggrepel)
ggplot(UIandS_sub_3v0, aes(x=avg_log2FC, y=-log10(p_val_adj), col=DEG, label=Label)) +
  geom_point(size = 3) + 
  theme_minimal() +
  geom_text_repel() +
  scale_colour_manual(name="Upregulated in:",values = mycolors) + 
  labs(x="Log2FC", y="-log(adj. p-value)") +
  geom_vline(xintercept=c(-0.6, 0.6), col="#b8300f") +
  geom_hline(yintercept=-log10(0.05), col="#b8300f") +
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(legend.title = element_text(colour="black", size=20, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 18)) +
  theme(axis.title = element_text(colour="black",size = 20))+
  theme(axis.text = element_text(colour="black", size = 10))
UIandS_sub_3v0 <- FindMarkers(UIandS_sub, ident.1 = 3, ident.2 = 0)

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
anchors <- FindTransferAnchors(reference = UIandS_sub, query = s87_D1, normalization.method = "SCT")
#Found 147 anchors

predictions.assay <- TransferData(anchorset = anchors, refdata = UIandS_sub$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = s87_D1[["pca"]], dims = 1:10)

s87_D1[["predictions"]] <- predictions.assay
DefaultAssay(s87_D1) <- "predictions"
SpatialFeaturePlot(s87_D1, features = c('0','1','2','3','4','5','6'), alpha = c(0.5, 1), ncol=7)

