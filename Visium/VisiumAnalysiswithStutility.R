#install.packages("devtools")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.11")

#BiocManager::install("spdep",force = TRUE)
devtools::install_github("jbergenstrahle/STUtility")

library(STutility)
library(Matrix)
library(magrittr)
library(dplyr)
library(ggplot2)


#make infoTable#
#path to 87_B1
s87_B1 = "/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-087-B1"
sample_A <- paste0(s87_B1, "/filtered_feature_bc_matrix.h5")
spotfiles_A <- paste0(s87_B1, "/spatial/tissue_positions_list.csv")
imges_A <- paste0(s87_B1, "/spatial/tissue_hires_image.png")
json_A <- paste0(s87_B1, "/spatial/scalefactors_json.json")

#path to 87_D1
s87_D1 = "/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-087-D1"
sample_B <- paste0(s87_D1, "/filtered_feature_bc_matrix.h5")
spotfiles_B <- paste0(s87_D1, "/spatial/tissue_positions_list.csv")
imges_B <- paste0(s87_D1, "/spatial/tissue_hires_image.png")
json_B <- paste0(s87_D1, "/spatial/scalefactors_json.json")

#path to 88_B1
s88_B1 = "/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-088-B1"
sample_C <- paste0(s88_B1, "/filtered_feature_bc_matrix.h5")
spotfiles_C <- paste0(s88_B1, "/spatial/tissue_positions_list.csv")
imges_C <- paste0(s88_B1, "/spatial/tissue_hires_image.png")
json_C <- paste0(s88_B1, "/spatial/scalefactors_json.json")

#path to 88_C1
s88_C1 = "/Desktop/Faiz/Visium/Faiz_Maryam__V10A06-088-C1"
sample_D <- paste0(s88_C1, "/filtered_feature_bc_matrix.h5")
spotfiles_D <- paste0(s88_C1, "/spatial/tissue_positions_list.csv")
imges_D <- paste0(s88_C1, "/spatial/tissue_hires_image.png")
json_D <- paste0(s88_C1, "/spatial/scalefactors_json.json")


#paste all together
infoTable<-data.frame("samples" = c(sample_A,sample_B,sample_C,sample_D), "spotfiles" = c(spotfiles_A,spotfiles_B,spotfiles_C,spotfiles_D),
                      "imgs" = c(imges_A,imges_B,imges_C,imges_D), "json"=c(json_A,json_B,json_C,json_D), stringsAsFactors = FALSE)
rownames(infoTable)<-c("D2","D10","sham","D21")
ids <- rownames(infoTable)
infoTable<-data.frame("samples" = c(sample_A,sample_B,sample_C,sample_D), "spotfiles" = c(spotfiles_A,spotfiles_B,spotfiles_C,spotfiles_D),
                      "imgs" = c(imges_A,imges_B,imges_C,imges_D), "json"=c(json_A,json_B,json_C,json_D), ids, stringsAsFactors = FALSE)

se <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 100, 
                     min.gene.spots = 5,
                     min.spot.count = 100,
                     platform =  "Visium")

st.object <- GetStaffli(se)
#An object of class Staffli 
#8169 spots across 4 samples.
##normalize
se <- SCTransform(se)


#quality control/threshold plots#
##QC##
p1 <- ggplot() +
  geom_histogram(data = se[[]], aes(nFeature_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per spot")
p2 <- ggplot() +
  geom_histogram(data = se[[]], aes(nCount_RNA), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per spots")
gene_attr <- data.frame(nUMI = Matrix::rowSums(se@assays$RNA@counts), 
                        nSpots = Matrix::rowSums(se@assays$RNA@counts > 0))

p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")

p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total spots per gene")

(p1 - p2)/(p3 - p4)

##adjusting/transforming images for comparability##
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)
se <- MaskImages(object = se,thresholding = TRUE,iso.blur = 4)
ImagePlot(se, method = "raster")
#rotate images so they are all the same
transforms <- c(list("1" = list("angle" = 120)),list("2" = list("angle" = 184,"mirror.x" = T)),
                list("3" = list("angle" = 90,"mirror.x" = T,"shift.x"=-30)),list("4" = list("angle" = 90, "shift.x"=50)))
se.rotate <- WarpImages(se, transforms)
ImagePlot(se.rotate, method = "raster",ncols = 3, indices=c(1,2,4), type="processed",annotate =F)
se.rotate <- SCTransform(se.rotate)
##more QC overlayed onto se.rotate images##
mt.genes <- grep(pattern = "^mt-", x = rownames(se.rotate), value = TRUE)
se.rotate$percent.mito <- (Matrix::colSums(se.rotate@assays$RNA@counts[mt.genes, ])/Matrix::colSums(se.rotate@assays$RNA@counts))*100
# Collect all genes coding for ribosomal proteins
rp.genes <- grep(pattern = "^Rpl|^Rps", x = rownames(se.rotate), value = TRUE)
se.rotate$percent.ribo <- (Matrix::colSums(se.rotate@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se.rotate@assays$RNA@counts))*100

ST.FeaturePlot(se.rotate, features = "percent.mito", cols = c("#ece1ee", "#e0435e", "#9f80a7", "#43061e", "#0c0000"), pt.size = 1.3, ncol = 2,add.alpha = TRUE)
ST.FeaturePlot(se.rotate, features = "percent.ribo", cols = c("#ece1ee", "#e0435e", "#9f80a7", "#43061e", "#0c0000"), pt.size = 1, ncol = 2,add.alpha = TRUE)

umi_data <- GetAssayData(object = se.rotate, slot = "counts", assay = "RNA")
dim(umi_data)
#[1] 14718  8169

# Calculate gene attributes
gene_attr <- data.frame(mean = rowMeans(umi_data),
                        detection_rate = rowMeans(umi_data > 0),
                        var = apply(umi_data, 1, var), 
                        row.names = rownames(umi_data)) %>%
  mutate(log_mean = log10(mean), log_var = log10(var))

# Obtain spot attributes from Seurat meta.data slot
spot_attr <- se[[c("nFeature_RNA", "nCount_RNA")]]

p1 <- ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha = 0.3, shape = 16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  ggtitle("Mean-variance relationship")

# add the expected detection rate under Poisson model
x = seq(from = -2, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
p2 <- ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha = 0.3, shape = 16) + 
  geom_line(data = poisson_model, color='red') +
  ggtitle("Mean-detection-rate relationship")

p1 - p2
##start to factorize/cluster##
set.seed(8)
se.rotate <- RunNMF(se.rotate, nfactors = 20)
se.rotate <- FindNeighbors(object = se.rotate, verbose = FALSE, reduction = "NMF", dims = 1:20)
se.rotate <- FindClusters(object = se.rotate, verbose = FALSE)
se.rotate <- RunUMAP(se.rotate, reduction = "NMF", dims = 1:20, n.neighbors = 10)

rm(se)

#look at isolated gene expression in any sample#
#These are just colour palettes I have generated that I like:):
library(scales)
pals <- c("#E7ECF0","#F5C245","#EB6359","#696267","#529085")
mar_pal<- c("#E0DDD6","#C1C5B6","#CA212C","#A01422","#49323E")
beach <-c("#DCE6E5","#B19C95","#1B6765","#345F54","#033238")
sunset <-c("#D7AE94","#B0857A","#FE8B48","#7D9BB3","#172426")
forest <- c("#B8CACB","#4B7C45","#595E4A","#3D4A43","#3A2838")
Mtn_sunset<-c("#e3cabb","#d8a19b","#ffb567","#bf433e","#32272d","#0b181e")
show_col(beach)
#and a custom theme that removes legends and makes the resulting figure more clean
customtheme <- theme(legend.position = "None", plot.title = element_blank())
#make the figure with:
FeatureOverlay(se.rotate, features = "Cd81", 
               sampleids = c(1,2,4),#here is where you select sample to look at
               add.alpha = TRUE, #sets it so only places with gene expression has dots
               pt.size = 0.8,
               pt.alpha = 0.6,
               cols = Mtn_sunset,#you can switch out the palette here
               ncols = 3,#this just dictates how many columns for the resulting figure
               sample.label = T,
               show.sb = FALSE) #& customtheme

##can also plot a list of genes over visium images
blood <- c("Hbb-bt","Hbb-bs","Tgfbi")# this can be as many genes as you want
se.rotate$blood <- (Matrix::colSums(se.rotate@assays$SCT@counts[blood, ]))
ST.FeaturePlot(se.rotate, features = "blood", cols = c("#ece1ee", "#e0435e", "#9f80a7", "#43061e", "#0c0000"), pt.size = 1.3, ncol = 2)
astro<-c("Gfap","Serpina3n","S100b","C4b","Aqp4","C3","Il33","Id3","Stat1")
se.rotate$astro <- (Matrix::colSums(se.rotate@assays$SCT@counts[astro, ]))
ST.FeaturePlot(se.rotate, features = "astro",
               cols = mar_pal,
               pt.size = 0.8, ncol = 4)
micro<-c("Il1b", "Hexb","Sall1","P2ry12","Trem2")
se.rotate$micro <- (Matrix::colSums(se.rotate@assays$SCT@counts[micro, ]))
ST.FeaturePlot(se.rotate, features = "micro",
               cols = beach,
               pt.size = 0.8, ncol = 2)


##identify spatially variable genes
library(spdep)
spatgenes <- CorSpatialGenes(se.rotate)

head(VariableFeatures(se.rotate))
V_feats<- VariableFeatures(se.rotate)
setwd('/Users/erica/Desktop/Faiz/Visium')
write.table(V_feats, "Visium_Vfeats.txt")

library(devtools)
##start to factorize/cluster with NMF##
se.rotate <- RunNMF(se.rotate, nfactors = 20)
#########
cscale <- c("lightgray", "mistyrose", "red", "darkred", "black")
bscale <- c("#f3a001","#ccccf4","#4c3277","#0077c5","#313853")
ST.DimPlot(se.rotate, 
           dims = c(1),#sets the factor to look at
           ncol = 4, # Sets the number of columns at dimensions level
           grid.ncol = 1, # Sets the number of columns at sample level
           reduction = "NMF", 
           pt.size = 1, 
           center.zero = F, 
           cols = bscale, 
           show.sb = FALSE)

ST.DimPlot(se.rotate, 
           dims = c(7),#sets the factor to look at
           ncol = 4, # Sets the number of columns at dimensions level
           grid.ncol = 1, # Sets the number of columns at sample level
           reduction = "NMF", 
           pt.size = 1, 
           center.zero = F, 
           cols = beach, 
           show.sb = FALSE)

FactorGeneLoadingPlot(se.rotate, factor =5, topn = 40, dark.theme = FALSE)
factor_geneList_5.10.19<- SummarizeAssocFeatures(
  se.rotate,
  dims = c(5,10,19),#need to format this as a list of multiple dims, here it is factor 5,10&9
  features.return = 50,
  features.use = NULL)

write.table(factor_geneList_5.10.19[[1]], "factors_5_10_19_genes.txt")

##further digging into cluster over time and DEG between cluster and border of cluster
#se.rotate <- RunUMAP(se.rotate, reduction = "NMF", dims = 1:20, n.neighbors = 10)
print(se.rotate[["umap"]])
se.rotate <- FindNeighbors(object = se.rotate, verbose = FALSE, reduction = "NMF", dims = 1:20)
se.rotate <- FindClusters(object = se.rotate, verbose = FALSE)
library(RColorBrewer)
n <- 19
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
F1<-c("Plp1", "Mbp", "Mobp", "Cnp", "Mal", "Fth1", "Trf", "Mag", "Cldn11", "Apod", 
      "Mog", "Cryab", "Qdpr", "Ptgds", "Tspan2", "Ppp1r14a", "Car2", "Plekhb1", "Sept4", "Dbndd2")

p1<- FeaturePlot(se.rotate, features = F1, reduction="NMF")
p2<-ST.FeaturePlot(object = se.rotate, features = "seurat_clusters", cols = col_vector, pt.size = 1, ncol = 2)
p1-p2
#to just look at clusters over visium data
FeatureOverlay(se.rotate, features = "seurat_clusters",reduction = "umap", sampleids = c(1,2,4), ncols = 3,pt.size = 1)

##isolating individual clusters##
se.rotate <- SetIdent(se.rotate, value = "seurat_clusters") #or SCT_snn_res.0.8
se.rotate <- RegionNeighbours(se.rotate, id = "12",keep.within.id = T, verbose = TRUE)
customtheme <- theme(legend.position = "None", plot.title = element_blank())
FeatureOverlay(se.rotate, features = "nbs_12", ncols = 4, sampleids = c(3,1,2,4), cols = c("#a1a1a1", "#675364"), pt.size = 1) & customtheme
se.rotate <- RegionNeighbours(se.rotate, id = "16",keep.within.id = T, verbose = TRUE)
FeatureOverlay(se.rotate, features = "nbs_16", ncols = 4, sampleids = c(3,1,2,4), cols = c("#070c3a", "#ffb916"), pt.size = 1)  & customtheme


###run DE test on these gene clusters^
library(magrittr)
library(dplyr)
se.rotate <- SetIdent(se.rotate, value = "nbs_12")
nbs_12.markers <- FindMarkers(se.rotate, ident.1 = "12", ident.2 = "nbs_12")
nbs_12.markers$gene <- rownames(nbs_12.markers)
se.subset <- SubsetSTData(se.rotate, expression = nbs_12 %in% c("12", "nbs_12"))
sorted.marks <- nbs_12.markers %>% top_n(n = 50, wt = abs(avg_log2FC))
sorted.marks <- sorted.marks[order(sorted.marks$avg_log2FC, decreasing = T), ]
customtheme <- theme(legend.position = "None", plot.title = element_blank())
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("darkturquoise", "darkorchid1"), disp.min = -2, disp.max = 2)
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("#070c3a", "#ffb916"))+ scale_fill_viridis(option='inferno')


write.table(sorted.marks, "nbs12_genes.txt")
#and now for cluster 4
se.rotate <- SetIdent(se.rotate, value = "nbs_16")
nbs_16.markers <- FindMarkers(se.rotate, ident.1 = "16", ident.2 = "nbs_16")
nbs_16.markers$gene <- rownames(nbs_16.markers)
se.subset <- SubsetSTData(se.rotate, expression = nbs_16 %in% c("16", "nbs_16"))
sorted.marks <- nbs_16.markers %>% top_n(n = 50, wt = abs(avg_log2FC))
sorted.marks <- sorted.marks[order(sorted.marks$avg_log2FC, decreasing = T), ]
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("forestgreen", "olivedrab3"), disp.min = -2, disp.max = 2)
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("#a1a1a1", "#675364"))+ scale_fill_viridis(option='inferno')


se.rotate <- SetIdent(se.rotate, value = "nbs_12")
nbs_12.markers <- FindMarkers(se.rotate, ident.1 = "12", ident.2 = "nbs_12")
nbs_12.markers$gene <- rownames(nbs_12.markers)
se.subset <- SubsetSTData(se.rotate, expression = nbs_12 %in% c("12", "nbs_12"))
sorted.marks <- nbs_12.markers %>% top_n(n = 50, wt = abs(avg_log2FC))
sorted.marks <- sorted.marks[order(sorted.marks$avg_log2FC, decreasing = T), ]
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("#a1a1a1", "#675364"))+ scale_fill_viridis(option='inferno')

