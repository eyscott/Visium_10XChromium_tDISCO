
NeuN <-c('2_2_2','2_4_8','3_2_8','3_3_2')
Gfap<-c('1_1_1','1_1_2','1_1_7',	'1_1_8')
Use<-c(NeuN,Gfap)
D10_Norm_geneMat_GvN<-D10_Norm_geneMat[ ,Use]
#convert zero to NA
D10_Norm_geneMat_GvN_NA<-na_if(D10_Norm_geneMat_GvN, 0)
D10_Norm_geneMat_GvN_lessNA<- D10_Norm_geneMat_GvN_NA[which(rowMeans(!is.na(D10_Norm_geneMat_GvN_NA)) > 0.2), ] 
D10_Norm_geneMat_GvN_lessNA[is.na(D10_Norm_geneMat_GvN_lessNA)] <- 0
write.table(D10_Norm_geneMat_GvN_lessNA, "D10_norm_GvN.txt")

D10_Norm_geneMat_GvN_lessNA<-as.data.frame(D10_Norm_geneMat_GvN_lessNA)
D10_Norm_geneMat_GvN_lessNA<- tibble::rownames_to_column(D10_Norm_geneMat_GvN_lessNA, "Gene_name")
D10_Norm_geneMat_GvN_melt<-melt(D10_Norm_geneMat_GvN_lessNA,id.vars = 'Gene_name')

D10_Norm_geneMat_GvN_melt_sum<- ddply(D10_Norm_geneMat_GvN_melt, c("Gene_name"), summarise,
                                   sum = sum(value),sd=sd(value), sem = sd(value)/sqrt(length(value)))

D10_exp_NORM_GvN_summarized <- merge(D10_Norm_geneMat_GvN_lessNA,D10_Norm_geneMat_GvN_melt_sum, by="Gene_name")
#Remove Gm genes
D10_exp_NORM_GvN_summarized_noGm <- D10_exp_NORM_GvN_summarized[!grepl("^Gm", D10_exp_NORM_GvN_summarized$Gene_name),]#now 35747

sub <- subset(D10_exp_NORM_GvN_summarized_noGm, c(sum > 25 & sem > 25)) 
rownames(sub)<-sub$Gene_name
data_mat<- as.matrix(sub[,-c(1,10:12)])

hr <- hclust(as.dist(1-cor(t(data_mat), method="pearson")),
             method="average")

lmat = rbind(c(0,4),c(0,3),c(2,1))
lwid = c(0.5,4)
lhei = c(1,0.5,4)

library(gplots)
library(viridis)
options("scipen"=100, "digits"=4)
col_breaks <- c(seq(0,50,1),seq(51,150,10),seq(151,500,50),seq(501,1000,100))
pdf(file='D10_GvN.pdf', width=5, height=10,bg="white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
heatmap.2(data_mat,   
          trace="none", 
          labRow = hr$labels,
          margins =c(5,5),    
          col=viridis,      
          breaks=col_breaks,   
          labCol = c("NeuN+","NeuN+","NeuN+","NeuN+","Gfap+","Gfap+", "Gfap+","Gfap+"),
          dendrogram="none", 
          Colv=F,
          Rowv=as.dendrogram(hr),
          cexCol = 1,
          cexRow = 0.5,
          adjCol = c(0.9,0.6),
          adjRow = c(0.35,0.2),
          hclustfun = hclust,
          keysize = 5,
          key.xlab = "TPM",
          key.title = NA,
          key.par = list(mar=c(1,3,1,6)),
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)            
dev.off() 

##Annotated above genes and selected astrocytic and neuronal markers
astro<-c('Slc1a3','Slc1a2','Gfap','S100b','Aldh1l1','Col9a3','Aldoa','Tpi1','Aqp4')
Neuro<-c('Blank','Dlg4','Dcx','Sypl2','Psd2','Grin3b','Slc17a8','Slc17a1','Cluap1','Tut1','Mbd1','Cabs1','Spem1','Mir5110','Hes7','Tango6','Ltbr','Htr2a')
all<-c(astro,Neuro)

D10_GvN_submat<- D10_Norm_geneMat_GvN_lessNA[match(all,D10_Norm_geneMat_GvN_lessNA$Gene_name),]
D10_GvN_submat_mat<-as.matrix(D10_GvN_submat[ ,-1])

row_labels<-data.frame(cbind(D10_GvN_submat$Gene_name,D10_GvN_submat))
row_labels<-D10_GvN_submat$Gene_name

col_breaks <- c(seq(0,50,1),seq(51,150,10),seq(151,500,50),seq(501,1000,100))
pdf(file='D10_GvN_2.pdf', width=5, height=8,bg="white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
heatmap.2(as.matrix(D10_GvN_submat),    # data matrix
          trace="none", 
          labRow = row_labels,
          margins =c(5,5),  
          col=viridis,    
          breaks=col_breaks, 
          labCol = c("NeuN","NeuN","NeuN","NeuN","Gfap","Gfap","Gfap","Gfap"),
          dendrogram="none",  
          Colv=F,
          Rowv=F,
          cexCol = 1.1,
          cexRow = 0.8,
          adjCol = c(0.9,0.6),
          adjRow = c(0.2,0.2),
          hclustfun = hclust,
          keysize = 5,
          key.xlab = "TPM",
          key.title = NA,
          key.par = list(mar=c(1,3,1,6)), 
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)               
dev.off() 
