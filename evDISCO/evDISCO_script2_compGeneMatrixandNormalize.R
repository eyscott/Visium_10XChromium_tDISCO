##compile all D10 data
library(reshape2)
library(plyr)
library(dplyr)

#round 1
#Get gene length
gene_L_1 <- read.table("FC_D10_annot.txt", header = T, stringsAsFactors = F)

####Get expression data
D10_exp_1 <- read.table("FC_D10.txt", header = T, stringsAsFactors = F)

#readcounts
sum_counts_1<-data.frame(Sums=colSums(D10_exp_1[,-1]))
rownames(sum_counts_1) <- paste("1",rownames(sum_counts_1),sep="_")
####round 2.5####
gene_L_2 <- read.table("FC_D10_annot.txt", header = T, stringsAsFactors = F)
####Get expression data
D10_exp_2 <- read.table("FC_D10.txt", header = T, stringsAsFactors = F)
#readcounts
sum_counts_2<-data.frame(Sums=colSums(D10_exp_2[,-1]))
rownames(sum_counts_2) <- paste("2",rownames(sum_counts_2),sep="_")

####round 3####
#Get gene length
gene_L_3 <- read.table("FC_D10_annot.txt", header = T, stringsAsFactors = F)
D10_exp_3 <- read.table("FC_D10.txt", header = T, stringsAsFactors = F)
#readcounts
sum_counts_3<-data.frame(Sums=colSums(D10_exp_3[,-1]))
rownames(sum_counts_3) <- paste("3",rownames(sum_counts_3),sep="_")

##Compile everything##
sum_counts_all<-rbind(sum_counts_1,sum_counts_2,sum_counts_3)
sum_counts_all <- rbind(data.frame(id="Round1",sum_counts_1),
                          data.frame(id="Round2",sum_counts_2),
                          data.frame(id="Round3",sum_counts_3))
library(ggplot2)
library(scales)
library(wesanderson)

# Basic violin plot
pdf(file='Sif1A_readcoutns.pdf', width=6, height=4,bg="white")
ggplot(sum_counts_all, aes(x=id,y=Sums,fill=id)) + 
  geom_violin() + ylab("Read counts") + xlab(" ") +
  scale_fill_manual(values=(wes_palette(n=3, name="Royal1"))) +
  theme(legend.position="none") +geom_dotplot(binaxis='y', stackdir='center',
                                              position=position_dodge(1))
dev.off()

##now compile expression data##
rownames(D10_exp_1)->D10_exp_1$Gene_id
rownames(D10_exp_2)->D10_exp_2$Gene_id
rownames(D10_exp_3)->D10_exp_3$Gene_id

D10_exp_1_melt <- melt(D10_exp_1,id.vars= c("Gene_id"))
D10_exp_2_melt <- melt(D10_exp_2,id.vars= c("Gene_id"))
D10_exp_3_melt <- melt(D10_exp_3,id.vars= c("Gene_id"))

D10_exp_melt_all <- rbind(data.frame(id="Round1",D10_exp_1_melt),
                        data.frame(id="Round2",D10_exp_2_melt),
                        data.frame(id="Round3",D10_exp_3_melt))
D10_exp_RAW<-dcast(D10_exp_melt_all, Gene_id~variable+id,value.var='value',
                   fun.aggregate = mean, sep = "\t")
write.table(D10_exp_RAW, "D10_RAW_mat.txt")

##normalize for gene length
gene_L_all<-rbind(gene_L_1,gene_L_2,gene_L_3)
D10_exp_melt_all_L <- merge(D10_exp_melt_all,gene_L_all, by.x="Gene_id", by.y="GeneID")
D10_exp_melt_all_L$Length_kb <- (D10_exp_melt_all_L$Length)/1000
D10_exp_melt_all_L$value_L <- (D10_exp_melt_all_L$value/D10_exp_melt_all_L$Length_kb)
#get rid of data we don't need anymore
D10_exp_melt_all_L_red <- D10_exp_melt_all_L[ ,c("Gene_id","id", "variable","value_L")]

#normalize for library depth
D10_exp_melt_all_L_red_presum<- ddply(D10_exp_melt_all_L_red, c("variable"), summarise,
                              Count_sum = sum(value_L))
D10_exp_melt_all_L_red_presum$RPM <- ((D10_exp_melt_all_L_red_presum$Count_sum)/1000000)
D10_exp_melt_all_L_RPM <- merge(D10_exp_melt_all_L_red, D10_exp_melt_all_L_red_presum, by="variable")
D10_exp_melt_all_L_RPM$value_L_RPM <- (D10_exp_melt_all_L_RPM$value_L/D10_exp_melt_all_L_RPM$RPM)
##get rid of some of the temp columns before making gene matrix
D10_exp_melt_all_L_RPM_red <- D10_exp_melt_all_L_RPM[ ,c("Gene_id","id", "variable", "value_L_RPM")]

D10_exp_NORM <- dcast(D10_exp_melt_all_L_RPM_red, Gene_id~variable+id,value.var='value_L_RPM',
                      fun.aggregate = mean, sep = "\t")

#change Gene_ids to gene name with Biomart
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "mmusculus_gene_ensembl")

annot <- getBM(attributes = c('ensembl_gene_id','external_gene_name','gene_biotype'),
               values = row.names(D10_exp_NORM), 
               mart = ensembl)
D10_exp_NORM_annot <- merge(D10_exp_NORM,annot, by.x='Gene_id', by.y="ensembl_gene_id")
D10_exp_NORM_annot_unique <- distinct(D10_exp_NORM_annot, external_gene_name, .keep_all= TRUE)
D10_exp_NORM_annot_unique_nopseudo <- D10_exp_NORM_annot_unique[!grepl("pseudogene|misc", D10_exp_NORM_annot_unique$gene_biotype),]
rownames(D10_exp_NORM_annot_unique_nopseudo)<-D10_exp_NORM_annot_unique_nopseudo$external_gene_name
D10_Norm_geneMat<-as.matrix(D10_exp_NORM_annot_unique_nopseudo[ ,c(2:34)])
write.table(D10_Norm_geneMat, "D10_Norm_mat.txt")

#generating statistic metrics
genes_per_cell <- data.frame(Matrix::colSums(D10_Norm_geneMat>0))       
genes_per_cell_split <- t(data.frame(strsplit(as.character(rownames(genes_per_cell)),"_",fixed = T)))
genes_per_cell_Round <- cbind(genes_per_cell,genes_per_cell_split[ ,4])
colnames(genes_per_cell_Round)<-c('GeneCount', 'Round')

pdf(file='Sif1B_GeneCounts.pdf', width=6, height=4,bg="white")
ggplot(genes_per_cell_Round, aes(x=Round,y=GeneCount,group=Round)) + 
  geom_violin(aes(fill=Round)) + scale_fill_manual(values=(wes_palette(n=3, name="Royal1"))) +
  ylab("Gene counts") + xlab(" ") +
  theme(legend.position="none") + geom_dotplot(binaxis='y', stackdir='center',
                                               position=position_dodge(1),dotsize = 0.6, color='white')
dev.off()

#look at gene content
D10_exp_NORM_prep <- D10_exp_NORM_annot_unique[ ,c(2:36)]
D10_exp_NORM_prepP <- D10_exp_NORM_prep[!grepl("^IG|*_gene$|TEC|*_pseudogene$|ribozyme|pseudogene", D10_exp_NORM_prep$gene_biotype),]#now 35747
D10_exp_NORM_prepP_melt<-melt(D10_exp_NORM_prepP, id.vars = c('external_gene_name','gene_biotype'))
D10_exp_NORM_prepP_melt_split <- t(data.frame(strsplit(as.character(D10_exp_NORM_prepP_melt$variable),"_",fixed = T)))
D10_exp_NORM_prepP_melt_Round <- cbind(D10_exp_NORM_prepP_melt,D10_exp_NORM_prepP_melt_split[ ,4])
colnames(D10_exp_NORM_prepP_melt_Round)<-c("Gene_name","gene_biotype","variable","value","Round")
library(RColorBrewer)
n <-20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(D10_exp_NORM_prepP_melt_Round,aes(x=factor(1),weight=value,fill=gene_biotype)) + 
  geom_bar(width=1,alpha = 1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="RNA type")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  theme(axis.text = element_text(colour="black", size = 8)) +
  scale_fill_manual(values = col_vector) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + facet_grid(cols = vars(Round))
