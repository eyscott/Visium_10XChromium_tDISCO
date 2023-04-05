A<-c("2_1_2",'2_1_4','2_1_6','2_4_4')
B<-c('1_1_1','1_1_2','1_1_7','1_1_8','2_4_6','3_2_1','3_2_5','3_3_5')
C<-c('2_4_7','2_2_4','2_2_6','2_2_8')
D<-c('3_1_2','3_1_4','3_2_2','3_2_4','3_2_6','3_3_8')

Fade<-c(A,B,C,D)

D10_Norm_geneMat_Fade<-D10_Norm_geneMat[ ,Fade]

lncRNA<-c('Break','Gm44432','Gm44039','Gm36787','Gm29946','Gm31763','Gm49953','Gm49076','Gm32172','4930519H02Rik','1110002E22Rik')
iv<-c('Break1','Chst3',"Gan",'Tprkb','Acap3','Gnai1','Klf13')
i<-c('Break2','Hpcal4','Basp1','Frrs1l','Penk','Nudcd3','Furin','Taok3','Dnajb2','Dennd5a','Egln2','Cfl1','B4galt6',
          'Ntrk2','Ndst1','Mapk1','Mapk10','Ccni','Timm17a','Bscl2','Spred1','Ubr4','Ywhaq','Rap1b','Basp1')
ii<-c('Break3','Fam102a','Cant1','Exoc3','Ice1','Rcor1','Mphosph9','Siah1a','Senp8','Gigyf2','Tmem119',
              'Zfp664','Zfp26','Olig2','Asxl1','Vcp','Rubcn','Toe1')
core<- c('Break8','Fkbp11','Sdhb','Gltp','Ltc4s','Dffb','Gm20661','Zfp799')
compilation<-c(lncRNA,peripheral,Not_injury,Injury)
compilation_2<-c(peripheral,Injury,core)
D10_Norm_geneMat_Fade_compilation<-D10_Norm_geneMat_Fade[rownames(D10_Norm_geneMat_Fade) %in% compilation, ] 
D10_Norm_geneMat_Fade_compilation2<-D10_Norm_geneMat_Fade[rownames(D10_Norm_geneMat_Fade) %in% compilation_2, ] 

D10_Norm_geneMat_Fade_compilation_ordered<- D10_Norm_geneMat_Fade_compilation[match(compilation,rownames(D10_Norm_geneMat_Fade_compilation)),]
D10_Norm_geneMat_Fade_compilation2_ordered<- D10_Norm_geneMat_Fade_compilation2[match(compilation_2,rownames(D10_Norm_geneMat_Fade_compilation2)),]

col_breaks <- c(seq(0,10,0.5),seq(11,50,1),seq(51,100,5),seq(101,500,5))
pdf(file='D10_Fade_comp.pdf', width=6, height=8,bg="white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
heatmap.2(D10_Norm_geneMat_Fade_compilation_ordered,  
          trace="none", 
          labRow = rownames(D10_Norm_geneMat_Fade_compilation_ordered),
          margins =c(5,5),   
          col=viridis,     
          breaks=col_breaks,   
          labCol =c(seq(1,22,1)),
          dendrogram="none",   
          Colv=F,
          Rowv=F,
          srtCol =90,
          cexCol = 2,
          cexRow = 0.8,
          adjCol = c(0.8,0.2),  
          adjRow = c(0.20,0.2),
          keysize = 5,
          key.xlab = "TPM",
          key.title = NA,
          key.par = list(mar=c(1,3,1,6)), 
          densadj = 0.4,
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)             
dev.off()
