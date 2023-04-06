library(reshape2)
library(plyr)
library(dplyr)
allData <- read.table('D10_LFQ_untarg_04_09.txt', header=T, stringsAsFactors = F)
#convert zero to NA
allData_NA<-na_if(allData, 0)
allData_lessNA<- allData_NA[which(rowMeans(!is.na(allData_NA)) > 0.1), ] ##brings to 316
##replace NAs with zeros
allData_lessNA[is.na(allData_lessNA)] <- 0
##match with gene names
names<-read.table("uniprot4_725_LFQ_mod.txt", header=T,stringsAsFactors = F)
allData_lessNA_names <-merge(allData_lessNA, names, by.x='Protein_ID',by.y='From')
rownames(allData_lessNA_names)<- make.names(allData_lessNA_names$To, unique=TRUE)
#removeProtein ID column
allData_Simpl<-allData_lessNA_names[ ,-c(1,22:24)]
library(tibble)
allData_se_prep <- tibble::rownames_to_column(allData_Simpl, "name")
allData_unique <- make_unique(allData_lessNA_names, "To", "Protein_ID", delim = "/t")

data_d10_design<-read.table('D10_design.txt',header=T, stringsAsFactors = F)

allData_se <- make_se(allData_unique, c(2:21), data_d10_design)
plot_frequency(allData_se)
allData_se_filt <- filter_missval(allData_se, thr = 5)
library(wesanderson)
plot_frequency(allData_se_filt)
plot_numbers(allData_se_filt)
plot_coverage(allData_se_filt)
plot_missval(allData_se_filt)
plot_detect(allData_se_filt)

# Get a logical vector
MNAR <- names(allData_se_filt)

# Perform a mixed imputation
mixed_imputation <- impute(
  allData_se_filt, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "zero") # imputation function for MNAR
# The wrapper function performs the full analysis

DE_analysis <- function(se) {
  se %>% 
    test_diff(., type = "manual", test = c("D_vs_B")) %>%
    add_rejections(., alpha = 0.05, lfc = 1) %>% 
    get_results()
}
allData_se_filt_results <- DE_analysis(allData_se_filt)
# Test manually defined comparisons
allData_se_manual <- test_diff(allData_se, type = "manual", 
                              test = c("D_vs_B"))

plot_pca(allData_se_manual, x = 1, y = 2, n = 317, point_size = 4)
