setwd('F:/tcga')
getwd()

install.packages("F:/tcga/MNMO_0.1.1.tar.gz",repos = NULL, type="source")
library(MNMO)
library(psych)
data(network)
data(protein_info)
data(mir_gene_network)
###############BRCA DATA
data(brca_normal_s)
data(brca_cancer_s)
data(brca_normal_si)
data(brca_cancer_si)
data(brca_G_mut)
#######################or using other data
#brca_normal_s <- read.table('normal_sample.txt',sep='\t',TRUE)
#brca_cancer_s <- read.table('cancer_sample.txt',sep='\t',TRUE)
#brca_normal_si <- read.table('normal_sample_mirna.txt',sep='\t',TRUE)
#brca_cancer_si <- read.table('cancer_sample_mirna.txt',sep='\t',TRUE)
#network <- read.csv("edges_4.csv" ,TRUE)
#network <- network[,-1]
#protein_info <- read.csv("protein_info.csv" ,TRUE)
#brca_G_mut <- read.table('brca_mut.txt',sep='\t',TRUE)
#rownames(brca_G_mut) <- brca_G_mut[,1]
#brca_G_mut<-brca_G_mut[,-1]
#mir_gene_network <- read.table('mir_gene_network.txt')
degs_05 <- cal_degs(brca_normal_s,brca_cancer_s,0.5)
mirna_value <- cal_degs(brca_normal_si,brca_cancer_si,0.5)
construct_net(brca_G_mut,network,protein_info,degs_05)

layer2_gene_list <- read.csv("layer2_gene_list.csv")
layer2_gene_list <- layer2_gene_list[,1]
layer2_network <- read.csv("layer2_network.csv")
gene_expscore_2 <-read.csv("gene_expscore_2.csv")
gene_expscore_3 <-read.csv("gene_expscore_3.csv")
brca_list <- MNMO(layer2_gene_list,layer2_network,brca_cancer_s,brca_G_mut,
                  mirna_value,mir_gene_network,gene_expscore_2,gene_expscore_3)

