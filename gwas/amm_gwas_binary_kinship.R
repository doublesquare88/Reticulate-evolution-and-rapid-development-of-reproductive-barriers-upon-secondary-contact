#modified from gwas_binary_kinship.RMD
#Dabao Lu 27.11.2021

library("qqman")
#install.packages("/Users/dabaosl/Dropbox/UiO/Kurs/GWAS_physalia/introduction_to_gwas/software/gMatrix_0.2.tar.gz", repos = NULL, type = "source")
library("gMatrix") ## own-made package
library("data.table")
source("/Users/dabaosl/Dropbox/UiO/Kurs/GWAS_physalia/introduction_to_gwas/software/gwas.r")
source("/Users/dabaosl/Dropbox/UiO/Kurs/GWAS_physalia/introduction_to_gwas/software/emma.r")

### Read in data  ###
setwd("/Users/dabaosl/Dropbox/UiO/Phd/Data_analysis/GWAS")
snpMatrix <- fread("illumina_1_2_3_4_hf_excl_DP3_GQ20_nomono_noindel_nodik_noTF_miss20_biallelic_maf0.01_nomissing_var_id.raw", header = TRUE)
#snpMatrix <- fread("illumina_1_2_3_4_hf_excl_DP3_GQ20_nomono_noindel_nodik_noTF_miss20_biallelic_nomissing_crosses_subset_var_id.raw", header = TRUE)
#str(snpMatrix)
SNP_INFO <- fread("illumina_1_2_3_4_hf_excl_DP3_GQ20_nomono_noindel_nodik_noTF_miss20_biallelic_maf0.01_nomissing_var_id.map")
#SNP_INFO <- fread("illumina_1_2_3_4_hf_excl_DP3_GQ20_nomono_noindel_nodik_noTF_miss20_biallelic_nomissing_crosses_subset_var_id.map")
#str(SNP_INFO)
names(SNP_INFO) <- c("Chr","SNP","cM","Pos")

#phenotypes <- read.table("cross_table_TA-1003-17-M1_phenotypes.txt", header = TRUE)
#phenotypes = read.table("cross_table_TA-1000-15-M2_phenotypes.txt", header = TRUE)
#phenotypes = read.table("cross_table_TA-1005-8-M1_phenotypes.txt", header = TRUE)
#phenotypes = read.table("cross_table_TA-1019-1-M1_phenotypes.txt", header = TRUE)
#phenotypes = read.table("cross_table_TA-1020-2-M2_phenotypes.txt", header = TRUE)
phenotypes = read.table("compatibility_all_NorAm_NB.txt", header=TRUE)

phenotypes <- phenotypes[phenotypes$id %in% snpMatrix$IID,]
phenotypes <- phenotypes[,c(1,3)]
str(phenotypes)

### Convert to Matrices ###
X <- as.matrix(snpMatrix[,-c(1:6)])
colnames(X) <- gsub("\\_[A-Z]{1}$","",colnames(X)) #to remove non numeric characters?
rownames(X) <- snpMatrix$IID

## Kinship matrix
K <- gVanRaden.2(X)
#heatmap(K,col=rev(heat.colors(75)))

str(phenotypes)
cross <- names(phenotypes[2])
cross
#a <- paste("phenotypes$",cross, sep="")
#a
Y <- as.matrix(phenotypes$compatibility_NorAm_NB)
dim(Y)
rownames(Y) <- phenotypes$id
head(Y)
### GWAS  ###
res <- amm_gwas(Y = Y, X = X, K = K, m = 1, use.SNP_INFO = TRUE)
str(res)
rownames(res) <- res$SNP
head(res)

#write out ordered res to file
res.sorted <- res[order(res$Pval,)]
str(res.sorted)
write.table(res.sorted, file="res_amm_gwas_all_inferred_all.txt")
#top1000
res.sorted.top1000 <- res.sorted[1:1000,]
write.table(res.sorted.top1000, file="res_amm_gwas_all_inferred_sorted_top1000.txt")
#bonferroni:
0.05/ncol(X)
#8.519745e-08
res.sorted.bonferroni <- filter(res.sorted, Pval < 8.519745e-08)
dim(res.sorted.bonferroni)
write.table(res.sorted.bonferroni, file="res_amm_gwas_all_inferred_sorted_bonferroni.txt")

#scf6
scf6 <- subset(res, res$Chr=="Scaffold06")
head(scf6)
scf6.sorted <- scf6[order(scf6$Pval,)]
head(scf6.sorted,20)

#scf9:
scf9 <- subset(res, res$Chr=="Scaffold09")
head(scf9)
scf9.sorted <- scf9[order(scf9$Pval,)]
head(scf9.sorted,30)

#scf5:
scf5 <- subset(res, res$Chr=="Scaffold05")
head(scf5)
scf5.sorted <- scf5[order(scf5$Pval,)]
head(scf5.sorted,30)



## rename scaffold to numeric
res$Chr <- sub("Scaffold0","",res$Chr) 
res$Chr <- sub("Scaffold","",res$Chr)
res$Chr <- as.integer(res$Chr)

##  Get SNPs with top P-values  
my_vars <- c("Chr", "Pos", "Pval")
res2 <- res[,..my_vars]
str(res2)
top10 <- head(order(res2$Pval),10)
top10
res2[top10,]
top20 <- head(order(res2$Pval),20)
top20
res2[top20,]
top30 <- head(order(res2$Pval),30)
top30
res2[top30,]

## bonferroni:
bf_log10 <- -log10(0.05/nrow(X))
bf_log10
#shouldn´t it be columns, number of SNPs instead?
0.05/ncol(X)
  #8.519745e-08

### PLOT results  ###
gwasResults <- res[,c("SNP","Chr","Pos","Pval")]
names(gwasResults) <- c("SNP","CHR","BP","P")
str(gwasResults)

#qq(gwasResults$P)
dataset="illumina_1_2_3_4_hf_excl_DP3_GQ20_nomono_noindel_nodik_noTF_miss20_biallelic_nomissing_crosses_subset"
trait_label=cross
png(paste(dataset,trait_label,"qq_amm.png",sep="_"))
qq(gwasResults$P)
dev.off()

png(paste(dataset,trait_label,"manhattan_amm.png",sep="_"))
manhattan(gwasResults, suggestiveline =FALSE , col = c("red","blue"),
          genomewideline = bf_log10, logp = TRUE, width = 800, height = 600, res = 100)
dev.off()

#manhattan(gwasResults, suggestiveline = FALSE, col = c("red","blue"))

z=qnorm(gwasResults$P/2)
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
lambda
