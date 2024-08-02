#script to make summary table for paper with gene function, p-values, enrichment tests, significance in gwas etc.

##for each gene:
#1) identify if there are SNPs falling within it
#2) count the number of SNPs
#3) record their p-values, then have it sorted after p-value
#4) separate columns for pfams and gwas

setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Incompatibility_cross_analysis")
library(tidyverse)
library(stringr)
library(stringi)
library(stats)

genome <- read.table("TA10106M1_final_MAKER.putative_function_domains.gff", sep="\t")
str(genome) #223110 obs, a lot of these are duplicated
genome.genes <- filter(genome, V3=="gene") 
str(genome.genes) #9298 obs, should represent unique genes
head(genome.genes)
pfam <- grep("Pfam", genome.genes$V9, fixed = TRUE) 
str(pfam) #4422 genes with pfam
head(genome.genes, 15)
genome.genes.pfam <- genome.genes[pfam,] #get genes with pfam annotation
colnames(genome.genes.pfam) <- c("seqid","source","type","start","end","score", "strand","phase","attributes") 
str(genome.genes.pfam)
#some genes have multiple pfams, extract these and put in new columns:
str_extract_all(genome.genes.pfam[1,9],"Pfam:PF[0-9]{5}")
pfams <- str_extract_all(genome.genes.pfam$attributes,"Pfam:PF[0-9]{5}") #extracts pfams with regular expression
length(pfams) #matches number of genes with pfams (4422)
lengths(pfams) #number of pfams in each gene
max(lengths(pfams)) #13
pfam.df <- as.data.frame(t(stri_list2matrix(pfams))) #convert to dataframe, each pfam for each gene is a column
str(pfam.df) 
genome.genes.pfam.split <- cbind(genome.genes.pfam, pfam.df)
str(genome.genes.pfam.split)
dim(genome.genes.pfam.split)

#Do the same for GOs:
go <- grep("Ontology_term", genome.genes$V9, fixed = TRUE) 
str(go) #3141 genes with go
head(genome.genes, 15)
genome.genes.go <- genome.genes[go,] #get genes with pfam annotation
colnames(genome.genes.go) <- c("seqid","source","type","start","end","score", "strand","phase","attributes") 
str(genome.genes.go)
dim(genome.genes.go)

str(genome.genes.go)
gos <- str_extract_all(genome.genes.go$attributes,"GO:[0-9]{7}") #extracts gos with regular expression
gos
length(gos) #matches number of genes with gos (3141)
lengths(gos) #number of pfams in each gene
max(lengths(gos)) #17
gos.df <- as.data.frame(t(stri_list2matrix(gos))) #convert to dataframe, each go for each gene is a column
str(gos.df) 
genome.genes.go.split <- cbind(genome.genes.go, gos.df)
str(genome.genes.go.split)
dim(genome.genes.go.split)

#######################
#                     #
#   GWAS              #
#                     #
#######################
###read in SNPs in GWAS analysis that are significant after bonferroni correction
gwas.all.bonferroni <- read.table("GWAS/output_amm/res_amm_gwas_all_inferred_sorted_bonferroni.txt")
str(gwas.all.bonferroni)
head(gwas.all.bonferroni,20)
#gwas.all.bonferroni <- head(gwas.all.bonferroni, 50)

#Find pfams for significant gwas genes
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
p.list <- c()
p.list.num <- c()
p.list.num.sort <- c()
lowest.p.list <- c()
for(gene in 1:nrow(genome.genes.pfam.split)) {
  gene.scf <- genome.genes.pfam[gene,1]
  gene.start <- genome.genes.pfam[gene,4]
  gene.end <- genome.genes.pfam[gene,5]
  gene.anno <- genome.genes.pfam[gene,9]
  snp <- filter(gwas.all.bonferroni, Chr==gene.scf & Pos > gene.start & Pos < gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")  
  } else {
    print("snp in gene")
    print(snp$Pval)
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    lowest.p <- c(snp$Pval[1])
    lowest.p.list <- c(lowest.p, lowest.p.list)
    p.list <- c(paste(snp$Pval, collapse = ' '), p.list)
    pval.sorted <- sort(snp$Pval)
    p.list.num <- c(snp$Pval, p.list)
    p.list.num.sort <- c(sort(snp$Pval), p.list.num.sort)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}

#can you assume that first p-value in each sublist will always be the lowest because input snp gwas list is sorted in this order?
#seems to be valid here
lowest.p.list
length(lowest.p.list)
p.list.num
p.list.num.sort
test <- sort(p.list.num)
str(test)

length(gene.list) #407 genes that contain significant GWAS snps
sum(num.snps) #874 snps (out of a total 1144 snps within genes)
genes.in.snps <- cbind(scf.list, gene.list, end.list, num.snps, lowest.p.list, p.list, gene.info)
genes.in.snps <- as.data.frame(genes.in.snps)
genes.in.snps$lowest.p.list <- as.numeric(genes.in.snps$lowest.p.list)
genes.in.snps <- genes.in.snps[order(genes.in.snps$lowest.p.list, decreasing = F),]
head(genes.in.snps)

#parse out pfam name:
genes.in.snps.pfams <- str_extract_all(genes.in.snps$gene.info,"PF[0-9]{5}") #extracts pfams
length(genes.in.snps.pfams) #matches number of genes in snps (407)
pfams.string <-  sapply(genes.in.snps.pfams, paste,collapse = " ")
pfams.string
length(pfams.string)
sum(lengths(genes.in.snps.pfams)) #total 903 pfams
str(genes.in.snps.pfams)
genes.in.snps.pfams.df <- cbind(genes.in.snps, pfams.string)
str(genes.in.snps.pfams.df)
#write_delim(genes.in.snps.pfams, "test_gwas.txt")

#parse out gene id:
genes.in.snps.geneid <- str_extract_all(genes.in.snps$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.snps.geneid)
genes.in.snps.geneid <- unlist(genes.in.snps.geneid)
genes.in.snps.geneid

genes.in.snps.pfams.df <- cbind(genes.in.snps.geneid, genes.in.snps, pfams.string)
#write_delim(genes.in.snps.pfams.df, "test_gwas.txt")

#No also get GOs for the gwas genes:
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
p.list <- c()
p.list.num <- c()
p.list.num.sort <- c()
lowest.p.list <- c()
for(gene in 1:nrow(genome.genes.go)) {
  gene.scf <- genome.genes.go[gene,1]
  gene.start <- genome.genes.go[gene,4]
  gene.end <- genome.genes.go[gene,5]
  gene.anno <- genome.genes.go[gene,9]
  snp <- filter(gwas.all.bonferroni, Chr==gene.scf & Pos > gene.start & Pos < gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")  
  } else {
    print("snp in gene")
    print(snp$Pval)
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    lowest.p <- c(snp$Pval[1])
    lowest.p.list <- c(lowest.p, lowest.p.list)
    p.list <- c(paste(snp$Pval, collapse = ' '), p.list)
    pval.sorted <- sort(snp$Pval)
    p.list.num <- c(snp$Pval, p.list)
    p.list.num.sort <- c(sort(snp$Pval), p.list.num.sort)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}
length(scf.list)
length(gene.list)
length(num.snps)
length(p.list.num)
length(p.list)

length(gene.list) #326 go genes that contain significant GWAS snps
sum(num.snps) #723 snps 
genes.in.snps <- cbind(scf.list, gene.list, end.list, num.snps, lowest.p.list, p.list, gene.info)
genes.in.snps <- as.data.frame(genes.in.snps)
genes.in.snps$lowest.p.list <- as.numeric(genes.in.snps$lowest.p.list)
genes.in.snps <- genes.in.snps[order(genes.in.snps$lowest.p.list, decreasing = F),]
head(genes.in.snps)
#parse out go name:
genes.in.snps.go <- str_extract_all(genes.in.snps$gene.info,"GO:[0-9]{7}") #extracts pfams
length(genes.in.snps.go) #matches number of genes in snps (326)
gos.string <-  sapply(genes.in.snps.go, paste,collapse = " ")
gos.string
length(gos.string)
sum(lengths(genes.in.snps.go)) #total 1147 gos
str(genes.in.snps.go)
genes.in.snps.go.df <- cbind(genes.in.snps, gos.string)
str(genes.in.snps.go.df)
#parse out and add gene id (actually should have done this at the very beginning, but anyways):
genes.in.snps.geneid <- str_extract_all(genes.in.snps$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.snps.geneid) #326
genes.in.snps.geneid <- unlist(genes.in.snps.geneid)
genes.in.snps.geneid
genes.in.snps.go.df <- cbind(genes.in.snps.geneid, genes.in.snps, gos.string)
write_delim(genes.in.snps.go.df, "test_gwas_go.txt")

###Now combine pfam into one table:
dim(genes.in.snps.go.df)
dim(genes.in.snps.pfams.df)
str(genes.in.snps.go.df)
str(genes.in.snps.pfams.df)
genes.in.snps.pfams.df.test <- genes.in.snps.pfams.df[2:40,]

genes.in.snps.go.df$pfams.string <- "NA"
genes.in.snps.pfams.df$gos.string <- "NA"
genes.in.snps.go.df <- genes.in.snps.go.df[,1:9]
genes.in.snps.pfams.df <- genes.in.snps.pfams.df[,1:9]

genes.gwas.pfams.go <- merge(x = genes.in.snps.pfams.df, y = genes.in.snps.go.df, by = "genes.in.snps.geneid",  all = FALSE)
genes.gwas.pfams.go <- full_join(genes.in.snps.pfams.df,genes.in.snps.go.df, 
                                 by = "genes.in.snps.geneid")
str(genes.gwas.pfams.go)
dim(genes.gwas.pfams.go) #407
table(genes.in.snps.pfams.df$genes.in.snps.geneid %in% genes.gwas.pfams.go$genes.in.snps.geneid) #check that all pfams gene ids made it:
table(genes.gwas.pfams.go$genes.in.snps.geneid %in% genes.in.snps.pfams.df$genes.in.snps.geneid)
table(genes.in.snps.go.df$genes.in.snps.geneid %in% genes.gwas.pfams.go$genes.in.snps.geneid) #check that all pfams gene ids made it:
table(genes.gwas.pfams.go$genes.in.snps.geneid %in% genes.in.snps.go.df$genes.in.snps.geneid)

#extract some more central gene info:
notes <- str_extract_all(genes.gwas.pfams.go$gene.info.x,"Note=(.*?);") 
notes <- unlist(notes)
notes
genes.gwas.pfams.go$notes <- notes
str(genes.gwas.pfams.go)
head(genes.gwas.pfams.go)
write_csv2(genes.gwas.pfams.go, "gwas_pfam_go.csv")

#######################
#                     #   
#     SCOPA SNPs      #
#                     #
#######################

scopa_all <- read.table("SCOPA_crosses_result_simple_sorted.txt", header=T)
str(scopa_all) #641873 snps
#apply bonferroni correction: 
0.05/641873 #7.79e-08
scopa.all.bonferroni <- filter(scopa_all, P.value < 7.79e-08)
str(scopa.all.bonferroni) #6075 snps, but numerical scaffold names, convert these to match with annotation file:
scopa.all.bonferroni$Chromosome <- gsub("12", "twelve",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("11", "eleven",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("10", "ten",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("9", "Scaffold09",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("8", "Scaffold08",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("7", "Scaffold07",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("6", "Scaffold06",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("5", "Scaffold05",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("4", "Scaffold04",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("3", "Scaffold03",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("2", "Scaffold02",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("1", "Scaffold01",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("twelve", "Scaffold12",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("eleven", "Scaffold11",scopa.all.bonferroni$Chromosome)
scopa.all.bonferroni$Chromosome <- gsub("ten", "Scaffold10",scopa.all.bonferroni$Chromosome)
table(scopa.all.bonferroni$Chromosome)        
str(scopa.all.bonferroni)
head(scopa.all.bonferroni,100)

#Find pfams for significant scopa genes
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
p.list <- c()
p.list.num <- c()
p.list.num.sort <- c()
lowest.p.list <- c()
for(gene in 1:nrow(genome.genes.pfam.split)) {
  gene.scf <- genome.genes.pfam[gene,1]
  gene.start <- genome.genes.pfam[gene,4]
  gene.end <- genome.genes.pfam[gene,5]
  gene.anno <- genome.genes.pfam[gene,9]
  snp <- filter(scopa.all.bonferroni, Chromosome==gene.scf & Position > gene.start & Position < gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")  
  } else {
    print("snp in gene")
    print(snp$P.value)
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Position), num.snps)
    lowest.p <- c(snp$P.value[1])
    lowest.p.list <- c(lowest.p, lowest.p.list)
    p.list <- c(paste(snp$P.value, collapse = ' '), p.list)
    pval.sorted <- sort(snp$P.value)
    p.list.num <- c(snp$P.value, p.list)
    p.list.num.sort <- c(sort(snp$P.value), p.list.num.sort)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}

#can you assume that first p-value in each sublist will always be the lowest because input snp gwas list is sorted in this order?
#seems to be valid here
lowest.p.list
length(lowest.p.list)
p.list.num
p.list.num.sort
test <- sort(p.list.num)
str(test)

length(gene.list) #1580 genes that contain significant GWAS snps
sum(num.snps) #4399 snps 
genes.in.snps <- cbind(scf.list, gene.list, end.list, num.snps, lowest.p.list, p.list, gene.info)
genes.in.snps <- as.data.frame(genes.in.snps)
genes.in.snps$lowest.p.list <- as.numeric(genes.in.snps$lowest.p.list)
genes.in.snps <- genes.in.snps[order(genes.in.snps$lowest.p.list, decreasing = F),]
dim(genes.in.snps)

#parse out pfam name:
genes.in.snps.pfams <- str_extract_all(genes.in.snps$gene.info,"PF[0-9]{5}") #extracts pfams
length(genes.in.snps.pfams) #matches number of genes in snps (407)
pfams.string <-  sapply(genes.in.snps.pfams, paste,collapse = " ")
pfams.string
length(pfams.string)
sum(lengths(genes.in.snps.pfams)) #total 3179 pfams
length(genes.in.snps.pfams)
genes.in.snps.pfams.df <- cbind(genes.in.snps, pfams.string)
str(genes.in.snps.pfams.df)
#write_delim(genes.in.snps.pfams, "test_gwas.txt")

#parse out gene id:
genes.in.snps.geneid <- str_extract(genes.in.snps$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.snps.geneid) #1580
genes.in.snps.geneid <- unlist(genes.in.snps.geneid)
length(genes.in.snps.geneid)
genes.in.snps.pfams.df <- cbind(genes.in.snps.geneid, genes.in.snps, pfams.string)
#write_delim(genes.in.snps.pfams.df, "test_gwas.txt")

#No also get GOs for the scopa genes:
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
p.list <- c()
p.list.num <- c()
p.list.num.sort <- c()
lowest.p.list <- c()
for(gene in 1:nrow(genome.genes.go)) {
  gene.scf <- genome.genes.go[gene,1]
  gene.start <- genome.genes.go[gene,4]
  gene.end <- genome.genes.go[gene,5]
  gene.anno <- genome.genes.go[gene,9]
  snp <- filter(scopa.all.bonferroni, Chromosome==gene.scf & Position > gene.start & Position < gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")  
  } else {
    print("snp in gene")
    print(snp$P.value)
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Position), num.snps)
    lowest.p <- c(snp$P.value[1])
    lowest.p.list <- c(lowest.p, lowest.p.list)
    p.list <- c(paste(snp$P.value, collapse = ' '), p.list)
    pval.sorted <- sort(snp$P.value)
    p.list.num <- c(snp$P.value, p.list)
    p.list.num.sort <- c(sort(snp$P.value), p.list.num.sort)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}

length(gene.list) #1188 go genes that contain significant GWAS snps
sum(num.snps) #3485 snps 
genes.in.snps <- cbind(scf.list, gene.list, end.list, num.snps, lowest.p.list, p.list, gene.info)
genes.in.snps <- as.data.frame(genes.in.snps)
genes.in.snps$lowest.p.list <- as.numeric(genes.in.snps$lowest.p.list)
genes.in.snps <- genes.in.snps[order(genes.in.snps$lowest.p.list, decreasing = F),]
head(genes.in.snps)
#parse out go name:
genes.in.snps.go <- str_extract_all(genes.in.snps$gene.info,"GO:[0-9]{7}") #extracts gos
length(genes.in.snps.go) #matches number of genes in snps (326)
gos.string <-  sapply(genes.in.snps.go, paste,collapse = " ")
gos.string
length(gos.string)
sum(lengths(genes.in.snps.go)) #total 1188 gos
str(genes.in.snps.go)
genes.in.snps.go.df <- cbind(genes.in.snps, gos.string)
str(genes.in.snps.go.df)
#parse out and add gene id (actually should have done this at the very beginning, but anyways):
genes.in.snps.geneid <- str_extract(genes.in.snps$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.snps.geneid) #1188
genes.in.snps.geneid <- unlist(genes.in.snps.geneid)
genes.in.snps.geneid
genes.in.snps.go.df <- cbind(genes.in.snps.geneid, genes.in.snps, gos.string)
#write_delim(genes.in.snps.go.df, "test_gwas_go.txt")

###Now combine pfam and go into one table:
dim(genes.in.snps.go.df)
dim(genes.in.snps.pfams.df)
str(genes.in.snps.go.df)
str(genes.in.snps.pfams.df)
#genes.in.snps.pfams.df.test <- genes.in.snps.pfams.df[2:40,]
genes.scopa.pfams.go <- full_join(genes.in.snps.pfams.df,genes.in.snps.go.df, 
                                 by = "genes.in.snps.geneid")
str(genes.scopa.pfams.go)
dim(genes.scopa.pfams.go) #1580
table(genes.in.snps.pfams.df$genes.in.snps.geneid %in% genes.scopa.pfams.go$genes.in.snps.geneid) #check that all pfams gene ids made it:
table(genes.scopa.pfams.go$genes.in.snps.geneid %in% genes.in.snps.pfams.df$genes.in.snps.geneid)
table(genes.in.snps.go.df$genes.in.snps.geneid %in% genes.scopa.pfams.go$genes.in.snps.geneid) #check that all pfams gene ids made it:
table(genes.scopa.pfams.go$genes.in.snps.geneid %in% genes.in.snps.go.df$genes.in.snps.geneid)

#extract some more central gene info:
notes <- str_extract(genes.scopa.pfams.go$gene.info.x,"Note=(.*?);") 
notes <- unlist(notes)
length(notes)
genes.scopa.pfams.go$notes <- notes
str(genes.scopa.pfams.go)
head(genes.scopa.pfams.go)
write_csv2(genes.scopa.pfams.go, "scopa_pfam_go.csv")

###############################
#                             #
#         DeltaFst            #
#                             #
###############################

delta_fst.all <- read.table("Delta_df_all_outliers.txt", header=T)
str(delta_fst.all)
delta_fst.noram <- read.table("Delta_df_NoramA_NoramB.txt", header=T)
str(delta_fst.noram) #has one addtional comparison not in delta_fst.all: D.NoramAE_NoramBE.NoramAS_NoramBW

#add D.NoramAE_NoramBE.NoramAS_NoramBW to delta.fst.all:
delta_fst.noram.add <- delta_fst.noram[,c()]
delta_fst_all_test <- full_join(delta_fst.all, delta_fst.noram)
str(delta_fst_all_test) #1631 < 1522(fst all) + 197(noram) # this is good, indicates that overlaps were identified
delta_fst_all_test$outlier.D.NoramAE_NoramBE.NoramAS_NoramBW
delta_fst_all_test <- delta_fst_all_test %>% mutate(outlier.D.NoramAE_NoramBE.NoramAS_NoramBW = ifelse(is.na(outlier.D.NoramAE_NoramBE.NoramAS_NoramBW), "background", outlier.D.NoramAE_NoramBE.NoramAS_NoramBW))
delta_fst_all_test$outlier.D.NoramAE_NoramBE.NoramAS_NoramBW
delta_fst_all_test$D.NoramAE_NoramBE.NoramAS_NoramBW
delta_fst_all <- delta_fst_all_test
#delta_fst_all <- delta_fst_all[1:200,]
head(delta_fst_all)

genome.genes.pfam.sub <- genome.genes.pfam[1:200,]
#names(delta_fst_all)[max.col(delta_fst_all == "outlier")]
#names(delta_fst_all)[which(delta_fst_all == "outlier", arr.ind = T)[,"col"]]
#b <- names(delta_fst_all)[which(delta_fst_all == "outlier", arr.ind = T)[,"col"]]
b <- lapply(apply(delta_fst_all,1, function(x) which(x=="outlier")),names)
length(b)
lengths(b)
b
str(b)
num_outliers <- lengths(b)
c <- sapply(b, unique) #make sure there are no duplicates: 
length(c) #1631
lengths(c) #appears to be the same as b

outliers.string <-  sapply(b, paste,collapse = ":")
delta_fst_all <- cbind(delta_fst_all, outliers.string)
delta_fst_all <- cbind(delta_fst_all, num_outliers)
str(delta_fst_all)

#how to get most extreme fst outlier value and know which comparison it is?
d <- head(delta_fst_all) 
dim(d)
d <- d[,6:42]
colnames(d)[apply(d,1,which.max)]
d
#max(delta_fst_all_test)
#df[1, which.max(),]
colnames(d)[apply(d,1,which.min)]
#problem is that outlier are within each comparison, so what has been identified as a outlier might not be the
#max value across comparisons.

###PFAM
gene.list <- c()
scf.list <- c()
snp.list <- c()
gene.info <- c()
gene.list <- c()
end.list <- c()
outlier.list <- c()
num.outliers.list <- c()
delta_window_list <- c()
for(gene in 1:nrow(genome.genes.pfam)) {
  gene.scf <- genome.genes.pfam[gene,1]
  gene.start <- genome.genes.pfam[gene,4]
  gene.end <- genome.genes.pfam[gene,5]
  gene.anno <- genome.genes.pfam[gene,9]
  window <- filter(delta_fst_all, Scaffold==gene.scf & gene.end >= Start & gene.start <= End)
  if (dim(window)[1] == 0) {
    print("no overlap")  
  } else {
    print("overlap")
    #print(window)
    scf.list <- c(gene.scf, scf.list)
    gene.info <- c(gene.anno, gene.info)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    #delta_window <- window$Start
    print(window$outliers.string)
    outlier <- unique(window$outliers.string)
    outlier <- paste(window$outliers.string, collapse = ":")
    #outlier <- sapply(window$outliers.string, paste, collapse = " ")
    outlier.list <- c(outlier, outlier.list)
    num.outliers <- sum(window$num_outliers)
    num.outliers.list <- c(num.outliers, num.outliers.list)
    #outlier <- lapply(apply(window,1, function(x) which(x=="outlier")),names)
    #outlier <- apply(window,1, function(x) which(x=="outlier"))
    #outlier.list <- c(outlier, outlier.list)
  }
}

window
test <- head(delta_fst_all)
dim(test)
a <- lapply(test$outliers.string, paste, collapse = ":") 
class(test$outliers.string)
a <- test$outliers.string
class(a)
b <- paste(a, collapse = " ")
length(b)
paste(test$outliers.string, collapse = ":")

length(outlier.list) #1357
length(scf.list) #1357 
length(gene.list) #1357 genes with pfam that overlap with deltafst windows
length(num.outliers.list) #1357
length(end.list) #1357
length(gene.info) #1357
a <- cbind(scf.list, gene.list, end.list, gene.info, outlier.list, num.outliers.list)
dim(a)
genes.in.delta <- as.data.frame(a)
dim(genes.in.delta)
#now also parse out each individual pfam in these genes:
genes.in.delta.pfams <- str_extract_all(genes.in.delta$gene.info,"Pfam:PF[0-9]{5}") #extracts pfams
length(genes.in.delta.pfams) #matches number of genes in snps (1357)
pfams.string <-  sapply(genes.in.delta.pfams, paste,collapse = " ")
length(pfams.string)
sum(lengths(genes.in.delta.pfams)) #total 2351 pfams
genes.in.delta.pfams.df <- cbind(genes.in.delta, pfams.string)
#parse out gene id:
genes.in.delta.geneid <- str_extract(genes.in.delta$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.delta.geneid) #1357
genes.in.delta.geneid <- unlist(genes.in.delta.geneid)
length(genes.in.delta.geneid) #1357
genes.in.delta.pfams.df <- cbind(genes.in.delta.geneid, genes.in.delta, pfams.string)
str(genes.in.delta.pfams.df)

###GO
gene.list <- c()
end.list <- c()
scf.list <- c()
snp.list <- c()
outlier.windows <- c()
num.outliers.list <- c()
gene.info <- c()
outlier.list <- c()
for(gene in 1:nrow(genome.genes.go)) {
  gene.scf <- genome.genes.go[gene,1]
  gene.start <- genome.genes.go[gene,4]
  gene.end <- genome.genes.go[gene,5]
  gene.anno <- genome.genes.go[gene,9]
  window <- filter(delta_fst_all, Scaffold==gene.scf & gene.end >= Start & gene.start <= End)
  if (dim(window)[1] == 0) {
    print("no overlap")  
  } else {
    print("overlap")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    gene.info <- c(gene.anno, gene.info)
    outlier <- unique(window$outliers.string)
    outlier <- paste(window$outliers.string, collapse = ":")
    outlier.list <- c(outlier, outlier.list)
    num.outliers <- sum(window$num_outliers)
    num.outliers.list <- c(num.outliers, num.outliers.list)
  }
}
str(delta_fst_all)

length(scf.list) #974
length(gene.list) #974 genes with go that overlap with deltafst windows
length(gene.info) #974
length(outlier.list) #974
length(end.list) #974
length(num.outliers.list) #974
head(outlier.list)
a <- cbind(scf.list, gene.list, end.list, gene.info, outlier.list, num.outliers.list)

head(a)
genes.in.delta <- as.data.frame(a)
str(genes.in.delta)
#now also parse out each individual go in these 974 genes:
genes.in.delta.go <- str_extract_all(genes.in.delta$gene.info,"GO:[0-9]{7}") #extracts gos
length(genes.in.delta.go) #matches number of genes in snps (974)
go.string <-  sapply(genes.in.delta.go, paste,collapse = " ")
length(go.string)
sum(lengths(genes.in.delta.go)) #total 2992 gos
genes.in.delta.go.df <- cbind(genes.in.delta, go.string)
#parse out gene id:
genes.in.delta.geneid <- str_extract(genes.in.delta$gene.info,"ID=TA10106M1_[0-9]{5}") #extracts gene IDs
length(genes.in.delta.geneid) #902
genes.in.delta.geneid <- unlist(genes.in.delta.geneid)
length(genes.in.delta.geneid)
genes.in.delta.go.df <- cbind(genes.in.delta.geneid, genes.in.delta, go.string)
str(genes.in.delta.go.df)

#combine go and pfam:
dim(genes.in.delta.go.df)
dim(genes.in.delta.pfams.df)
str(genes.in.delta.go.df)
str(genes.in.delta.pfams.df)
#genes.in.snps.pfams.df.test <- genes.in.snps.pfams.df[2:40,]
genes.delta.pfams.go <- full_join(genes.in.delta.pfams.df,genes.in.delta.go.df, 
                                  by = "genes.in.delta.geneid")
str(genes.delta.pfams.go)
dim(genes.delta.pfams.go) #1357
table(genes.in.delta.pfams.df$genes.in.delta.geneid %in% genes.delta.pfams.go$genes.in.delta.geneid) #check that all pfams gene ids made it:
table(genes.delta.pfams.go$genes.in.delta.geneid %in% genes.in.delta.pfams.df$genes.in.delta.geneid)
table(genes.in.delta.go.df$genes.in.delta.geneid %in% genes.delta.pfams.go$genes.in.delta.geneid) #check that all pfams gene ids made it:
table(genes.delta.pfams.go$genes.in.delta.geneid %in% genes.in.delta.go.df$genes.in.delta.geneid)

#extract some more central gene info:
notes <- str_extract(genes.delta.pfams.go$gene.info.x,"Note=(.*?);") 
notes <- unlist(notes)
length(notes)
genes.delta.pfams.go$notes <- notes
str(genes.delta.pfams.go)
head(genes.delta.pfams.go)
dim(genes.delta.pfams.go)

genes.delta.pfams.go$gene.list.x <- as.numeric(genes.delta.pfams.go$gene.list.x)
genes.delta.pfams.go.ordered <- genes.delta.pfams.go[order(genes.delta.pfams.go$scf.list.x, genes.delta.pfams.go$gene.list.x),]
head(genes.delta.pfams.go.ordered, 20)

#try to extract and parse out each outlier comparison name as you did for pfam genes to deal with duplicates:
head(genes.delta.pfams.go.ordered$outlier.list.x)
old_outlier <- genes.delta.pfams.go.ordered$outlier.list.x

deltas <- str_extract_all(genes.delta.pfams.go.ordered$outlier.list.x,"outlier.D(.*?):")
length(deltas)
lengths(deltas)
my_delta <- lengths(deltas)
sum(lengths(deltas)) #2659
deltas.unique <- sapply(deltas, unique)
adj.num.outliers <- lengths(deltas.unique)
adj.num.outliers
sum(lengths(deltas.unique)) #1934
#get indices where there was only one window and you weren´t able to extract:
zero <- which(adj.num.outliers==0)
adj.num.outliers
a <- lapply(deltas.unique, paste, collapse = ":") 
b <- unlist(a)
str(b)

new_list <- list()
n <- 0
for(i in my_delta){
  #print(i)
  n <- n + 1
  if(i == 0){
    new <- old_outlier[n] 
  }
  else{
    new <- b[n]
  }
  new_list <- c(new, new_list)
}
 
str(new_list)  
length(new_list) 
lengths(new_list)
sum(lengths(new_list))
addendum <- unlist(new_list)
str(addendum) 
addendum[1]
adj.num.outliers

addendum.df <- cbind(addendum,adj.num.outliers)
write_csv2(as.data.frame(addendum.df), "delta_outliers_corrected_info.csv")

#write_csv2(genes.delta.pfams.go.ordered, "delta_pfam_go_end_counts_corrected.csv")
#write_csv2(genes.delta.pfams.go.ordered, "delta_pfam_go_end_counts.csv")

###Add pos selection info:
setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Relate/Detect_positive_selection")
pos.sel.noramA <- read.table("noramA.sele", header=T)
str(pos.sel.noramA)
head(sort(10^pos.sel.noramA$when_mutation_has_freq2),20)
quantile(10^pos.sel.noramA$when_mutation_has_freq2, probs=0.02)
hist(10^pos.sel.noramA$when_mutation_has_freq2)
#need to filter out entries for which p-values cannot be calculated that has been set to -10log=1
pos.sel.noramA.filtered <- filter(pos.sel.noramA, when_mutation_has_freq2 < 1)
dim(pos.sel.noramA.filtered)
hist(10^pos.sel.noramA.filtered$when_mutation_has_freq2)
pos.sel.noramA.filtered$pvalue <- 10^pos.sel.noramA.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.noramA.filtered$pvalue < 0.05/length(pos.sel.noramA.filtered$pvalue)) #4.908939e-07, only 1 snp qualifies for this...
pos.sel.noramA.filtered$BH <- p.adjust(pos.sel.noramA.filtered$pvalue, method = "BH") 
table(pos.sel.noramA.filtered$BH < 0.05) #2
pos.sel.noramA.filtered.BH.sign <- filter(pos.sel.noramA.filtered, BH < 0.05)
str(pos.sel.noramA.filtered.BH.sign)

pos.sel.noramB <- read.table("noramBE.sele", header=T)
head(sort(10^pos.sel.noramB$when_mutation_has_freq2),20)
quantile(10^pos.sel.noramB$when_mutation_has_freq2, probs=0.02)
pos.sel.noramB.filtered <- filter(pos.sel.noramB, when_mutation_has_freq2 < 1)
pos.sel.noramB.filtered$pvalue <- 10^pos.sel.noramB.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.noramB.filtered$pvalue < 0.05/length(pos.sel.noramB.filtered$pvalue)) #6
pos.sel.noramB.filtered$BH <- p.adjust(pos.sel.noramB.filtered$pvalue, method = "BH") 
table(pos.sel.noramB.filtered$BH < 0.05) #53
pos.sel.noramB.filtered.BH.sign <- filter(pos.sel.noramB.filtered, BH < 0.05)
str(pos.sel.noramB.filtered.BH.sign)

pos.sel.eurasia <- read.table("Eurasia.sele", header=T)
head(sort(10^pos.sel.eurasia$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.eurasia$when_mutation_has_freq2),20)
pos.sel.eurasia.filtered <- filter(pos.sel.eurasia, when_mutation_has_freq2 < 1)
pos.sel.eurasia.filtered$pvalue <- 10^pos.sel.eurasia.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.eurasia.filtered$pvalue < 0.05/length(pos.sel.eurasia.filtered$pvalue)) #37 snp qualifies for this...
pos.sel.eurasia.filtered$BH <- p.adjust(pos.sel.eurasia.filtered$pvalue, method = "BH") 
table(pos.sel.eurasia.filtered$BH < 0.05) #BH correction gives more than 495 significant snp..
pos.sel.eurasia.filtered.BH.sign <- filter(pos.sel.eurasia.filtered, BH < 0.05)
str(pos.sel.eurasia.filtered.BH.sign)

pos.sel.euro <- read.table("Euro.sele", header=T)
head(sort(10^pos.sel.euro$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.euro$when_mutation_has_freq2),20)
pos.sel.euro.filtered <- filter(pos.sel.euro, when_mutation_has_freq2 < 1)
pos.sel.euro$pvalue <- 10^pos.sel.euro$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.euro.filtered$pvalue < 0.05/length(pos.sel.euro.filtered$pvalue)) #0 
pos.sel.euro.filtered$BH <- p.adjust(pos.sel.euro.filtered$pvalue, method = "BH") 
table(pos.sel.euro.filtered$BH < 0.05) #0 snps
pos.sel.euro.filtered.BH.sign <- filter(pos.sel.euro.filtered, BH < 0.05)
str(pos.sel.euro.filtered.BH.sign)

pos.sel.asiaSW <- read.table("AsiaSW.sele", header=T)
head(sort(10^pos.sel.asiaSW$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.asiaSW$when_mutation_has_freq2),20)
pos.sel.asiaSW.filtered <- filter(pos.sel.asiaSW, when_mutation_has_freq2 < 1)
pos.sel.asiaSW.filtered$pvalue <- 10^pos.sel.asiaSW.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.asiaSW.filtered$pvalue < 0.05/length(pos.sel.asiaSW.filtered$pvalue)) #0 
pos.sel.asiaSW.filtered$BH <- p.adjust(pos.sel.asiaSW.filtered$pvalue, method = "BH") 
table(pos.sel.asiaSW.filtered$BH < 0.05) #0 snps
dim(pos.sel.asiaSW.filtered)
sort(pos.sel.asiaSW.filtered$pvalue)

pos.sel.asiaSC <- read.table("AsiaSichuan.sele", header=T)
head(sort(10^pos.sel.asiaSC$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.asiaSC$when_mutation_has_freq2),20)
pos.sel.asiaSC.filtered <- filter(pos.sel.asiaSC, when_mutation_has_freq2 < 1)
pos.sel.asiaSC.filtered$pvalue <- 10^pos.sel.asiaSC.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.asiaSC.filtered$pvalue < 0.05/length(pos.sel.asiaSC.filtered$pvalue)) #0 
pos.sel.asiaSC.filtered$BH <- p.adjust(pos.sel.asiaSC.filtered$pvalue, method = "BH") 
table(pos.sel.asiaSC.filtered$BH < 0.05) #0 snps
dim(pos.sel.asiaSC.filtered)
sort(pos.sel.asiaSC.filtered$pvalue)

pos.sel.asiaNE <- read.table("AsiaNE.sele", header=T)
head(sort(10^pos.sel.asiaNE$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.asiaNE$when_mutation_has_freq2),20)
pos.sel.asiaNE.filtered <- filter(pos.sel.asiaNE, when_mutation_has_freq2 < 1)
pos.sel.asiaNE.filtered$pvalue <- 10^pos.sel.asiaNE.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.asiaNE.filtered$pvalue < 0.05/length(pos.sel.asiaNE.filtered$pvalue)) #0 
pos.sel.asiaNE.filtered$BH <- p.adjust(pos.sel.asiaNE.filtered$pvalue, method = "BH") 
table(pos.sel.asiaNE.filtered$BH < 0.05) #0 snps
dim(pos.sel.asiaNE.filtered)
sort(pos.sel.asiaNE.filtered$pvalue)

pos.sel.asiaN <- read.table("AsiaN.sele", header=T)
head(sort(10^pos.sel.asiaN$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.asiaN$when_mutation_has_freq2),20)
pos.sel.asiaN.filtered <- filter(pos.sel.asiaN, when_mutation_has_freq2 < 1)
pos.sel.asiaN.filtered$pvalue <- 10^pos.sel.asiaN.filtered$when_mutation_has_freq2 
#bonferroni correction of p-values:
table(pos.sel.asiaN.filtered$pvalue < 0.05/length(pos.sel.asiaN.filtered$pvalue)) #0 
pos.sel.asiaN.filtered$BH <- p.adjust(pos.sel.asiaN.filtered$pvalue, method = "BH") 
table(pos.sel.asiaN.filtered$BH < 0.05) #0 snps
dim(pos.sel.asiaN.filtered)
sort(pos.sel.asiaN.filtered$pvalue)
head(pos.sel.asiaN.filtered)

###################
#                 #
#     Delta       #
#                 #
###################
genes.delta.pfams.go.ordered <- read.csv2("delta_pfam_go_end_counts.csv")
str(genes.delta.pfams.go.ordered)

num.snps <- c(length(snp$Pos), num.snps)
lowest.p <- c(snp$Pval[1])
lowest.p.list <- c(lowest.p, lowest.p.list)
p.list <- c(paste(snp$Pval, collapse = ' '), p.list)
pval.sorted <- sort(snp$Pval)
p.list.num <- c(snp$Pval, p.list)
p.list.num.sort <- c(sort(snp$Pval), p.list.num.sort)

gene.list <- c()
end.list <- c()
scf.list <- c()
snp.list <- c()
p.list <- c()
for(gene in 1:nrow(genes.delta.pfams.go.ordered)){
  gene.scf <- genes.delta.pfams.go.ordered[gene,2]
  gene.start <- genes.delta.pfams.go.ordered[gene,3]
  gene.end <- genes.delta.pfams.go.ordered[gene,4]
  #snp <- filter(pos.sel.noramB.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.noramA.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.eurasia.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  snp <- filter(pos.sel.euro.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
      if (dim(snp)[1] == 0) {
    print("no snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- 0
    snp.list <- c(snp_gene, snp.list)
    pvalue <- NA
    p.list <- c(pvalue, p.list)
  } else {
    print("snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- length(snp$pos)
    snp.list <- c(snp_gene, snp.list)
    p.list <- c(paste(snp$BH, collapse = ' '), p.list)
  }
}

length(scf.list)
length(gene.list)
length(end.list)
length(snp.list)
length(p.list)
snp.list
p.list

noramA <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
noramB <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
eurasia <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
euro <- cbind(scf.list, gene.list, end.list, snp.list, p.list)

noramA <- as.data.frame(noramA)
colnames(noramA) <- c("scf","start", "end", "noramA_sign_pos_sel_snps", "noramA_p-value")
noramB <- as.data.frame(noramB)
colnames(noramB) <- c("scf","start", "end", "noramB_sign_pos_sel_snps", "noramB_p-value")
eurasia <- as.data.frame(eurasia)
colnames(eurasia) <- c("scf","start", "end", "eurasia_sign_pos_sel_snps", "eurasia_p-value")
euro <- as.data.frame(euro)
colnames(euro) <- c("scf","start", "end", "euro_sign_pos_sel_snps", "euro_p-value")

all_pops <- cbind(noramA, noramB, eurasia, euro)
str(all_pops)
all_pops$start <- as.numeric(all_pops$start)

delta_all_pops.ordered <- all_pops[order(all_pops$scf, all_pops$start),]
head(delta_all_pops.ordered)

setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Incompatibility_cross_analysis")
write_csv2(delta_all_pops.ordered, "delta_all_pops.ordered.csv")

###################
#                 #
#     GWAS        #
#                 #
###################
genes.gwas.pfams.go <- read.csv2("gwas_pfam_go.csv")
str(genes.gwas.pfams.go)

gene.list <- c()
end.list <- c()
scf.list <- c()
snp.list <- c()
p.list <- c()
gene.id.list <- c()
for(gene in 1:nrow(genes.gwas.pfams.go)){
  gene.scf <- genes.gwas.pfams.go[gene,2]
  gene.start <- genes.gwas.pfams.go[gene,3]
  gene.end <- genes.gwas.pfams.go[gene,4]
  gene.id <- genes.gwas.pfams.go[gene,1]
  snp <- filter(pos.sel.noramB.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.noramA.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.eurasia.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.euro.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- 0
    snp.list <- c(snp_gene, snp.list)
    pvalue <- NA
    p.list <- c(pvalue, p.list)
    gene.id.list <- c(gene.id, gene.id.list)
  } else {
    print("snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- length(snp$pos)
    snp.list <- c(snp_gene, snp.list)
    p.list <- c(paste(snp$BH, collapse = ' '), p.list)
    gene.id.list <- c(gene.id, gene.id.list)
  }
}

length(scf.list)
length(gene.list)
length(end.list)
length(snp.list)
length(p.list)
length(gene.id.list)
snp.list
p.list

noramB <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
noramA <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
eurasia <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
euro <- cbind(scf.list, gene.list, end.list, snp.list, p.list, gene.id.list)
str(euro)

noramA <- as.data.frame(noramA)
colnames(noramA) <- c("scf.noramA","start.noramA", "end.noramA", "noramA_sign_pos_sel_snps", "noramA_p-value")
noramB <- as.data.frame(noramB)
colnames(noramB) <- c("scf.noramB","start.noramB", "end.noramB", "noramB_sign_pos_sel_snps", "noramB_p-value")
eurasia <- as.data.frame(eurasia)
colnames(eurasia) <- c("scf.eurasia","start.eurasia", "end.eurasia", "eurasia_sign_pos_sel_snps", "eurasia_p-value")
euro <- as.data.frame(euro)
colnames(euro) <- c("scf.euro","start.euro", "end.euro", "euro_sign_pos_sel_snps", "euro_p-value", "genes.in.snps.geneid")

all_pops <- cbind(noramA, noramB, eurasia, euro)
str(all_pops)
head(all_pops)

str(genes.gwas.pfams.go)

table(genes.gwas.pfams.go$genes.in.snps.geneid %in% all_pops$genes.in.snps.geneid) #407
table(all_pops$genes.in.snps.geneid %in% genes.gwas.pfams.go$genes.in.snps.geneid) #407

my_ids <- as.vector(genes.gwas.pfams.go$genes.in.snps.geneid)
my_ids

all_pops.ordered <- all_pops %>% arrange(factor(genes.in.snps.geneid, levels = my_ids))
all_pops.ordered <-  all_pops[match(all_pops$genes.in.snps.geneid, my_ids), ]
all_pops$genes.in.snps.geneid

str(all_pops.ordered)
head(all_pops.ordered)

setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Incompatibility_cross_analysis")
write_csv2(all_pops.ordered, "gwas_all_pops.ordered.csv")

###################
#                 #
#     SCOPA       #
#                 #
###################
genes.scopa.pfams.go <- read.csv2("scopa_pfam_go.csv")
str(genes.scopa.pfams.go)

gene.list <- c()
end.list <- c()
scf.list <- c()
snp.list <- c()
p.list <- c()
gene.id.list <- c()
for(gene in 1:nrow(genes.scopa.pfams.go)){
  gene.scf <- genes.scopa.pfams.go[gene,2]
  gene.start <- genes.scopa.pfams.go[gene,3]
  gene.end <- genes.scopa.pfams.go[gene,4]
  gene.id <- genes.scopa.pfams.go[gene,1]
  #snp <- filter(pos.sel.noramB.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.noramA.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  #snp <- filter(pos.sel.eurasia.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  snp <- filter(pos.sel.euro.filtered.BH.sign, rs_id==gene.scf & pos >= gene.start & pos <= gene.end)
  if (dim(snp)[1] == 0) {
    print("no snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- 0
    snp.list <- c(snp_gene, snp.list)
    pvalue <- NA
    p.list <- c(pvalue, p.list)
    gene.id.list <- c(gene.id, gene.id.list)
  } else {
    print("snps in gene")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    end.list <- c(gene.end, end.list)
    snp_gene <- length(snp$pos)
    snp.list <- c(snp_gene, snp.list)
    p.list <- c(paste(snp$BH, collapse = ' '), p.list)
    gene.id.list <- c(gene.id, gene.id.list)
  }
}

length(scf.list)
length(gene.list)
length(end.list)
length(snp.list)
length(p.list)
length(gene.id.list)
snp.list
p.list

noramB <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
noramA <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
eurasia <- cbind(scf.list, gene.list, end.list, snp.list, p.list)
euro <- cbind(scf.list, gene.list, end.list, snp.list, p.list, gene.id.list)
str(euro)

noramA <- as.data.frame(noramA)
colnames(noramA) <- c("scf.noramA","start.noramA", "end.noramA", "noramA_sign_pos_sel_snps", "noramA_p-value")
noramB <- as.data.frame(noramB)
colnames(noramB) <- c("scf.noramB","start.noramB", "end.noramB", "noramB_sign_pos_sel_snps", "noramB_p-value")
eurasia <- as.data.frame(eurasia)
colnames(eurasia) <- c("scf.eurasia","start.eurasia", "end.eurasia", "eurasia_sign_pos_sel_snps", "eurasia_p-value")
euro <- as.data.frame(euro)
colnames(euro) <- c("scf.euro","start.euro", "end.euro", "euro_sign_pos_sel_snps", "euro_p-value", "genes.in.snps.geneid")

all_pops <- cbind(noramA, noramB, eurasia, euro)
str(all_pops)
head(all_pops)
str(genes.scopa.pfams.go)

table(genes.scopa.pfams.go$genes.in.snps.geneid %in% all_pops$genes.in.snps.geneid) #1580
table(all_pops$genes.in.snps.geneid %in% genes.scopa.pfams.go$genes.in.snps.geneid) #1580

my_ids <- as.vector(genes.scopa.pfams.go$genes.in.snps.geneid)
my_ids

all_pops.ordered <- all_pops %>% arrange(factor(genes.in.snps.geneid, levels = my_ids))
#all_pops.ordered <-  all_pops[match(all_pops$genes.in.snps.geneid, my_ids), ]
all_pops$genes.in.snps.geneid

str(all_pops.ordered)
head(all_pops.ordered)

setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Incompatibility_cross_analysis")
write_csv2(all_pops.ordered, "scopa_all_pops.ordered.csv")



