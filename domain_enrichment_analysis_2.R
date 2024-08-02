#R script to 
#1) extract pfam of genes within a given set of SNPs
#2) use fischers exact test to see if the pfam occurence in a given subset (SNPs or genomic windows)
#are significantly different from that of the whole genome

setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/Incompatibility_cross_analysis")
library(tidyverse)
library(stringr)
library(stringi)
library(stats)

###read in annotation:
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

###############################
#                             #
# Genes in sign. GWAS SNPs    # 
#                             #
###############################

###read in SNPs in GWAS analysis that are significant after bonferroni correction
gwas.all.bonferroni <- read.table("GWAS/output_amm/res_amm_gwas_all_inferred_sorted_bonferroni.txt")
str(gwas.all.bonferroni)

gwas.all.bonferroni.out <- gwas.all.bonferroni[,c(2,4)]
str(gwas.all.bonferroni.out)
gwas.all.bonferroni.out$end <- gwas.all.bonferroni.out$Pos + 1
gwas.all.bonferroni.out <- arrange(gwas.all.bonferroni.out, gwas.all.bonferroni.out$Chr)
str(gwas.all.bonferroni.out)
write_delim(gwas.all.bonferroni.out, file="gwas_bonferroni_snp_extn.tab", delim ="\t")

###get genes that contain significant gwas snps:
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
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
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}

length(scf.list) 
length(gene.list) #407 genes that contain significant GWAS snps
length(num.snps)
length(gene.info)
a <- cbind(scf.list, gene.list, num.snps, gene.info)
head(a)
genes.in.snps <- as.data.frame(a)
str(genes.in.snps)
b <- cbind(scf.list, gene.list, end.list)
b <- as.data.frame(b)
write_delim(b, "gwas_genes_pos.txt", delim = "\t")
#now also parse out each individual pfam in these 407 genes:
genes.in.snps.pfams <- str_extract_all(genes.in.snps$gene.info,"Pfam:PF[0-9]{5}") #extracts pfams
length(genes.in.snps.pfams) #matches number of genes in snps (407)
lengths(genes.in.snps.pfams) #number of pfams in each gene
str(genes.in.snps.pfams)
max(lengths(genes.in.snps.pfams)) #11
genes.in.snps.pfams.df <- as.data.frame(t(stri_list2matrix(genes.in.snps.pfams))) #convert to dataframe, each pfam for each gene is a column
str(genes.in.snps.pfams.df)
genes.in.snps.pfam.split <- cbind(genes.in.snps, genes.in.snps.pfams.df)
str(genes.in.snps.pfam.split)
dim(genes.in.snps.pfam.split)

###Do fischers exact test
#two classes: 1)pfam all genes and 2)pfam genes with snps
#occurences of each pfam in each class

#make a count of pfams:
#1) all genes (with pfam):
n_distinct(genome.genes.pfam.split$V1) #1675
sapply(genome.genes.pfam.split, function(x) n_distinct(x))
genome.pfams <- c(genome.genes.pfam.split$V1,
                      genome.genes.pfam.split$V2,
                      genome.genes.pfam.split$V3,
                      genome.genes.pfam.split$V4,
                      genome.genes.pfam.split$V5,
                      genome.genes.pfam.split$V6,
                      genome.genes.pfam.split$V7,
                      genome.genes.pfam.split$V8,
                      genome.genes.pfam.split$V9,
                      genome.genes.pfam.split$V10,
                      genome.genes.pfam.split$V11,
                      genome.genes.pfam.split$V12,
                      genome.genes.pfam.split$V13)
str(genome.pfams)
length(genome.pfams)
genome.pfams.counts <- table(genome.pfams)
genome.pfams.counts <- as.data.frame(genome.pfams.counts)
colnames(genome.pfams.counts) <- c("pfam","freq_all")
str(genome.pfams.counts) #2680 unique pfams in total 4422 genes

#2) gwas genes with pfam:
selected.genes.pfams <- c(genes.in.snps.pfam.split$V1,
                         genes.in.snps.pfam.split$V2,
                         genes.in.snps.pfam.split$V3,
                         genes.in.snps.pfam.split$V4,
                         genes.in.snps.pfam.split$V5,
                         genes.in.snps.pfam.split$V6,
                         genes.in.snps.pfam.split$V7,
                         genes.in.snps.pfam.split$V8,
                         genes.in.snps.pfam.split$V9,
                         genes.in.snps.pfam.split$V10,
                         genes.in.snps.pfam.split$V11)
length(selected.genes.pfams) 
selected.genes.pfams.counts <- table(selected.genes.pfams)
selected.genes.pfams.counts <- as.data.frame(selected.genes.pfams.counts)
colnames(selected.genes.pfams.counts) <- c("pfam", "freq_selected")
str(selected.genes.pfams.counts) #651 pfams in 407 gwas genes
head(selected.genes.pfams.counts)

#combine 1) and 2):
combined.table.all <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam", all.x = TRUE)
dim(combined.table.all) #2680 pfams
head(combined.table.all,20)
#replace NA with 0:
combined.table.all[is.na(combined.table.all)] = 0
head(combined.table.all,20)
fisher.table <- combined.table.all[,2:3]
rownames(fisher.table) <- combined.table.all$pfam
head(fisher.table)
fisher.table <- t(fisher.table)
dim(fisher.table)
test <- fisher.test(fisher.table, simulate.p.value=TRUE) #must use simulate to not exceed memory limitation
test #p-value=1

#Actually classes should be mutually exclusive so they can add up to total, so what you want is:
#1) all pfam in snp genes VERSUS all pfam not in snp genes
#2) all pfam outside snp genes VERSUS all pfam not outside snp genes
combined.table.snps <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam")
dim(combined.table.snps) #651: i.e. omits pfams not present in snp genes
head(combined.table.snps)
not.snps.table <- subset(combined.table.all, !(combined.table.all$pfam %in% combined.table.snps$pfam))
dim(not.snps.table) #2029 pfams not in snps genes, i.e. only outside snps
#pfam not found outside snp genes:
only.snps.table <- filter(combined.table.snps, freq_all == freq_selected)
dim(only.snps.table) #245 pfams only in snp genes, i.e. 245 pfams absent outside snps
table(genome.genes.pfam.split$V1 %in% only.snps.table$pfam) #81, does this equal to the number of genes?
table(genome.genes.pfam.split$V2 %in% only.snps.table$pfam) #62
table(genome.genes.pfam.split$V3 %in% only.snps.table$pfam) #50
table(genome.genes.pfam.split$V4 %in% only.snps.table$pfam) #23
table(genome.genes.pfam.split$V5 %in% only.snps.table$pfam) #20
table(genome.genes.pfam.split$V6 %in% only.snps.table$pfam) #10
table(genome.genes.pfam.split$V7 %in% only.snps.table$pfam) #7
table(genome.genes.pfam.split$V8 %in% only.snps.table$pfam) #3
table(genome.genes.pfam.split$V9 %in% only.snps.table$pfam) #2
table(genome.genes.pfam.split$V10 %in% only.snps.table$pfam) #0
81 + 62 + 50 + 23 + 20 + 20 + 10 + 7 + 3 + 2 #278 pfams, but some duplicates in here (one gene can have multiple pfams), therefore > 245
table(not.snps.table$pfam %in% only.snps.table$pfam) #appears to not be overlapping
table(only.snps.table$pfam  %in% not.snps.table$pfam) #appears to not be overlapping
pfams.present <- c(651,2680-245)
pfams.absent <- c(645,245)
fisher.combined.table <- cbind(pfams.present, pfams.absent)
rownames(fisher.combined.table) <- c("within_snp","outside_snp")
fisher.combined.table
test <- fisher.test(fisher.combined.table) 
test
test$p.value #3.431034e-174
#The above is only a presence absence of each pfam and doesn´t take into account counts of each pfam

#Try instead to test pfam by pfam with counts, and then correct for multiple testing:
head(combined.table.all,3)
#counts pfam present outside snp genes: freq_all - freq_selected
#counts pfam present inside snp genes: freq_selected
#counts pfams absent outside snp genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside snp genes: 4422 - freq_selected

pvalues <- c() #loop through each pfam and do an F-test for each:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 407-freq_selected
  absent_outside <- 4422-absent_within_snp-present_outside_snp-present_within_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction
filter(sorted.final.table, freq_outside==0) #get pfams that only within genes with GWAS snps

#as above, but make test one-sided:
pvalues <- c() #loop through each pfam and do an F-test for each:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 407-freq_selected
  absent_outside <- 4422-absent_within_snp-present_outside_snp-present_within_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  test <- fisher.test(cont_table, alternative = "greater")
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction, doesn´t seem to make a difference
filter(sorted.final.table, freq_outside==0) #get pfams that only within genes with GWAS snps

#test instead subset of snp genes versus all genes (but does this violate the assumption of exclusivity of categories?)
#or perhaps here the test is rather testing for the significance of the overlap?
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 407-freq_selected
  absent_outside <- 4422-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #1 significant after correction, pfam absent in snp genes:
#Pfam:PF13515 (Fusaric acid resistance protein-like): https://www.ebi.ac.uk/interpro/entry/pfam/PF13515/ 

#go even further by replacing "outside_snp_gene" with counts of total:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_all <- freq_all
  absent_within_snp <- 407-freq_selected
  absent_outside <- 4422-freq_all
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","all_genes")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #23 significant after correction
filter(sorted.final.table, BH < 0.05 & freq_outside==0) #1 pfam occuring only within gwas genes:
#Pfam:PF13185 (GAF domain): https://www.ebi.ac.uk/interpro/entry/pfam/PF13185/
filter(sorted.final.table, BH < 0.05 & freq_selected > 0) #7 significant pfams after correction that are present in gwas genes:
#Pfam:PF14027 (Questin oxidase-like): https://www.ebi.ac.uk/interpro/entry/pfam/PF14027/
#Pfam:PF13185 (GAF domain): https://www.ebi.ac.uk/interpro/entry/pfam/PF13185/
#Pfam:PF14033 (Protein of unknown function in bacteria and fungi): https://www.ebi.ac.uk/interpro/entry/pfam/PF14033/
#Pfam:PF13923 (Zinc finger): https://www.ebi.ac.uk/interpro/entry/pfam/PF13923/
#Pfam:PF13193 (AMP-binding enzyme C-terminal domain): https://www.ebi.ac.uk/interpro/entry/pfam/PF13193/
#Pfam:PF13774 (Regulated-SNARE-like domain, protein transport from ER to plasma membrane): https://www.ebi.ac.uk/interpro/entry/pfam/PF13774/ 
#Pfam:PF00083 (Glucose transporter): https://www.ebi.ac.uk/interpro/entry/pfam/PF00083/

#exclude pfams that are not present in snp genes? (will make BH correction softer as fewer tests are conducted)
pfam.snp.genes <- filter(combined.table.all, freq_selected > 0)
dim(pfam.snp.genes) #651 pfams present in 407 snp genes
head(pfam.snp.genes, 5)
a <- filter(pfam.snp.genes, freq_selected == freq_all)
dim(a) #245 pfams present only in snp genes, i.e. 651-245=406 pfams shared
a <- filter(pfam.snp.genes, freq_selected < freq_all)
dim(a) #406
str(genome.genes.pfam.split)
b <- filter(genome.genes.pfam.split, genome.genes.pfam.split$V1 %in% a$pfam)
dim(b)#2211 genes with shared pfam occuring in total, i.e. 2211-407=1804 genes outside snps genes carrying shared pfams?
#in addition 251 pfams uniquely present in 81 snp genes (ref line 162-163)
pvalues <- c()
for(i in 1:nrow(pfam.snp.genes)){
  freq_all <- pfam.snp.genes[i,2]
  freq_selected <- pfam.snp.genes[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 407-freq_selected
  absent_outside <- 2211-absent_within_snp-present_outside_snp-present_within_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(pfam.snp.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #still nothing, but Pfam:PF01213 gets close to sign.
#Pfam:PF01213 (cyclase-associated protein family, actin-binding, involved in fruit body formation in slime molds,
# in yeast involved in cyclase activation and vesicle trafficking and endocytosis)
# https://www.ebi.ac.uk/interpro/entry/pfam/PF01213/
filter(sorted.final.table, freq_outside==0)

#try to not count genes, but instead use pfams:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 651-freq_selected
  absent_outside <- 2680-absent_within_snp-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05)

#Pfam:PF14027 (Questin oxidase-like): https://www.ebi.ac.uk/interpro/entry/pfam/PF14027/
#Pfam:PF14033 (functionally uncharacterised): https://www.ebi.ac.uk/interpro/entry/pfam/PF14033/

#or actually what should be more correct is to use counts of pfams:
sum(combined.table.all$freq_all) #7167 total counts of pfams
sum(combined.table.all$freq_selected) #903 total counts of pfams in gwas genes
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 903-freq_selected
  absent_outside <- 7167-absent_within_snp-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0

#only use pfams in gwas genes:
str(pfam.snp.genes)
sum(pfam.snp.genes$freq_all) #3536 total counts of pfams
sum(pfam.snp.genes$freq_selected) #903 total counts of pfams in gwas genes
pvalues <- c()
for(i in 1:nrow(pfam.snp.genes)){
  freq_all <- pfam.snp.genes[i,2]
  freq_selected <- pfam.snp.genes[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 903-freq_selected
  absent_outside <- 3536-absent_within_snp-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(pfam.snp.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05)
#2 Pfam:PF01213  (Adenylate cyclase associated (CAP) N termin; actin binding proteins): https://www.ebi.ac.uk/interpro/entry/pfam/PF01213/
#3 Pfam:PF13774   (Longin; transport of proteins from ER to plasma membrane): https://www.ebi.ac.uk/interpro/entry/pfam/PF13774/


#try to quantiy the overlap between pfams in gwas snp genes and all genes:
#1. pfams that are provided in both lists: 651
#2. pfams that are in gwas snps but not outside: 245
#3. pfams that are outside and not in gwas snps: 2029
#4. pfams that are absent from both (?): 0  -> so this is probably not the right way
presents <- c(651,245)
absents <- c(2029, 0)
overlap.F <- rbind(presents, absents) 
overlap.F
fisher.test(overlap.F)
#p-value < 2.2e-16

#######################################
#                                     #
# genes in SCOPA snps                 #
#                                     #  
#######################################

#read in scopa snps
scopa_noram <- read.table("SCOPA_crosses_NorAm_num_scf_gen_result_simple_sorted.txt", header=T)
str(scopa_noram)
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

###get genes that contain significant scopa snps:
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
end.list <- c()
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
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    gene.info <- c(gene.anno, gene.info)
    end.list <- c(gene.end, end.list)
  }
}

length(scf.list) 
length(gene.list) #1580 genes that contain significant GWAS snps
length(num.snps)
length(gene.info)
a <- cbind(scf.list, gene.list, num.snps, gene.info)
head(a)
genes.in.snps <- as.data.frame(a)
str(genes.in.snps)

b <- cbind(scf.list, gene.list, end.list)
b <- as.data.frame(b)
write_delim(b, "scopa_genes_pos.txt", delim = "\t")

#now also parse out each individual pfam in these 1580 genes:
genes.in.snps.pfams <- str_extract_all(genes.in.snps$gene.info,"Pfam:PF[0-9]{5}") #extracts pfams
length(genes.in.snps.pfams) #matches number of genes in snps (1580)
lengths(genes.in.snps.pfams) #number of pfams in each gene
genes.in.snps.pfams.df <- as.data.frame(t(stri_list2matrix(genes.in.snps.pfams))) #convert to dataframe, each pfam for each gene is a column
str(genes.in.snps.pfams.df)
genes.in.snps.pfam.split <- cbind(genes.in.snps, genes.in.snps.pfams.df)
str(genes.in.snps.pfam.split)
dim(genes.in.snps.pfam.split)

###Do fischers exact test
#two classes: 1)pfam all genes and 2)pfam genes with snps
#occurences of each pfam in each class

#make a count of pfams:
#1) all genes (with pfam):
n_distinct(genome.genes.pfam.split$V1) 
sapply(genome.genes.pfam.split, function(x) n_distinct(x))
genome.pfams <- c(genome.genes.pfam.split$V1,
                  genome.genes.pfam.split$V2,
                  genome.genes.pfam.split$V3,
                  genome.genes.pfam.split$V4,
                  genome.genes.pfam.split$V5,
                  genome.genes.pfam.split$V6,
                  genome.genes.pfam.split$V7,
                  genome.genes.pfam.split$V8,
                  genome.genes.pfam.split$V9,
                  genome.genes.pfam.split$V10,
                  genome.genes.pfam.split$V11,
                  genome.genes.pfam.split$V12,
                  genome.genes.pfam.split$V13)
str(genome.pfams)
length(genome.pfams)
genome.pfams.counts <- table(genome.pfams)
genome.pfams.counts <- as.data.frame(genome.pfams.counts)
colnames(genome.pfams.counts) <- c("pfam","freq_all")
str(genome.pfams.counts) #2680 unique pfams in total 4422 genes

#2) scopa genes with pfam:
selected.genes.pfams <- c(genes.in.snps.pfam.split$V1,
                          genes.in.snps.pfam.split$V2,
                          genes.in.snps.pfam.split$V3,
                          genes.in.snps.pfam.split$V4,
                          genes.in.snps.pfam.split$V5,
                          genes.in.snps.pfam.split$V6,
                          genes.in.snps.pfam.split$V7,
                          genes.in.snps.pfam.split$V8,
                          genes.in.snps.pfam.split$V9,
                          genes.in.snps.pfam.split$V10,
                          genes.in.snps.pfam.split$V11)
length(selected.genes.pfams) 
selected.genes.pfams.counts <- table(selected.genes.pfams)
selected.genes.pfams.counts <- as.data.frame(selected.genes.pfams.counts)
colnames(selected.genes.pfams.counts) <- c("pfam", "freq_selected")
str(selected.genes.pfams.counts) #1692 pfams in 1580 scopa genes
head(selected.genes.pfams.counts)

#combine 1) and 2):
combined.table.all <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam", all.x = TRUE)
dim(combined.table.all) #2680 pfams
head(combined.table.all,20)
combined.table.all[is.na(combined.table.all)] = 0 #replace NA with 0
head(combined.table.all,20)
fisher.table <- combined.table.all[,2:3]
rownames(fisher.table) <- combined.table.all$pfam
head(fisher.table)
fisher.table <- t(fisher.table)
dim(fisher.table)
test <- fisher.test(fisher.table, simulate.p.value=TRUE) #must use simulate to not exceed memory limitation
test #p-value=1

#Actually classes should be mutually exclusive so they can add up to total, so what you want is:
#1) all pfam in snp genes VERSUS all pfam not in snp genes
#2) all pfam outside snp genes VERSUS all pfam not outside snp genes
combined.table.snps <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam")
dim(combined.table.snps) #1692: i.e. omits pfams not present in snp genes
head(combined.table.snps)
not.snps.table <- subset(combined.table.all, !(combined.table.all$pfam %in% combined.table.snps$pfam))
dim(not.snps.table) #988 pfams not in snps genes, i.e. only outside snps
#pfam not found outside snp genes:
only.snps.table <- filter(combined.table.snps, freq_all == freq_selected)
dim(only.snps.table) #969 pfams only in snp genes, i.e. 969 pfams absent outside snps
table(genome.genes.pfam.split$V1 %in% only.snps.table$pfam) #477, does this equal to the number of genes?
table(genome.genes.pfam.split$V2 %in% only.snps.table$pfam) #303
table(genome.genes.pfam.split$V3 %in% only.snps.table$pfam) #200
table(genome.genes.pfam.split$V4 %in% only.snps.table$pfam) #100
table(genome.genes.pfam.split$V5 %in% only.snps.table$pfam) #54
table(genome.genes.pfam.split$V6 %in% only.snps.table$pfam) #23
table(genome.genes.pfam.split$V7 %in% only.snps.table$pfam) #14
table(genome.genes.pfam.split$V8 %in% only.snps.table$pfam) #9
table(genome.genes.pfam.split$V9 %in% only.snps.table$pfam) #4
table(genome.genes.pfam.split$V10 %in% only.snps.table$pfam) #0
477 + 303 + 200 + 100 + 54 + 23 + 14 + 9 + 4 #1184 pfams, but some duplicates in here, therefore > 969
table(not.snps.table$pfam %in% only.snps.table$pfam) #appears to not be overlapping
table(only.snps.table$pfam  %in% not.snps.table$pfam) #appears to not be overlapping

#Try to test pfam by pfam with counts, and then correct for multiple testing:
head(combined.table.all,3)
#counts pfam present outside snp genes: freq_all - freq_selected
#counts pfam present inside snp genes: freq_selected
#counts pfams absent outside snp genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside snp genes: 4422 - freq_selected

pvalues <- c() #loop through each pfam and do a Fischers exact test:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 1580-freq_selected
  absent_outside <- 4422-1580-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #none significant after correction

#test instead subset of snp genes versus all genes (but does this violate the assumption of exclusivity of categories?)
#replace "outside_snp_genes" where pfam is absent with all genes (4422)
# but shouldn´t "outside_snp_genes" where pfam is present also be replaced with all genes?

pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 1580-freq_selected
  absent_outside <- 4422-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #7 pfams significant
filter(sorted.final.table, BH < 0.05 & freq_selected > 0) #6 pfams significant that are present in scopa snp genes
#Pfam:PF14027 (Questin oxidase-like): https://www.ebi.ac.uk/interpro/entry/pfam/PF14027/
#Pfam:PF04675 (DNA ligase N terminus): https://www.ebi.ac.uk/interpro/entry/pfam/PF04675/
#Pfam:PF13515 (Fusaric acid resistance protein-like): https://www.ebi.ac.uk/interpro/entry/pfam/PF13515/
#Pfam:PF13185
#Pfam:PF13696
#Pfam:PF12449

#go all the way by also replacing "outside_snp_gene" for pfam present with counts of total:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_all <- freq_all
  absent_within_snp <- 1580-freq_selected
  absent_all <- 4422-freq_all
  genes_pfam_present <- c(present_within_snp,present_all)
  genes_pfam_absent <- c(absent_within_snp,absent_all)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","all_genes")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction

#exclude pfams that are not present in snp genes? (will make BH correction softer as fewer tests are conducted)
#but this is probably "unfair"...
pfam.snp.genes <- filter(combined.table.all, freq_selected > 0)
dim(pfam.snp.genes) #1692 pfams present in 1580 snp genes
head(pfam.snp.genes, 5)
a <- filter(pfam.snp.genes, freq_selected == freq_all)
dim(a) #969 pfams present only in snp genes
a <- filter(pfam.snp.genes, freq_selected < freq_all)
dim(a) #723 pfams shared snp genes and genes outside snps
str(genome.genes.pfam.split)
b <- filter(genome.genes.pfam.split, genome.genes.pfam.split$V1 %in% a$pfam)
dim(b)#3020 genes with shared pfam occuring in total, i.e. 3020-1580=1440 genes outside snps genes carrying shared pfams?
pvalues <- c()
for(i in 1:nrow(pfam.snp.genes)){
  freq_all <- pfam.snp.genes[i,2]
  freq_selected <- pfam.snp.genes[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 1580-freq_selected
  absent_outside <- 1440-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(pfam.snp.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #8 pfams significant after correction:
#Pfam:PF13949 (ALIX V-shaped domain binding to HIV): https://www.ebi.ac.uk/interpro/entry/pfam/PF13949/
#Pfam:PF00380 (Ribosomal protein S9/S16): https://www.ebi.ac.uk/interpro/entry/pfam/PF00380/ 
#Pfam:PF01769 (Divalent cation transporter): https://www.ebi.ac.uk/interpro/entry/pfam/PF01769/ 
#Pfam:PF11711 (Inner membrane protein import complex subunit Tim54 -> import to mitochondria): https://www.ebi.ac.uk/interpro/entry/pfam/PF11711/
#Pfam:PF00204 (DNA gyrase B): https://www.ebi.ac.uk/interpro/entry/pfam/PF00204/
#Pfam:PF13877 (Potential Monad-binding region of RPAP3 -> regulate APOPTOSIS, contain TPR-repeats towards the terminus: https://www.ebi.ac.uk/interpro/entry/pfam/PF13877/
#Pfam:PF13907 (Domain of unknown function, found at the C-terminus of chromodomain-helicase-DNA-binding proteins): https://www.ebi.ac.uk/interpro/entry/pfam/PF13907/
#Pfam:PF13862 (BCCIP; In fungi involved in nuclear export, actin cytoskeleton organisation and vesicular transport, homologue in human may promote cell cycle arrest): https://www.ebi.ac.uk/interpro/entry/pfam/PF13862/

#try to do with pfam counts instead of genes:
sum(combined.table.all$freq_all) #7167 total counts of pfams
sum(combined.table.all$freq_selected) #3177 total counts of pfams in scopa genes
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 3177-freq_selected
  absent_outside <- 7167-3177-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05)
#Pfam:PF14033 (functionally uncharacterised): https://www.ebi.ac.uk/interpro/entry/pfam/PF14033/
#Pfam:PF00350 (dynamin, endocytosis in the eukaryotic cell): https://www.ebi.ac.uk/interpro/entry/pfam/PF00350/

#counts of pfam that occur in scopa genes:
  #note here that probably it´s not fair to infer anything about those genes that don´t occur outside scopa genes 
  #that is pfams that don´t occur within scopa genes
as you´ve eliminated their counterpart, t
str(pfam.snp.genes)
sum(pfam.snp.genes$freq_all) #5795 total counts of pfams
sum(pfam.snp.genes$freq_selected) #3177 total counts of pfams in gwas genes
pvalues <- c()
for(i in 1:nrow(pfam.snp.genes)){
  freq_all <- pfam.snp.genes[i,2]
  freq_selected <- pfam.snp.genes[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 3177-freq_selected
  absent_outside <- 5795-absent_within_snp-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(t(cont_table), alternative="greater")
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(pfam.snp.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0

#Pfam:PF00380 (Ribosomal protein S9/S16; small ribosomal subunit) https://www.ebi.ac.uk/interpro/entry/pfam/PF00380/
#Pfam:PF01769 (Divalent cation transporter, Mg+) https://www.ebi.ac.uk/interpro/entry/pfam/PF01769/
#Pfam:PF11711 (Inner membrane protein import complex subunit Tim54, import into mitochondrion) https://www.ebi.ac.uk/interpro/entry/pfam/PF11711/ 
#Pfam:PF00204 (DNA gyrase B, ribosomal S5 domain 2-like fold) https://www.ebi.ac.uk/interpro/entry/pfam/PF00204/
#Pfam:PF13877 (Potential Monad-binding region of RPAP3, involved in regulating apoptosis, TPR-repeats towards the N_terminus) https://www.ebi.ac.uk/interpro/entry/pfam/PF13877/
#Pfam:PF13862 (BRCA2 and CDKN1A-interacting protein; nuclear export, actin cytoskeleton organisation and vesicular transport,Its homologue in human and mouse, BCCIP, may promote cell cycle arrest) https://www.ebi.ac.uk/interpro/entry/pfam/PF13862/
#Pfam:PF13907 (Domain of unknown function; C-terminus of chromodomain-helicase-DNA-binding proteins): https://www.ebi.ac.uk/interpro/entry/pfam/PF13907/
#Pfam:PF06831 (Formamidopyrimidine-DNA glycosylase H2TH domain; DNA repair enzyme) https://www.ebi.ac.uk/interpro/entry/pfam/PF06831/
#Pfam:PF12689 (Acid Phosphatase) https://www.ebi.ac.uk/interpro/entry/pfam/PF12689/


#######################################
#                                     #
# genes in deltafst outlier windows   #
#                                     #  
#######################################

delta_fst.all <- read.table("Delta_df_all_outliers.txt", header=T)
str(delta_fst.all)

delta_fst.noram <- read.table("Delta_df_NoramA_NoramB.txt", header=T)
str(delta_fst.noram)

wd40.delta.fst.noram <-left_join(pfam004000.gene, delta_fst.noram, by = c("V1"="Scaffold")) %>% filter(V4 >= Start, V5 <= End)
str(wd40.delta.fst.noram)

pfam05729.delta.fst.all <- filter(delta_fst.all, Scaffold=="Scaffold06" & Start > 1439976 & End < 1439976 + 11000)

gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
for(gene in 1:nrow(genome.genes.pfam.split)) {
  gene.scf <- genome.genes.pfam[gene,1]
  gene.start <- genome.genes.pfam[gene,4]
  gene.end <- genome.genes.pfam[gene,5]
  gene.anno <- genome.genes.pfam[gene,9]
  snp <- filter(delta_fst.all, Scaffold==gene.scf & gene.end >= Start & gene.start <= End)
  if (dim(snp)[1] == 0) {
    print("no overlap")  
  } else {
    print("overlap")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    gene.info <- c(gene.anno, gene.info)
  }
}

length(scf.list) 
length(gene.list) #1260 genes with pfam that overlap with deltafst windows
length(num.snps)
length(gene.info)
a <- cbind(scf.list, gene.list, num.snps, gene.info)
head(a)
genes.in.delta <- as.data.frame(a)
str(genes.in.delta)
#now also parse out each individual pfam in these 407 genes:
genes.in.delta.pfams <- str_extract_all(genes.in.delta$gene.info,"Pfam:PF[0-9]{5}") #extracts pfams
length(genes.in.delta.pfams) #matches number of genes in snps (1260)
lengths(genes.in.delta.pfams) #number of pfams in each gene
max(lengths(genes.in.delta.pfams)) #13
genes.in.delta.pfams.df <- as.data.frame(t(stri_list2matrix(genes.in.delta.pfams))) #convert to dataframe, each pfam for each gene is a column
str(genes.in.delta.pfams.df)
genes.in.delta.pfam.split <- cbind(genes.in.delta, genes.in.delta.pfams.df)
str(genes.in.delta.pfam.split)
dim(genes.in.delta.pfam.split)

###Do fischers exact test
#two classes: 1)pfam all genes and 2)pfam genes with snps
#occurences of each pfam in each class

#make a count of pfams:
#1) all genes (with pfam):
n_distinct(genome.genes.pfam.split$V1) 
sapply(genome.genes.pfam.split, function(x) n_distinct(x))
genome.pfams <- c(genome.genes.pfam.split$V1,
                  genome.genes.pfam.split$V2,
                  genome.genes.pfam.split$V3,
                  genome.genes.pfam.split$V4,
                  genome.genes.pfam.split$V5,
                  genome.genes.pfam.split$V6,
                  genome.genes.pfam.split$V7,
                  genome.genes.pfam.split$V8,
                  genome.genes.pfam.split$V9,
                  genome.genes.pfam.split$V10,
                  genome.genes.pfam.split$V11,
                  genome.genes.pfam.split$V12,
                  genome.genes.pfam.split$V13)
str(genome.pfams)
length(genome.pfams)
genome.pfams.counts <- table(genome.pfams)
genome.pfams.counts <- as.data.frame(genome.pfams.counts)
colnames(genome.pfams.counts) <- c("pfam","freq_all")
str(genome.pfams.counts) #2680 unique pfams in total 4422 genes

#2) deltafst genes with pfam:
selected.genes.pfams <- c(genes.in.delta.pfam.split$V1,
                          genes.in.delta.pfam.split$V2,
                          genes.in.delta.pfam.split$V3,
                          genes.in.delta.pfam.split$V4,
                          genes.in.delta.pfam.split$V5,
                          genes.in.delta.pfam.split$V6,
                          genes.in.delta.pfam.split$V7,
                          genes.in.delta.pfam.split$V8,
                          genes.in.delta.pfam.split$V9,
                          genes.in.delta.pfam.split$V10,
                          genes.in.delta.pfam.split$V11)
length(selected.genes.pfams) 
selected.genes.pfams.counts <- table(selected.genes.pfams)
selected.genes.pfams.counts <- as.data.frame(selected.genes.pfams.counts)
colnames(selected.genes.pfams.counts) <- c("pfam", "freq_selected")
str(selected.genes.pfams.counts) #1161 pfams in 1260 deltafst genes
head(selected.genes.pfams.counts)

#combine 1) and 2):
combined.table.all <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam", all.x = TRUE)
dim(combined.table.all) #2680 pfams
head(combined.table.all,20)
combined.table.all[is.na(combined.table.all)] = 0 #replace NA with 0
head(combined.table.all,20)
fisher.table <- combined.table.all[,2:3]
rownames(fisher.table) <- combined.table.all$pfam
head(fisher.table)
fisher.table <- t(fisher.table)
dim(fisher.table)
test <- fisher.test(fisher.table, simulate.p.value=TRUE) #must use simulate to not exceed memory limitation
test #p-value=1

#Actually classes should be mutually exclusive so they can add up to total, so what you want is:
#1) all pfam in snp genes VERSUS all pfam not in snp genes
#2) all pfam outside snp genes VERSUS all pfam not outside snp genes
combined.table.delta <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam")
dim(combined.table.delta) #1161: i.e. omits pfams not present in deltafst genes
head(combined.table.delta)
not.delta.table <- subset(combined.table.all, !(combined.table.all$pfam %in% combined.table.delta$pfam))
dim(not.delta.table) #1519 pfams not in snps genes, i.e. only outside snps
#pfam not found outside snp genes:
only.delta.table <- filter(combined.table.delta, freq_all == freq_selected)
dim(only.delta.table) #488 pfams only in deltafst genes, i.e. 488 pfams absent outside deltafst
table(genome.genes.pfam.split$V1 %in% only.delta.table$pfam) #248, does this equal to the number of genes?
table(genome.genes.pfam.split$V2 %in% only.delta.table$pfam) #129
table(genome.genes.pfam.split$V3 %in% only.delta.table$pfam) #77
table(genome.genes.pfam.split$V4 %in% only.delta.table$pfam) #36
table(genome.genes.pfam.split$V5 %in% only.delta.table$pfam) #22
table(genome.genes.pfam.split$V6 %in% only.delta.table$pfam) #11
table(genome.genes.pfam.split$V7 %in% only.delta.table$pfam) #4
table(genome.genes.pfam.split$V8 %in% only.delta.table$pfam) #4
table(genome.genes.pfam.split$V9 %in% only.delta.table$pfam) #3
table(genome.genes.pfam.split$V10 %in% only.delta.table$pfam) #0

table(not.delta.table$pfam %in% only.delta.table$pfam) #appears to not be overlapping
table(only.delta.table$pfam  %in% not.delta.table$pfam) #appears to not be overlapping

#Try to test pfam by pfam with counts, and then correct for multiple testing:
head(combined.table.all,3)
#counts pfam present outside snp genes: freq_all - freq_selected
#counts pfam present inside snp genes: freq_selected
#counts pfams absent outside snp genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside snp genes: 4422 - freq_selected

pvalues <- c() #loop through each pfam and do a Fischers exact test:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 1260-freq_selected
  absent_outside_delta <- 4422-1260-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_delta_gene","outside_delta_gene")
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #none significant after correction

#test instead subset of snp genes versus all genes (but does this violate the assumption of exclusivity of categories?)
#replace "outside_snp_genes" where pfam is absent with all genes (4422)
# but shouldn´t "outside_snp_genes" where pfam is present also be replaced with all genes?
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 1260-freq_selected
  absent_outside_delta <- 4422-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 pfams significant

#go all the way by also replacing "outside_snp_gene" for pfam present with counts of total:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_all <- freq_all
  absent_within_delta <- 1580-freq_selected
  absent_all <- 4422-freq_all
  genes_pfam_present <- c(present_within_delta, present_all)
  genes_pfam_absent <- c(absent_within_delta, absent_all)
  cont_table <- as.data.frame(cbind(genes_pfam_present, genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","all_genes")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction

#exclude pfams that are not present in snp genes? (will make BH correction softer as fewer tests are conducted)
#but this is probably "unfair"...
pfam.delta.genes <- filter(combined.table.all, freq_selected > 0)
dim(pfam.delta.genes) #1161 pfams present in 1260 delta genes
head(pfam.delta.genes, 5)
a <- filter(pfam.delta.genes, freq_selected == freq_all)
dim(a) #488 pfams present only in snp genes
a <- filter(pfam.delta.genes, freq_selected < freq_all)
dim(a) #673 pfams shared snp genes and genes outside snps
str(genome.genes.pfam.split)
b <- filter(genome.genes.pfam.split, genome.genes.pfam.split$V1 %in% a$pfam)
dim(b)#3005 genes with shared pfam occuring in total, i.e. 3005-1260=1745 genes outside delta genes carrying shared pfams?

pvalues <- c()
for(i in 1:nrow(pfam.delta.genes)){
  freq_all <- pfam.delta.genes[i,2]
  freq_selected <- pfam.delta.genes[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 1260-freq_selected
  absent_outside_delta <- 1745-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_delta_gene","outside_delta_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
pfam.delta.genes$freq_outside <- pfam.delta.genes$freq_all - pfam.delta.genes$freq_selected
final.table <- cbind(pfam.delta.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 pfams significant after correction:

#try with pfam counts instead of gene counts:
sum(combined.table.all$freq_all) #7167 total counts of pfams
sum(combined.table.all$freq_selected) #2174 total counts of pfams in deltafst
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 2174-freq_selected
  absent_outside <- 7167-2174-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table, alternative="less")
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0

#only use pfams in gwas genes:
str(pfam.delta.genes)
sum(pfam.delta.genes$freq_all) #5105 total counts of pfams
sum(pfam.delta.genes$freq_selected) #2174 total counts of pfams in gwas genes
pvalues <- c()
for(i in 1:nrow(pfam.delta.genes)){
  freq_all <- pfam.delta.genes[i,2]
  freq_selected <- pfam.delta.genes[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 2174-freq_selected
  absent_outside <- 5105-absent_within_snp-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.delta.genes$freq_outside <- pfam.delta.genes$freq_all - pfam.delta.genes$freq_selected
final.table <- cbind(pfam.delta.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05)


###########################
#                         #  
#   DeltaFst Noram        #
#                         #
###########################

delta_fst.noram <- read.table("Delta_df_NoramA_NoramB.txt", header=T)
str(delta_fst.noram)
gene.list <- c()
scf.list <- c()
snp.list <- c()
num.snps <- c()
gene.info <- c()
for(gene in 1:nrow(genome.genes.pfam.split)) {
  gene.scf <- genome.genes.pfam[gene,1]
  gene.start <- genome.genes.pfam[gene,4]
  gene.end <- genome.genes.pfam[gene,5]
  gene.anno <- genome.genes.pfam[gene,9]
  snp <- filter(delta_fst.noram, Scaffold==gene.scf & gene.end >= Start & gene.start <= End)
  if (dim(snp)[1] == 0) {
    print("no overlap")  
  } else {
    print("overlap")
    scf.list <- c(gene.scf, scf.list)
    gene.list <- c(gene.start, gene.list)
    num.snps <- c(length(snp$Pos), num.snps)
    gene.info <- c(gene.anno, gene.info)
  }
}

length(scf.list) 
length(gene.list) #185 genes with pfam that overlap with deltafst windows
length(num.snps)
length(gene.info)
a <- cbind(scf.list, gene.list, num.snps, gene.info)
head(a)
genes.in.delta <- as.data.frame(a)
str(genes.in.delta)
#now also parse out each individual pfam in these 185 genes:
genes.in.delta.pfams <- str_extract_all(genes.in.delta$gene.info,"Pfam:PF[0-9]{5}") #extracts pfams
length(genes.in.delta.pfams) #matches number of genes in snps (185)
lengths(genes.in.delta.pfams) #number of pfams in each gene
max(lengths(genes.in.delta.pfams)) #8
genes.in.delta.pfams.df <- as.data.frame(t(stri_list2matrix(genes.in.delta.pfams))) #convert to dataframe, each pfam for each gene is a column
str(genes.in.delta.pfams.df)
genes.in.delta.pfam.split <- cbind(genes.in.delta, genes.in.delta.pfams.df)
str(genes.in.delta.pfam.split)
dim(genes.in.delta.pfam.split)

###Do fischers exact test
#two classes: 1)pfam all genes and 2)pfam genes with snps
#occurences of each pfam in each class

#make a count of pfams:
#1) all genes (with pfam):
n_distinct(genome.genes.pfam.split$V1) 
sapply(genome.genes.pfam.split, function(x) n_distinct(x))
genome.pfams <- c(genome.genes.pfam.split$V1,
                  genome.genes.pfam.split$V2,
                  genome.genes.pfam.split$V3,
                  genome.genes.pfam.split$V4,
                  genome.genes.pfam.split$V5,
                  genome.genes.pfam.split$V6,
                  genome.genes.pfam.split$V7,
                  genome.genes.pfam.split$V8,
                  genome.genes.pfam.split$V9,
                  genome.genes.pfam.split$V10,
                  genome.genes.pfam.split$V11,
                  genome.genes.pfam.split$V12,
                  genome.genes.pfam.split$V13)
str(genome.pfams)
length(genome.pfams)
genome.pfams.counts <- table(genome.pfams)
genome.pfams.counts <- as.data.frame(genome.pfams.counts)
colnames(genome.pfams.counts) <- c("pfam","freq_all")
str(genome.pfams.counts) #2680 unique pfams in total 4422 genes

#2) deltafst genes with pfam:
selected.genes.pfams <- c(genes.in.delta.pfam.split$V1,
                          genes.in.delta.pfam.split$V2,
                          genes.in.delta.pfam.split$V3,
                          genes.in.delta.pfam.split$V4,
                          genes.in.delta.pfam.split$V5,
                          genes.in.delta.pfam.split$V6,
                          genes.in.delta.pfam.split$V7,
                          genes.in.delta.pfam.split$V8,
                          genes.in.delta.pfam.split$V9,
                          genes.in.delta.pfam.split$V10,
                          genes.in.delta.pfam.split$V11)
length(selected.genes.pfams) 
selected.genes.pfams.counts <- table(selected.genes.pfams)
selected.genes.pfams.counts <- as.data.frame(selected.genes.pfams.counts)
colnames(selected.genes.pfams.counts) <- c("pfam", "freq_selected")
str(selected.genes.pfams.counts) #254 pfams in 370 deltafst genes
head(selected.genes.pfams.counts)

#combine 1) and 2):
combined.table.all <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam", all.x = TRUE)
dim(combined.table.all) #2680 pfams
head(combined.table.all,20)
combined.table.all[is.na(combined.table.all)] = 0 #replace NA with 0
head(combined.table.all,20)
fisher.table <- combined.table.all[,2:3]
rownames(fisher.table) <- combined.table.all$pfam
head(fisher.table)
fisher.table <- t(fisher.table)
dim(fisher.table)
test <- fisher.test(fisher.table, simulate.p.value=TRUE) #must use simulate to not exceed memory limitation
test #p-value=1

#Actually classes should be mutually exclusive so they can add up to total, so what you want is:
#1) all pfam in snp genes VERSUS all pfam not in snp genes
#2) all pfam outside snp genes VERSUS all pfam not outside snp genes
combined.table.delta <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam")
dim(combined.table.delta) #254: i.e. omits pfams not present in deltafst genes
head(combined.table.delta)
not.delta.table <- subset(combined.table.all, !(combined.table.all$pfam %in% combined.table.delta$pfam))
dim(not.delta.table) #2426 pfams not in snps genes, i.e. only outside snps
#pfam not found outside snp genes:
only.delta.table <- filter(combined.table.delta, freq_all == freq_selected)
dim(only.delta.table) #67 pfams only in deltafst genes, i.e. 67 pfams absent outside deltafst
table(genome.genes.pfam.split$V1 %in% only.delta.table$pfam) #23, does this equal to the number of genes?
table(genome.genes.pfam.split$V2 %in% only.delta.table$pfam) #16
table(genome.genes.pfam.split$V3 %in% only.delta.table$pfam) #15
table(genome.genes.pfam.split$V4 %in% only.delta.table$pfam) #8
table(genome.genes.pfam.split$V5 %in% only.delta.table$pfam) #2
table(genome.genes.pfam.split$V6 %in% only.delta.table$pfam) #3
table(genome.genes.pfam.split$V7 %in% only.delta.table$pfam) #1
table(genome.genes.pfam.split$V8 %in% only.delta.table$pfam) #0

table(not.delta.table$pfam %in% only.delta.table$pfam) #appears to not be overlapping
table(only.delta.table$pfam  %in% not.delta.table$pfam) #appears to not be overlapping

#Try to test pfam by pfam with counts, and then correct for multiple testing:
head(combined.table.all,3)
#counts pfam present outside snp genes: freq_all - freq_selected
#counts pfam present inside snp genes: freq_selected
#counts pfams absent outside snp genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside snp genes: 4422 - freq_selected

pvalues <- c() #loop through each pfam and do a Fischers exact test:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 185-freq_selected
  absent_outside_delta <- 4422-185-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_delta_gene","outside_delta_gene")
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction

#test instead subset of snp genes versus all genes (but does this violate the assumption of exclusivity of categories?)
#replace "outside_snp_genes" where pfam is absent with all genes (4422)
# but shouldn´t "outside_snp_genes" where pfam is present also be replaced with all genes?
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 185-freq_selected
  absent_outside_delta <- 4422-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 pfams significant

#go all the way by also replacing "outside_snp_gene" for pfam present with counts of total:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_delta <- freq_selected
  present_all <- freq_all
  absent_within_delta <- 1580-freq_selected
  absent_all <- 4422-freq_all
  genes_pfam_present <- c(present_within_delta, present_all)
  genes_pfam_absent <- c(absent_within_delta, absent_all)
  cont_table <- as.data.frame(cbind(genes_pfam_present, genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","all_genes")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #9 significant after correction

#exclude pfams that are not present in snp genes? (will make BH correction softer as fewer tests are conducted)
#but this is probably "unfair"...
pfam.delta.genes <- filter(combined.table.all, freq_selected > 0)
dim(pfam.delta.genes) #254 pfams present in 185 delta genes
head(pfam.delta.genes, 5)
a <- filter(pfam.delta.genes, freq_selected == freq_all)
dim(a) #67 pfams present only in snp genes
a <- filter(pfam.delta.genes, freq_selected < freq_all)
dim(a) #187 pfams shared snp genes and genes outside snps
str(genome.genes.pfam.split)
b <- filter(genome.genes.pfam.split, genome.genes.pfam.split$V1 %in% a$pfam)
dim(b)#3005 genes with shared pfam occuring in total, i.e. 3005-1260=1745 genes outside delta genes carrying shared pfams?

pvalues <- c()
for(i in 1:nrow(pfam.delta.genes)){
  freq_all <- pfam.delta.genes[i,2]
  freq_selected <- pfam.delta.genes[i,3]
  present_within_delta <- freq_selected
  present_outside_delta <- freq_all-freq_selected
  absent_within_delta <- 1260-freq_selected
  absent_outside_delta <- 1745-present_outside_delta
  genes_pfam_present <- c(present_within_delta,present_outside_delta)
  genes_pfam_absent <- c(absent_within_delta,absent_outside_delta)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_delta_gene","outside_delta_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
pfam.delta.genes$freq_outside <- pfam.delta.genes$freq_all - pfam.delta.genes$freq_selected
final.table <- cbind(pfam.delta.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #39 pfams significant after correction:

#try with pfam counts instead of gene counts:
sum(combined.table.all$freq_all) #7167 total counts of pfams
sum(combined.table.all$freq_selected) #327 total counts of pfams in deltafst
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within_snp <- freq_selected
  present_outside_snp <- freq_all-freq_selected
  absent_within_snp <- 327-freq_selected
  absent_outside <- 7167-327-present_within_snp-present_outside_snp
  genes_pfam_present <- c(present_within_snp,present_outside_snp)
  genes_pfam_absent <- c(absent_within_snp,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.snp.genes$freq_outside <- pfam.snp.genes$freq_all - pfam.snp.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0




###########################
#                         #
#     DeltaFst Asia       #
#                         #
###########################