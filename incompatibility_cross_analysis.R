setwd("/Users/dabaosl/Dropbox (UiO)/UiO/Phd/Data_analysis/")
library(tidyverse)

gwas.all.inferred <- read.table("Incompatibility_cross_analysis/res_amm_gwas_all_inferred_sorted_top1000.txt")
#gwas.all.inferred <- read.table("GWAS/output_amm/res_amm_gwas_all_inferred_sorted_top1000.txt")
str(gwas.all.inferred)
tail(gwas.all.inferred)

delta_fst.all <- read.table("Incompatibility_cross_analysis/Delta_df_all_outliers.txt", header=T)
#delta_fst.all <- read.table("Genome_scans/Delta_df_all_outliers.txt", header=T)
str(delta_fst.all)

delta_fst.noram <- read.table("Incompatibility_cross_analysis/Delta_df_NoramA_NoramB.txt", header=T)
#delta_fst.noram <- read.table("Genome_scans/Delta_df_NoramA_NoramB.txt", header=T)
str(delta_fst.noram) 

pos.sel.noramA <- read.table("Incompatibility_cross_analysis/noramA.sele", header=T)
#pos.sel.noramA <- read.table("Relate/Detect_positive_selection/noramA.sele", header=T)
str(pos.sel.noramA)
head(sort(10^pos.sel.noramA$when_mutation_has_freq2),20)
quantile(10^pos.sel.noramA$when_mutation_has_freq2, probs=0.02)
hist(10^pos.sel.noramA$when_mutation_has_freq2)
#need to filter out entries for which p-values cannot be calculated that has been set to -10log=1
pos.sel.noramA.filtered <- filter(pos.sel.noramA, when_mutation_has_freq2 < 1)
hist(10^pos.sel.noramA.filtered$when_mutation_has_freq2)

pos.sel.noramB <- read.table("Incompatibility_cross_analysis/noramBE.sele", header=T)
#pos.sel.noramB <- read.table("Relate/Detect_positive_selection/noramBE.sele", header=T)
head(sort(10^pos.sel.noramB$when_mutation_has_freq2),20)
quantile(10^pos.sel.noramB$when_mutation_has_freq2, probs=0.02)
pos.sel.noramB.filtered <- filter(pos.sel.noramB, when_mutation_has_freq2 < 1)

pos.sel.eurasia <- read.table("Incompatibility_cross_analysis/Eurasia.sele", header=T)
#pos.sel.eurasia <- read.table("Relate/Detect_positive_selection/Eurasia.sele", header=T)
head(sort(10^pos.sel.eurasia$when_mutation_has_freq2),20)
mean(sort(10^pos.sel.eurasia$when_mutation_has_freq2),20)
pos.sel.eurasia.filtered <- filter(pos.sel.eurasia, when_mutation_has_freq2 < 1)

pos.sel.euro <- read.table("Incompatibility_cross_analysis/Euro.sele", header=T)
#pos.sel.euro <- read.table("Relate/Detect_positive_selection/Euro.sele", header=T)
head(sort(10^pos.sel.euro$when_mutation_has_freq2),20)
pos.sel.euro.filtered <- filter(pos.sel.euro, when_mutation_has_freq2 < 1)

go.delta.all <- read.table("Incompatibility_cross_analysis/Delta_GOenrichment.tab", header=T, sep="\t")
#go.delta.all <- read.table("GeneOntology_analysis/Output_20220824/Delta_GOenrichment.tab", header=T, sep="\t")
str(go.delta.all)
go.delta.noram <- read.table("Incompatibility_cross_analysis/DeltadfNoramANoramB.txt_GOenrichment.tab", header=T, sep="\t")
#go.delta.noram <- read.table("GeneOntology_analysis/Output_20220824/DeltadfNoramANoramB.txt_GOenrichment.tab", header=T, sep="\t")
str(go.delta.noram) 

scopa.crosses.noram <- read.table("Incompatibility_cross_analysis/SCOPA_crosses_NorAm_num_scf_gen_result_simple_sorted.txt", header=T)
#scopa.crosses.noram <- read.table("GWAS/SCOPA_crosses_NorAm_num_scf_gen_result_simple_sorted.txt", header=T)
str(scopa.crosses.noram)
scopa.crosses.subset <- read.table("Incompatibility_cross_analysis/SCOPA_crosses_result_simple_sorted.txt", header=T)
#scopa.crosses.subset <- read.table("GWAS/SCOPA_crosses_result_simple_sorted.txt", header=T)
str(scopa.crosses.subset)

#scopa.asia <- read.table("GWAS/", header=T) result not available

caspases <- read.table("Incompatibility_cross_analysis/caspases_TA10106M1_final_MAKER.putative_function_domains.gff", sep="\t")
str(caspases)

rcd1 <- read.table("Incompatibility_cross_analysis/rcd_TA10106M1_final_MAKER.putative_function_domains.gff", sep="\t")
str(rcd1)

pfam05729 <- read.table("Incompatibility_cross_analysis/PF05729_TA10106M1_final_MAKER.putative_function_domains.gff",sep="\t")
str(pfam05729)

pfam00400 <- read.table("Incompatibility_cross_analysis/WD40_PF00400_TA10106M1_final_MAKER.putative_function_domains.gff", sep="\t")
str(pfam00400)
head(pfam00400,10)
pfam004000.gene <- filter(pfam00400, V3=="gene")
str(pfam004000.gene)
table(pfam004000.gene$V1)
write.table(pfam004000.gene, file="pfam004000_gene.txt")
str(pfam004000.gene)

het_s <- read.table("Incompatibility_cross_analysis/het_s_podospora_tblastn_hits.txt")
str(het_s)
het_s
het_S <- read.table("Incompatibility_cross_analysis/het_bigS_podospora_tblastn_hits.txt")
str(het_S)
het_S
het_c <- read.table("Incompatibility_cross_analysis/het_c_podospora_tblastn_hits.txt")
str(het_c)
het_c
het_c2 <- read.table("Incompatibility_cross_analysis/het_c2_podospora_tblastn_hits.txt")
str(het_c2)
het_c2
un24 <- read.table("Incompatibility_cross_analysis/un24_neurospora_tblastn_hits.txt") 
str(un24)  
un24
blast.hits <- rbind(het_s, het_S, het_c, het_c2, un24)
str(blast.hits)
colnames(blast.hits) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
str(blast.hits)
blast.hits.sign <- filter(blast.hits, evalue < 1e-05)
blast.hits.sign

het_r <- read.table("Incompatibility_cross_analysis/het_r_podospora_tblastn_hits.txt")
str(het_r)
het_r.sign <- filter(het_r, V11 < 1e-05)
dim(het_r.sign) #784

het_d <- read.table("Incompatibility_cross_analysis/het_d_podospora_tblastn_hits.txt")
str(het_d)
het_d.sign <- filter(het_d, V11 < 1e-05)
dim(het_d.sign) #346

het_e <- read.table("Incompatibility_cross_analysis/het_e_podospora_tblastn_hits.txt")
str(het_e)
het_e.sign <- filter(het_e, V11 < 1e-05)
dim(het_e.sign) #171

het_6 <- read.table("Incompatibility_cross_analysis/het_6_neurospora_tetrasperma_tblastn_hits.txt")
str(het_6)
het_6.sign <- filter(het_6, V11 < 1e-05)
dim(het_6.sign) #0

blast.hits.more

####check blast hits first:
#in deltafst all:
str(delta_fst.all)

blast.test <- filter(delta_fst.all, Scaffold=="Scaffold01" & Start > 35000)
blast.test
blast.test <- filter(blast.test, End < 35000 + 11000)
blast.test
dim(blast.test)

blast.hits.sign[1,]
blast1 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 2809158)
blast1 <- filter(blast1, End < 2809158 + 11000)
blast1 #1 hit
#entire gene
blast1 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 2809704)
blast1 <- filter(blast1, End < 2809704 + 11000)
blast1 #0 hit

blast.hits.sign[2,]
blast2 <- filter(delta_fst.all, Scaffold=="Scaffold04" & Start > 1470389)
blast2 <- filter(blast2, End < 1470389 + 11000 )
blast2 #0 hits
#entire gene
blast1 <- filter(delta_fst.all, Scaffold=="Scaffold04" & Start > 1470202)
blast1 <- filter(blast1, End < 1470202 + 11000)
blast1 #0 hit

blast.hits.sign[3,]
blast3 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 2809158)
blast3 <- filter(blast3, End < 2809158 + 11000)
blast3 #1 hit (same as blast hit1)

blast.hits.sign[4,]
blast4 <- filter(delta_fst.all, Scaffold=="Scaffold04" & Start > 1470389)
blast4 <- filter(blast4, End < 1470389 + 11000)
blast4 #0 hits

blast.hits.sign[5,]
blast5 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 3221045)
blast5 <- filter(blast5, End < 3221045 + 11000)
blast5 #0 rows

blast.hits.sign[6,]
blast6 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 3222693)
blast6 <- filter(blast6, End < 3222693 + 11000)
blast6 #0 rows

blast.hits.sign[7,]
blast7 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 3223523)
blast7 <- filter(blast7, End < 3223523 + 11000)
blast7 #0 rows

blast.hits.sign[8,]
blast8 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 3223414)
blast8 <- filter(blast8, End < 3223414 + 11000)
blast8 #O rows
#entire gene:
blast8 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 3214538)
blast8 <- filter(blast1, End < 3214538 + 11000)
blast8 #0 hit

blast.hits.sign[9,]
blast9 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1888969)
blast9 <- filter(blast9, End < 1888969 + 11000)
blast9 #0 rows

blast.hits.sign[10,]
blast10 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1887921)
blast10 <- filter(blast10, End < 1887921 + 11000)
blast10 #0 rows

blast.hits.sign[11,]
blast11 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1888724)
blast11 <- filter(blast11, End < 1888724 + 11000)
blast11 #0 rows

blast.hits.sign[12,]
blast12 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1887727)
blast12 <- filter(blast12, End < 1887727 + 11000)
blast12 #0 rows

blast13 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1887494)
blast13 <- filter(blast13, End < 1887494 + 11000)
blast13  #0 rows
#entire gene:
blast13 <- filter(delta_fst.all, Scaffold=="Scaffold07" & Start > 1887259)
blast13 <- filter(blast13, End < 1887259 + 11000)
blast13  

#in deltafst noram:
blast1 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 2809158)
blast1 <- filter(blast1, End < 2809158 + 11000)
blast1 #0 hit
#entire gene
blast1 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 2809704)
blast1 <- filter(blast1, End < 2809704 + 11000)
blast1 #0 hit

blast2 <- filter(delta_fst.noram, Scaffold=="Scaffold04" & Start > 1470389)
blast2 <- filter(blast2, End < 1470389 + 11000 )
blast2 #0 hits
#entire gene
blast2 <- filter(delta_fst.noram, Scaffold=="Scaffold04" & Start > 1470202)
blast2 <- filter(blast1, End < 1470202 + 11000)
blast2 #0 hit

blast3 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 2809158)
blast3 <- filter(blast3, End < 2809158 + 11000)
blast3 #0 hits

blast4 <- filter(delta_fst.noram, Scaffold=="Scaffold04" & Start > 1470389)
blast4 <- filter(blast4, End < 1470389 + 11000)
blast4 #0 hits

blast5 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 3221045)
blast5 <- filter(blast5, End < 3221045 + 11000)
blast5 #0 rows

blast6 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 3222693)
blast6 <- filter(blast6, End < 3222693 + 11000)
blast6 #0 rows

blast7 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 3223523)
blast7 <- filter(blast7, End < 3223523 + 11000)
blast7 #0 rows

blast8 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 3223414)
blast8 <- filter(blast8, End < 3223414 + 11000)
blast8 #O rows
#entire gene:
blast8 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 3214538)
blast8 <- filter(blast1, End < 3214538 + 11000)
blast8 #0 hit

blast9 <- filter(delta_fst.noram, Scaffold=="Scaffold07" & Start > 1888969)
blast9 <- filter(blast9, End < 1888969 + 11000)
blast9 #0 rows

blast10 <- filter(delta_fst.noram, Scaffold=="Scaffold07" & Start > 1887921)
blast10 <- filter(blast10, End < 1887921 + 11000)
blast10 #0 rows

blast11 <- filter(delta_fst.noram, Scaffold=="Scaffold07" & Start > 1888724)
blast11 <- filter(blast11, End < 1888724 + 11000)
blast11 #0 rows

blast12 <- filter(delta_fst.noram, Scaffold=="Scaffold07" & Start > 1887727)
blast12 <- filter(blast12, End < 1887727 + 11000)
blast12 #0 rows
#entire gene:
blast13 <- filter(delta_fst.noram, Scaffold=="Scaffold07" & Start > 1887259)
blast13 <- filter(blast13, End < 1887259 + 11000)
blast13 

##in positive selection:
#NoramA
str(pos.sel.noramA)
blast.hits.sign[1,]
blast1 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast1) #6 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
blast1 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 2809704 & pos < 2810125)
dim(blast1) #4 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[2,]
blast2 <- filter(pos.sel.noramA, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast2) #2 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast2 <- filter(pos.sel.noramA, rs_id=="Scaffold04" & pos > 1470202 & pos < 1471894)
dim(blast2)
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[3,]
blast3 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast3) #6 SNPs
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[4,]
blast4 <- filter(pos.sel.noramA, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast4) #2 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[5,]
blast5 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3221045 & pos < 3222646)
dim(blast5) #5 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[6,]
blast6 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3222693 & pos < 3223364)
dim(blast6) #10 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
blast6$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[7,]
blast7 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3223414  & pos < 3223470)
dim(blast7) #0 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[8,]
blast8 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3223414 & pos < 3223470)
dim(blast8) #0 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast8 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3214538 & pos < 3227385)
dim(blast8) #75 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
table(blast8$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

blast.hits.sign[9,]
blast9 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1888969 & pos < 1890018)
dim(blast9) #2 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
blast9 #1 SNP with p-value 0.048
blast9$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[10,]
blast10 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1887921 & pos < 1888673 )
dim(blast10) #1 SNPs
blast10$p.value.freq2 <- 10^blast10$when_mutation_has_freq2 
blast10$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[11,]
blast11 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1888724 & pos < 1888948 )
dim(blast11) #0 SNPs
blast11$p.value.freq2 <- 10^blast11$when_mutation_has_freq2 
blast11$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[12,]
blast12 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1887727 & pos < 1887867)
dim(blast12) #0 SNPs
blast12$p.value.freq2 <- 10^blast12$when_mutation_has_freq2 
blast12$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[13,]
blast13 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1887494 & pos < 1887673)
dim(blast13) #0 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast13 <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1887259 & pos < 1892493)
dim(blast13) #8 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
table(blast13$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

##NoramB
str(pos.sel.noramB)
blast.hits.sign[1,]
blast1 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast1) #5 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1 #1 where p-value is 0.02
blast1$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast1 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 2809704 & pos < 2810125)
dim(blast1) #2 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[2,]
blast2 <- filter(pos.sel.noramB, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast2) #2 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2 #1 SNP with p-value 0.027
blast2$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast2 <- filter(pos.sel.noramB, rs_id=="Scaffold04" & pos > 1470202 & pos < 1471894)
dim(blast2)
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[3,]
blast3 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast3) #5 SNPs
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3 #1 SNP p-value 0.027
blast3$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[4,]
blast4 <- filter(pos.sel.noramB, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast4) #2 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #1 SNP p-value 0.027
blast4$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[5,]
blast5 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3221045 & pos < 3222646)
dim(blast5) #5 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5 
blast5$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[6,]
blast6 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3222693 & pos < 3223364)
dim(blast6) #3 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
blast6
blast6$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[7,]
blast7 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3223414  & pos < 3223470)
dim(blast7) #0 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7
blast7$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[8,]
blast8 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3223414 & pos < 3223470)
dim(blast8) #0 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8
blast8$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast8 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3214538 & pos < 3227385)
dim(blast8) #72 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
table(blast8$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

blast.hits.sign[9,]
blast9 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1888969 & pos < 1890018)
dim(blast9) #2 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
blast9 #2 SNPs with p-value 0.027
blast9$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[10,]
blast10 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1887921 & pos < 1888673 )
dim(blast10) #1 SNPs
blast10$p.value.freq2 <- 10^blast10$when_mutation_has_freq2 
blast10 #1 SNP p-value 0.027
blast10$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[11,]
blast11 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1888724 & pos < 1888948 )
dim(blast11) #2 SNPs
blast11$p.value.freq2 <- 10^blast11$when_mutation_has_freq2 
blast11
blast11$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[12,]
blast12 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1887727 & pos < 1887867)
dim(blast12) #0 SNPs
blast12$p.value.freq2 <- 10^blast12$when_mutation_has_freq2 
blast12
blast12$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[13,]
blast13 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1887494 & pos < 1887673)
dim(blast13) #0 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
blast13$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast13 <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1887259 & pos < 1892493)
dim(blast13) #83 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
table(blast13$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

#Eurasia
str(pos.sel.eurasia)
blast.hits.sign[1,]
blast1 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast1) #3 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1 #1 where p-value is 0.02
blast1$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast1 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 2809704 & pos < 2810125)
dim(blast1) #0 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[2,]
blast2 <- filter(pos.sel.eurasia, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast2) #3 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2 #0 significant
blast2$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast2 <- filter(pos.sel.eurasia, rs_id=="Scaffold04" & pos > 1470202 & pos < 1471894)
dim(blast2)
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[3,]
blast3 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast3) #3 SNPs
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3 #1 SNP p-value 0.02
blast3$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[4,]
blast4 <- filter(pos.sel.eurasia, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast4) #3 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #0 significant
blast4$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[5,]
blast5 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3221045 & pos < 3222646)
dim(blast5) #15 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5  #3 significant
blast5$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[6,]
blast6 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3222693 & pos < 3223364)
dim(blast6) #8 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
blast6 #1 somehow significant
blast6$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[7,]
blast7 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3223414  & pos < 3223470)
dim(blast7) #0 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7
blast7$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[8,]
blast8 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3223414 & pos < 3223470)
dim(blast8) #0 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8
blast8$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast8 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3214538 & pos < 3227385)
dim(blast8) #96 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
table(blast8$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
#4

blast.hits.sign[9,]
blast9 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1888969 & pos < 1890018)
dim(blast9) #26 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
blast9 #1 SNPs with p-value 0.02
blast9$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[10,]
blast10 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1887921 & pos < 1888673 )
dim(blast10) #18 SNPs
blast10$p.value.freq2 <- 10^blast10$when_mutation_has_freq2 
blast10 #1 SNP p-value 0.03
blast10$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[11,]
blast11 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1888724 & pos < 1888948 )
dim(blast11) #5 SNPs
blast11$p.value.freq2 <- 10^blast11$when_mutation_has_freq2 
blast11
blast11$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[12,]
blast12 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1887727 & pos < 1887867)
dim(blast12) #5 SNPs
blast12$p.value.freq2 <- 10^blast12$when_mutation_has_freq2 
blast12
blast12$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[13,]
blast13 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1887494 & pos < 1887673)
dim(blast13) #0 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
blast13$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast13 <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1887259 & pos < 1892493)
dim(blast13) #83 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
table(blast13$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

#Euro
str(pos.sel.euro)
blast.hits.sign[1,]
blast1 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast1) #3 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
blast1 #1 where p-value is 0.02
blast1$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast1 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 2809704 & pos < 2810125)
dim(blast1) #0 SNPs

blast.hits.sign[2,]
blast2 <- filter(pos.sel.euro, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast2) #3 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2 #0 significant
blast2$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast2 <- filter(pos.sel.euro, rs_id=="Scaffold04" & pos > 1470202 & pos < 1471894)
dim(blast2)
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[3,]
blast3 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 2809158 & pos < 2809970)
dim(blast3) #3 SNPs
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3 #1 SNP p-value 0.02
blast3$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[4,]
blast4 <- filter(pos.sel.euro, rs_id=="Scaffold04" & pos > 1470389 & pos < 1470730)
dim(blast4) #3 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #0 significant
blast4$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[5,]
blast5 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3221045 & pos < 3222646)
dim(blast5) #15 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5  #3 significant
blast5$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[6,]
blast6 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3222693 & pos < 3223364)
dim(blast6) #8 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
blast6 #1 somehow significant
blast6$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[7,]
blast7 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3223414  & pos < 3223470)
dim(blast7) #0 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7
blast7$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[8,]
blast8 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3223414 & pos < 3223470)
dim(blast8) #0 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8
blast8$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene:
blast8 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3214538 & pos < 3227385)
dim(blast8) #75 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
table(blast8$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

blast.hits.sign[9,]
blast9 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1888969 & pos < 1890018)
dim(blast9) #26 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
blast9 #1 SNPs with p-value 0.02
blast9$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[10,]
blast10 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1887921 & pos < 1888673 )
dim(blast10) #18 SNPs
blast10$p.value.freq2 <- 10^blast10$when_mutation_has_freq2 
blast10 #1 SNP p-value 0.03
blast10$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[11,]
blast11 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1888724 & pos < 1888948 )
dim(blast11) #5 SNPs
blast11$p.value.freq2 <- 10^blast11$when_mutation_has_freq2 
blast11
blast11$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[12,]
blast12 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1887727 & pos < 1887867)
dim(blast12) #5 SNPs
blast12$p.value.freq2 <- 10^blast12$when_mutation_has_freq2 
blast12
blast12$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)

blast.hits.sign[13,]
blast13 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1887494 & pos < 1887673)
dim(blast13) #0 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
blast13$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#entire gene
blast13 <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1887259 & pos < 1892493)
dim(blast13) #42 SNPs
blast13$p.value.freq2 <- 10^blast13$when_mutation_has_freq2 
blast13
table(blast13$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
#
##in GWAS
str(gwas.all.inferred)
tail(gwas.all.inferred)

blast.hits.sign[1,]
blast1 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 2809158 & Pos < 2809970)
(blast1)
dim(blast1) #0 SNPs
#try entire gene:
blast1 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 2809704 & Pos < 2810125)
table(blast1$P.value < 1e-5)
blast1

blast.hits.sign[2,]
blast2 <- filter(gwas.all.inferred, Chr=="Scaffold04" & Pos > 1470389 & Pos < 1470730)
dim(blast2) #0 snps
#try entire gene:
blast2 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold04", Pos > 1470202 & Pos < 1471894)
table(blast2$P.value < 1e-5)
blast2

blast.hits.sign[3,]
blast3 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 2809158 & Pos < 2809970)
dim(blast3) #0 SNPs

blast.hits.sign[4,]
blast4 <- filter(gwas.all.inferred, Chr=="Scaffold04" & Pos > 1470389 & Pos < 1470730)
dim(blast4) #0 SNPs

blast.hits.sign[5,]
blast5 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 3221045 & Pos < 3222646)
dim(blast5) #0 SNPs

blast.hits.sign[6,]
blast6 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 3222693 & Pos < 3223364)
dim(blast6) #0 SNPs

blast.hits.sign[7,]
blast7 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 3223414  & Pos < 3223470)
dim(blast7) #0 SNPs

blast.hits.sign[8,]
blast8 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 3223414 & Pos < 3223470)
dim(blast8) #0 SNPs
#try entire gene
blast5678  <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 3214538 & Pos < 3227385)
table(blast5678$P.value < 1e-5)

blast.hits.sign[9,]
blast9 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1888969 & Pos < 1890018)
dim(blast9) #0 SNPs

blast.hits.sign[10,]
blast10 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1887921 & Pos < 1888673 )
dim(blast10) #0 SNPs

blast.hits.sign[11,]
blast11 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1888724 & Pos < 1888948 )
dim(blast11) #0 SNPs

blast.hits.sign[12,]
blast12 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1887727 & Pos < 1887867)
dim(blast12) #0 SNPs

blast.hits.sign[13,]
blast13 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1887494 & Pos < 1887673)
dim(blast13) #0 SNPs
#entire gene
blast9to13 <- filter(gwas.all.inferred, Chr=="Scaffold07" & Pos > 1887259 & Pos < 1892493)
table(blast9to13$P.value < 1e-5)

head(gwas.all.inferred)
blast.test <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 1122517  & Pos < 1122600)
dim(blast.test) #works in theory


##SCOPA
#SCOPA noram
blast.hits.sign[1,]
blast1 <- filter(scopa.crosses.noram, Chromosome== 5, Position > 2809158 & Position < 2809970)
dim(blast1) #20 SNPs
blast1 #3 with low p-value
table(blast1$P.value < 1e-5)

blast.hits.sign[2,]
blast2 <- filter(scopa.crosses.noram, Chromosome==4 & Position > 1470389 & Position < 1470730)
dim(blast2) #11 snps
blast2 #1 somewhat low
table(blast2$P.value < 1e-5)

blast.hits.sign[3,]
blast3 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 2809158 & Position < 2809970)
dim(blast3) #3 SNPs
blast3 #3 low p-values
table(blast3$P.value < 1e-5)

blast.hits.sign[4,]
blast4 <- filter(scopa.crosses.noram, Chromosome==4 & Position > 1470389 & Position < 1470730)
dim(blast4) #4 SNPs
blast4 #1 somewhat low
table(blast4$P.value < 1e-5)

blast.hits.sign[5,]
blast5 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 3221045 & Position < 3222646)
dim(blast5) #70 SNPs
blast5 #1 quite low and a couple somewhat low
table(blast5$P.value < 1e-5)

blast.hits.sign[6,]
blast6 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 3222693 & Position < 3223364)
dim(blast6) #41 SNPs
blast6 #3 quite low
table(blast6$P.value < 1e-5)

blast.hits.sign[7,]
blast7 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 3223414  & Position < 3223470)
dim(blast7) #3 SNPs
blast7
table(blast6$P.value < 1e-5)

blast.hits.sign[8,]
blast8 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 3223414 & Position < 3223470)
dim(blast8) #3 SNPs
blast8
table(blast8$P.value < 1e-5)

blast.hits.sign[9,]
blast9 <- filter(scopa.crosses.noram, Chromosome==7 & Position > 1888969 & Position < 1890018)
dim(blast9) #57 SNPs
blast9 #6 quite low
table(blast9$P.value < 1e-5)

blast.hits.sign[10,]
blast10 <- filter(scopa.crosses.noram, Chromosome==7 & Position > 1887921 & Position < 1888673 )
dim(blast10) #56 SNPs
blast10 #1 quite low, a couple somewhat low
table(blast10$P.value < 1e-5)

blast.hits.sign[11,]
blast11 <- filter(scopa.crosses.noram, Chromosome==7 & Position > 1888724 & Position < 1888948 )
dim(blast11) #16 SNPs
blast11 #1 quite low, a couple somewhat low
table(blast11$P.value < 1e-5)

blast.hits.sign[12,]
blast12 <- filter(scopa.crosses.noram, Chromosome==7 & Position > 1887727 & Position < 1887867)
dim(blast12) #12 SNPs
blast12 #2 quite low
table(blast12$P.value < 1e-5)

blast.hits.sign[13,]
blast13 <- filter(scopa.crosses.noram, Chromosome==7 & Position > 1887494 & Position < 1887673)
dim(blast13) #9 SNPs
blast13 # a few somewhat low
table(blast13$P.value < 1e-5)

str(scopa.crosses.noram)
head(scopa.crosses.noram)
blast.test <- filter(scopa.crosses.noram, Chromosome==1 & Position > 2000510  & Position < 2000518)
blast.test
dim(blast.test) #works in theory

#SCOPA all crosses subset
blast.hits.sign[1,]
blast1 <- filter(scopa.crosses.subset, Chromosome== 5, Position > 2809158 & Position < 2809970)
dim(blast1) #20 SNPs
blast1 #5 with low p-value
#try entire gene:
blast1 <- filter(scopa.crosses.subset, Chromosome== 5, Position > 2809704 & Position < 2810125)
table(blast1$P.value < 1e-5)
blast1

blast.hits.sign[2,]
blast2 <- filter(scopa.crosses.subset, Chromosome==4 & Position > 1470389 & Position < 1470730)
dim(blast2) #11 snps
blast2 #1 very low, 7 quite low
#try entire gene:
blast2 <- filter(scopa.crosses.subset, Chromosome== 4, Position > 1470202 & Position < 1471894)
table(blast2$P.value < 1e-5)
blast2

blast.hits.sign[3,]
blast3 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 2809158 & Position < 2809970)
dim(blast3) #20 SNPs
blast3 #5 low p-values

blast.hits.sign[4,]
blast4 <- filter(scopa.crosses.subset, Chromosome==4 & Position > 1470389 & Position < 1470730)
dim(blast4) #11 SNPs
blast4 #3 very low, 5 quite low

blast.hits.sign[5,]
blast5 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 3221045 & Position < 3222646)
dim(blast5) #70 SNPs
blast5 #11 SNPs very low and many somewhat low

blast.hits.sign[6,]
blast6 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 3222693 & Position < 3223364)
dim(blast6) #41 SNPs
blast6 #1 very low, a few somewhat low

blast.hits.sign[7,]
blast7 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 3223414  & Position < 3223470)
dim(blast7) #3 SNPs
blast7

blast.hits.sign[8,]
blast8 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 3223414 & Position < 3223470)
dim(blast8) #3 SNPs
blast8 

#5,6,7,8 is same gene, could be reduced to one search:
blast5678  <- filter(scopa.crosses.subset, Chromosome==5 & Position > 3214538 & Position < 3227385)
table(blast5678$P.value < 1e-5)

blast.hits.sign[9,]
blast9 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1888969 & Position < 1890018)
dim(blast9) #57 SNPs
blast9 #3 very low, many somewhat low

blast.hits.sign[10,]
blast10 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1887921 & Position < 1888673 )
dim(blast10) #56 SNPs
blast10 #2 very low many somewhat low

blast.hits.sign[11,]
blast11 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1888724 & Position < 1888948 )
dim(blast11) #16 SNPs
blast11 #1 quite low, a couple somewhat low

blast.hits.sign[12,]
blast12 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1887727 & Position < 1887867)
dim(blast12) #12 SNPs
blast12 #1 very low , two somewhat low

blast.hits.sign[13,]
blast13 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1887494 & Position < 1887673)
dim(blast13) #9 SNPs
blast13 # 2 very low, a few somewhat low

#9,10,11,12,13 is same gene:
blast9to13 <- filter(scopa.crosses.subset, Chromosome==7 & Position > 1887259 & Position < 1892493)
table(blast9to13$P.value < 1e-5)
blast9to13

str(scopa.crosses.noram)
head(scopa.crosses.noram)
blast.test <- filter(scopa.crosses.noram, Chromosome==1 & Position > 2000510  & Position < 2000518)
blast.test
dim(blast.test) #works in theory

###check domain
pfam05729
#deltafst all
pfam05729.delta.fst.all <- filter(delta_fst.all, Scaffold=="Scaffold06" & Start > 1439976 & End < 1439976 + 11000)
pfam05729.delta.fst.all #0
#deltafst noram
pfam05729.delta.fst.noram <- filter(delta_fst.noram, Scaffold=="Scaffold06" & Start > 1439976 & End < 1439976 + 11000)
pfam05729.delta.fst.noram #0
#pos.sel noramA
pfam05729.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold06" & pos > 1439976  & pos < 1445370)
dim(pfam05729.pos.sel.noramA) #4 SNPs
pfam05729.pos.sel.noramA$pvalue <- 10^pfam05729.pos.sel.noramA$when_mutation_has_freq2 
pfam05729.pos.sel.noramA #3 somewhat significant
pfam05729.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#pos sel noramB
pfam05729.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold06" & pos > 1439976  & pos < 1445370)
dim(pfam05729.pos.sel.noramB) #12 SNPs
pfam05729.pos.sel.noramB$pvalue <- 10^pfam05729.pos.sel.noramB$when_mutation_has_freq2 
pfam05729.pos.sel.noramB #3 or 4 somewhat significant
pfam05729.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#pos sel Eurasia
pfam05729.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold06" & pos > 1439976  & pos < 1445370)
dim(pfam05729.pos.sel.eurasia) #17 SNPs
pfam05729.pos.sel.eurasia$pvalue <- 10^pfam05729.pos.sel.eurasia$when_mutation_has_freq2 
pfam05729.pos.sel.eurasia #several somewhat significant
pfam05729.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05)
#pos sel Euro:
pfam05729.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold06" & pos > 1439976  & pos < 1445370)
dim(pfam05729.pos.sel.euro) #8 SNPs
pfam05729.pos.sel.euro$pvalue <- 10^pfam05729.pos.sel.euro$when_mutation_has_freq2 
pfam05729.pos.sel.euro 
pfam05729.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#gwas
pfam05729.gwas <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold06", Pos > 1439976 & Pos < 1445370)
dim(pfam05729.gwas) #0 SNPs
#scopa crosses noram
pfam05729.scopa.noram <- filter(scopa.crosses.noram, Chromosome==6 & Position > 1439976  & Position < 1445370)
dim(pfam05729.scopa.noram) #45 SNPs
pfam05729.scopa.noram #some somewhat significant
#scopa crosses subset
pfam05729.scopa.subset <- filter(scopa.crosses.subset, Chromosome==6 & Position > 1439976  & Position < 1445370)
dim(pfam05729.scopa.subset) #45 SNPs
pfam05729.scopa.subset #many somewhat significant

#pfam004000
str(pfam004000.gene)
str(delta_fst.noram)
#delta.fst.noram
wd40.delta.fst.noram <-left_join(pfam004000.gene, delta_fst.noram, by = c("V1"="Scaffold")) %>% filter(V4 >= Start, V5 <= End)
str(wd40.delta.fst.noram)
#4 genes
head(wd40.delta.fst.noram,10)
# Scaffold05 maker gene 3535366 3540285  .  +  .
# Scaffold05 maker gene 3710714 3711293  .  -  .
# Scaffold10 maker gene 2080746 2088064  .  +  .
# Scaffold01 maker gene 1970115 1972286  .  -  .

#delta.fst.all
wd40.delta.fst.all <-left_join(pfam004000.gene, delta_fst.all, by = c("V1"="Scaffold")) %>% filter(V4 >= Start, V5 <= End)
str(wd40.delta.fst.all)
#24 genes
head(wd40.delta.fst.all,24)
wd40.delta.fst.all[3,]
#1  Scaffold03 maker gene 1180329 1186588  .  -  .
#2  Scaffold03 maker gene 2277944 2280931  .  +  .
#3  Scaffold03 maker gene 3356359 3357989  .  +  .
#4  Scaffold03 maker gene 4488488 4494271  .  +  .
#5  Scaffold06 maker gene  289366  291920  .  -  .
#6  Scaffold05 maker gene 1895702 1897213  .  -  .
#7  Scaffold05 maker gene 2261246 2262456  .  +  .
#8  Scaffold05 maker gene 2261246 2262456  .  +  .
#9  Scaffold05 maker gene 3535366 3540285  .  +  .
#10 Scaffold05 maker gene 4883313 4884898  .  +  .
#11 Scaffold12 maker gene 1300977 1305868  .  +  .
#12 Scaffold10 maker gene  374220  374936  .  +  .
#13 Scaffold10 maker gene  383820  384862  .  -  .
#14 Scaffold01 maker gene 2130910 2132724  .  -  .
#15 Scaffold01 maker gene 3309095 3313763  .  +  .
#16 Scaffold07 maker gene 1882961 1887017  .  -  .
#17 Scaffold07 maker gene 2495402 2499338  .  +  .
#18 Scaffold02 maker gene  510389  513008  .  +  .
#19 Scaffold02 maker gene 2185478 2192354  .  +  .
#20 Scaffold02 maker gene 2794003 2795777  .  +  .
#21 Scaffold02 maker gene 2822363 2824510  .  -  .
#22 Scaffold04 maker gene 2352891 2354877  .  +  .
#23 Scaffold04 maker gene 4515977 4517169  .  -  .
#24 Scaffold04 maker gene 4515977 4517169  .  -  .

#7 and 8 redundant, 23 and 24 redundant

tail(wd40.delta.fst.all,4)
#gwas
wd40.gwas <-left_join(pfam004000.gene, gwas.all.inferred, by = c("V1"="Chr")) %>% filter(Pos  >= V4, Pos <= V5)

#23 genes
head(wd40.gwas,23)
#1  Scaffold03 maker gene 2639627 2645103  .  -  .
#2  Scaffold03 maker gene 3166781 3170573  .  +  .
#3  Scaffold05 maker gene 1895702 1897213  .  -  .
#4  Scaffold05 maker gene 1895702 1897213  .  -  .
#5  Scaffold05 maker gene 2261246 2262456  .  +  .
#6  Scaffold05 maker gene 3535366 3540285  .  +  .
#7  Scaffold05 maker gene 3535366 3540285  .  +  .
#8  Scaffold05 maker gene 3535366 3540285  .  +  .
#9  Scaffold05 maker gene 3535366 3540285  .  +  .
#10 Scaffold05 maker gene 3535366 3540285  .  +  .
#11 Scaffold05 maker gene 3535366 3540285  .  +  .
#12 Scaffold10 maker gene 1377956 1391649  .  -  .
#13 Scaffold10 maker gene 1377956 1391649  .  -  .
#14 Scaffold01 maker gene  344473  349967  .  -  .
#15 Scaffold01 maker gene 1636850 1655600  .  -  .
#16 Scaffold01 maker gene 1970115 1972286  .  -  .
#17 Scaffold07 maker gene 2495402 2499338  .  +  .
#18 Scaffold07 maker gene 2495402 2499338  .  +  .
#19 Scaffold07 maker gene 2495402 2499338  .  +  .
#20 Scaffold02 maker gene 1348871 1351951  .  +  .
#21 Scaffold02 maker gene 2185478 2192354  .  +  .
#22 Scaffold02 maker gene 2185478 2192354  .  +  .
#23 Scaffold02 maker gene 2185478 2192354  .  +  .

#wd40 noram
wd40.gwas.1 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 3535366 & Pos < 3540285)
dim(wd40.gwas.1)
wd40.gwas.1 #some quite significant
wd40.gwas.2 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 3710714 & Pos < 3711293)
wd40.gwas.2 #0 SNPs
wd40.gwas.3 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold10", Pos > 2080746 & Pos < 2088064)
wd40.gwas.3 #0 SNPs
wd40.gwas.4 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold01", Pos > 1970115 & Pos < 1972286)
wd40.gwas.4 #1 SNPs

#wd40 all:
wd40.gwas.1 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold03", Pos > 1180329 & Pos < 1186358)
dim(wd40.gwas.1) #0
wd40.gwas.2 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold03", Pos > 2277944 & Pos < 2280931)
wd40.gwas.2 #0 SNPs
wd40.gwas.3 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold03", Pos > 3356359 & Pos < 3357989)
wd40.gwas.3 #0 SNPs
wd40.gwas.4 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold03", Pos > 4488488 & Pos < 4494271)
wd40.gwas.4 #0 SNPs
wd40.gwas.5 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold06", Pos > 289366 & Pos < 291920)
wd40.gwas.5 #0 SNPs
wd40.gwas.6 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 1895702 & Pos < 1897213)
wd40.gwas.6
wd40.gwas.7 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 2261246 & Pos < 2262456)
wd40.gwas.7
wd40.gwas.8 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 3535366 & Pos < 3540285)
wd40.gwas.8
wd40.gwas.9 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold05", Pos > 4883313 & Pos < 4884898)
wd40.gwas.9
wd40.gwas.10 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold12", Pos > 1300977 & Pos < 1305868)
wd40.gwas.10
wd40.gwas.11 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold10", Pos > 374220 & Pos < 374936)
wd40.gwas.11
wd40.gwas.12 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold10", Pos > 383820 & Pos < 384862)
wd40.gwas.12
wd40.gwas.13 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold01", Pos > 2130910 & Pos < 2132724)
wd40.gwas.13
wd40.gwas.14 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold01", Pos > 3309095 & Pos < 3313763)
wd40.gwas.14
wd40.gwas.15 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold07", Pos > 1882961 & Pos < 1887017)
wd40.gwas.15
wd40.gwas.16 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold07", Pos > 2495402 & Pos < 2499338)
wd40.gwas.16
wd40.gwas.17 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold02", Pos > 510389 & Pos < 513008)
wd40.gwas.17
wd40.gwas.18 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold02", Pos > 2185478 & Pos < 2192354)
wd40.gwas.18
wd40.gwas.19 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold02", Pos > 2794003 & Pos < 2795777)
wd40.gwas.19
wd40.gwas.20 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold02", Pos > 2822363 & Pos < 2824510)
wd40.gwas.20
wd40.gwas.20 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold04", Pos > 2352891 & Pos < 2354877)
wd40.gwas.20
wd40.gwas.20 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold04", Pos > 4515977 & Pos < 4517169)
wd40.gwas.20

#scopa crosses noram
scopa.crosses.noram.chr <- replace(scopa.crosses.noram$Chromosome,scopa.crosses.noram$Chromosome==1, "Scaffold01")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==2, "Scaffold02")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==3, "Scaffold03")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==4, "Scaffold04")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==5, "Scaffold05")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==6, "Scaffold06")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==7, "Scaffold07")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==8, "Scaffold08")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==9, "Scaffold09")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==10, "Scaffold10")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==11, "Scaffold11")
scopa.crosses.noram.chr <- replace(scopa.crosses.noram.chr,scopa.crosses.noram.chr==12, "Scaffold12")
scopa.crosses.noram.chr
scopa.crosses.noram$Chromosome <- scopa.crosses.noram.chr
str(scopa.crosses.noram)
wd40.scopa.noram <- left_join(pfam004000.gene, scopa.crosses.noram, by = c("V1"="Chromosome")) %>% filter(Position  >= V4, Position <= V5)
str(wd40.scopa.noram)
#18160 SNPs, must investigate for specific genes:
wd40.scopa.noram.1 <- filter(scopa.crosses.noram, Chromosome=="Scaffold05" & Position > 3535366  & Position < 3540285)
wd40.scopa.noram.1 #some somewhat significant
wd40.scopa.noram.2 <- filter(scopa.crosses.noram, Chromosome=="Scaffold05" & Position > 3710714  & Position < 3711293)
wd40.scopa.noram.2 #0 SNPs
wd40.scopa.noram.3 <- filter(scopa.crosses.noram, Chromosome=="Scaffold10" & Position > 2080746  & Position < 2088064)
dim(wd40.scopa.noram.3) #224 SNPs
wd40.scopa.noram.4 <- filter(scopa.crosses.noram, Chromosome=="Scaffold01" & Position > 1970115  & Position < 1972286)
dim(wd40.scopa.noram.4) #82 SNPs

#scopa crosses subset
pfam05729.scopa.subset <- filter(scopa.crosses.subset, Chromosome==6 & Position > 1439976  & Position < 1445370)
dim(pfam05729.scopa.subset) #45 SNPs
#scopa crosses subset
scopa.crosses.subset.chr <- replace(scopa.crosses.subset$Chromosome,scopa.crosses.subset$Chromosome==1, "Scaffold01")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==2, "Scaffold02")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==3, "Scaffold03")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==4, "Scaffold04")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==5, "Scaffold05")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==6, "Scaffold06")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==7, "Scaffold07")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==8, "Scaffold08")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==9, "Scaffold09")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==10, "Scaffold10")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==11, "Scaffold11")
scopa.crosses.subset.chr <- replace(scopa.crosses.subset.chr,scopa.crosses.subset.chr==12, "Scaffold12")
scopa.crosses.subset.chr
scopa.crosses.subset$Chromosome <- scopa.crosses.subset.chr
wd40.scopa.crosses.subset <- left_join(pfam004000.gene, scopa.crosses.subset, by = c("V1"="Chromosome")) %>% filter(Position  >= V4, Position <= V5)
str(wd40.scopa.crosses.subset)
#18160 SNPs
#look closer:
#noram
wd40.scopa.crosses.subset.1 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 3535366  & Position < 3540285)
wd40.scopa.crosses.subset.1 #some quite significant
table(wd40.scopa.crosses.subset.1$P.value < 1e-5)
wd40.scopa.crosses.subset.2 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 3710714  & Position < 3711293)
wd40.scopa.crosses.subset.2 #0 SNPs
wd40.scopa.crosses.subset.3 <- filter(scopa.crosses.subset, Chromosome=="Scaffold10" & Position > 2080746  & Position < 2088064)
dim(wd40.scopa.crosses.subset.3) #224 SNPs
wd40.scopa.crosses.subset.3
table(wd40.scopa.crosses.subset.3$P.value < 1e-5)
wd40.scopa.crosses.subset.4 <- filter(scopa.crosses.subset, Chromosome=="Scaffold01" & Position > 1970115  & Position < 1972286)
dim(wd40.scopa.crosses.subset.4) #82 SNPs
wd40.scopa.crosses.subset.4
table(wd40.scopa.crosses.subset.4$P.value < 1e-5) #16

#all
wd40.scopa.1 <- filter(scopa.crosses.subset, Chromosome =="Scaffold03" & Position > 1180329 & Position < 1186358)
dim(wd40.scopa.1) #261
table(wd40.scopa.1$P.value < 1e-5)
wd40.scopa.1$P.value
wd40.scopa.2 <- filter(scopa.crosses.subset, Chromosome=="Scaffold03"& Position > 2277944 & Position < 2280931)
table(wd40.scopa.2$P.value < 1e-5)
wd40.scopa.2 #17 SNPs
wd40.scopa.3 <- filter(scopa.crosses.subset, Chromosome=="Scaffold03"& Position > 3356359 & Position < 3357989)
table(wd40.scopa.3$P.value < 1e-5)
wd40.scopa.3 #10 SNPs
wd40.scopa.4 <- filter(scopa.crosses.subset, Chromosome=="Scaffold03" & Position > 4488488 & Position < 4494271)
table(wd40.scopa.4$P.value < 1e-5)
wd40.scopa.4 #44 SNPs
wd40.scopa.5 <- filter(scopa.crosses.subset, Chromosome=="Scaffold06" & Position > 289366 & Position < 291920)
table(wd40.scopa.5$P.value < 1e-5)
wd40.scopa.5 #14 SNPs
wd40.scopa.6 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 1895702 & Position < 1897213)
table(wd40.scopa.6$P.value < 1e-5)
wd40.scopa.6
wd40.scopa.7 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 2261246 & Position < 2262456)
table(wd40.scopa.7$P.value < 1e-5)
wd40.scopa.7
wd40.scopa.8 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 3535366 & Position < 3540285)
table(wd40.scopa.8$P.value < 1e-5)
wd40.scopa.8
wd40.scopa.9 <- filter(scopa.crosses.subset, Chromosome=="Scaffold05" & Position > 4883313 & Position < 4884898)
table(wd40.scopa.9$P.value < 1e-5)
wd40.scopa.9
wd40.scopa.10 <- filter(scopa.crosses.subset, Chromosome=="Scaffold12" & Position > 1300977 & Position < 1305868)
table(wd40.scopa.10$P.value < 1e-5)
wd40.scopa.10
wd40.scopa.11 <- filter(scopa.crosses.subset, Chromosome=="Scaffold10" & Position > 374220 & Position < 374936)
table(wd40.scopa.11$P.value < 1e-5)
wd40.scopa.11
wd40.scopa.12 <- filter(scopa.crosses.subset, Chromosome=="Scaffold10" & Position > 383820 & Position < 384862)
wd40.scopa.12
wd40.scopa.13 <- filter(scopa.crosses.subset, Chromosome=="Scaffold01"& Position > 2130910 & Position < 2132724)
table(wd40.scopa.13$P.value < 1e-5)
wd40.scopa.13
wd40.scopa.14 <- filter(scopa.crosses.subset, Chromosome=="Scaffold01" & Position > 3309095 & Position < 3313763)
table(wd40.scopa.14$P.value < 1e-5)
wd40.scopa.14
wd40.scopa.15 <- filter(scopa.crosses.subset, Chromosome=="Scaffold07" & Position > 1882961 & Position < 1887017)
table(wd40.scopa.15$P.value < 1e-5)
wd40.scopa.15
wd40.scopa.16 <- filter(scopa.crosses.subset, Chromosome=="Scaffold07" & Position > 2495402 & Position < 2499338)
table(wd40.scopa.16$P.value < 1e-5)
wd40.scopa.16
wd40.scopa.17 <- filter(scopa.crosses.subset, Chromosome=="Scaffold02" & Position > 510389 & Position < 513008)
wd40.scopa.17
table(wd40.scopa.17$P.value < 1e-5)
wd40.scopa.18 <- filter(scopa.crosses.subset, Chromosome=="Scaffold02" & Position > 2185478 & Position < 2192354)
table(wd40.scopa.18$P.value < 1e-5)
wd40.scopa.18
wd40.scopa.19 <- filter(scopa.crosses.subset, Chromosome=="Scaffold02" & Position > 2794003 & Position < 2795777)
table(wd40.scopa.19$P.value < 1e-5)
wd40.scopa.19
wd40.scopa.20 <- filter(scopa.crosses.subset, Chromosome=="Scaffold02" & Position > 2822363 & Position < 2824510)
table(wd40.scopa.20$P.value < 1e-5)
wd40.scopa.20
wd40.scopa.21 <- filter(scopa.crosses.subset, Chromosome=="Scaffold04" & Position > 2352891 & Position < 2354877)
table(wd40.scopa.21$P.value < 1e-5)
wd40.scopa.21
wd40.scopa.22 <- filter(scopa.crosses.subset, Chromosome=="Scaffold04" & Position > 4515977 & Position < 4517169)
table(wd40.scopa.22$P.value < 1e-5)
wd40.gwas.22

#selection:
#noram
wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.noramB) #20 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
#no significant
wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3710714  & pos < 3711293)
dim(wd40.pos.sel.noramB) 
#0 SNPs
wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold10" & pos > 2080746  & pos < 2088064)
dim(wd40.pos.sel.noramB) #56 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))
#14 snps
wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold01" & pos > 1970115  & pos < 1972286)
dim(wd40.pos.sel.noramB) #12 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))
wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.noramA) #5 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
#0 snps
wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3710714  & pos < 3711293)
dim(wd40.pos.sel.noramA) 
#0 SNPs
wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold10" & pos > 2080746  & pos < 2088064)
dim(wd40.pos.sel.noramA) #51 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
#0 snps
wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold01" & pos > 1970115  & pos < 1972286)
dim(wd40.pos.sel.noramA) #7 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
#0 snps
wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.eurasia) #30 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
#1
wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3710714  & pos < 3711293)
dim(wd40.pos.sel.eurasia) #0 SNPs
#0
wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold10" & pos > 2080746 & pos < 2088064)
dim(wd40.pos.sel.eurasia) #
#0
wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold01" & pos > 1970115  & pos < 1972286)
dim(wd40.pos.sel.eurasia) #17 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
#1
wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.euro) #18 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
#0
wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3710714  & pos < 3711293)
dim(wd40.pos.sel.euro) 
#0 SNPs
wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold10" & pos > 2080746  & pos < 2088064)
dim(wd40.pos.sel.euro) #39 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
#2
wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold01" & pos > 1970115  & pos < 1972286)
dim(wd40.pos.sel.euro) #2 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
#0

#all
#noramA
wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold03" & pos > 1180329  & pos < 1186588)
dim(wd40.pos.sel.noramA) #88 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold03" & pos > 2277944  & pos < 2280931)
dim(wd40.pos.sel.noramA) #15 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold03" & pos > 3356359  & pos < 3357989)
dim(wd40.pos.sel.noramA) #5 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold03" & pos > 4488488  & pos < 4494271)
dim(wd40.pos.sel.noramA) #35 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold06" & pos > 289366  & pos < 291920)
dim(wd40.pos.sel.noramA) #19 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 1895702  & pos < 1897213)
dim(wd40.pos.sel.noramA) #7 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 2261246  & pos < 2262456)
dim(wd40.pos.sel.noramA) #6 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.noramA) #30 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 4883313  & pos < 4884898)
dim(wd40.pos.sel.noramA) #16 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold12" & pos > 1300977  & pos < 1305868)
dim(wd40.pos.sel.noramA) #0 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold10" & pos > 374220  & pos < 374936)
dim(wd40.pos.sel.noramA) #0 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold10" & pos > 383820  & pos < 384862)
dim(wd40.pos.sel.noramA) #0 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold01" & pos > 2130910  & pos < 2132724)
dim(wd40.pos.sel.noramA) #9 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold01" & pos > 3309095  & pos < 3313763)
dim(wd40.pos.sel.noramA) #25 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 1882961  & pos < 1887017)
dim(wd40.pos.sel.noramA) #25 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold07" & pos > 2495402  & pos < 2499338)
dim(wd40.pos.sel.noramA) #29 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold02" & pos > 510389  & pos < 513008)
dim(wd40.pos.sel.noramA) #15 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold02" & pos > 2185478  & pos < 2192354)
dim(wd40.pos.sel.noramA) #33 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold02" & pos > 2794003  & pos < 2795777)
dim(wd40.pos.sel.noramA) #12 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold02" & pos > 2822363  & pos < 2824510)
dim(wd40.pos.sel.noramA) #14 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold04" & pos > 2352891  & pos < 2354877)
dim(wd40.pos.sel.noramA) #3 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold04" & pos > 4515977  & pos < 4517169)
dim(wd40.pos.sel.noramA) #0 SNPs
wd40.pos.sel.noramA$pvalue <- 10^wd40.pos.sel.noramA$when_mutation_has_freq2 
table(wd40.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))

#noramB
wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold03" & pos > 1180329  & pos < 1186588)
dim(wd40.pos.sel.noramB) #88 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold03" & pos > 2277944  & pos < 2280931)
dim(wd40.pos.sel.noramB) #15 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold03" & pos > 3356359  & pos < 3357989)
dim(wd40.pos.sel.noramB) #5 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold03" & pos > 4488488  & pos < 4494271)
dim(wd40.pos.sel.noramB) #35 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold06" & pos > 289366  & pos < 291920)
dim(wd40.pos.sel.noramB) #19 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 1895702  & pos < 1897213)
dim(wd40.pos.sel.noramB) #7 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 2261246  & pos < 2262456)
dim(wd40.pos.sel.noramB) #6 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.noramB) #30 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 4883313  & pos < 4884898)
dim(wd40.pos.sel.noramB) #16 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold12" & pos > 1300977  & pos < 1305868)
dim(wd40.pos.sel.noramB) #0 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold10" & pos > 374220  & pos < 374936)
dim(wd40.pos.sel.noramB) #0 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold10" & pos > 383820  & pos < 384862)
dim(wd40.pos.sel.noramB) #0 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold01" & pos > 2130910  & pos < 2132724)
dim(wd40.pos.sel.noramB) #9 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold01" & pos > 3309095  & pos < 3313763)
dim(wd40.pos.sel.noramB) #25 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 1882961  & pos < 1887017)
dim(wd40.pos.sel.noramB) #25 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold07" & pos > 2495402  & pos < 2499338)
dim(wd40.pos.sel.noramB) #29 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold02" & pos > 510389  & pos < 513008)
dim(wd40.pos.sel.noramB) #15 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold02" & pos > 2185478  & pos < 2192354)
dim(wd40.pos.sel.noramB) #33 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold02" & pos > 2794003  & pos < 2795777)
dim(wd40.pos.sel.noramB) #12 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold02" & pos > 2822363  & pos < 2824510)
dim(wd40.pos.sel.noramB) #14 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold04" & pos > 2352891  & pos < 2354877)
dim(wd40.pos.sel.noramB) #3 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold04" & pos > 4515977  & pos < 4517169)
dim(wd40.pos.sel.noramB) #0 SNPs
wd40.pos.sel.noramB$pvalue <- 10^wd40.pos.sel.noramB$when_mutation_has_freq2 
table(wd40.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))

#eurasia
wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold03" & pos > 1180329  & pos < 1186588)
dim(wd40.pos.sel.eurasia) #88 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold03" & pos > 2277944  & pos < 2280931)
dim(wd40.pos.sel.eurasia) #15 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold03" & pos > 3356359  & pos < 3357989)
dim(wd40.pos.sel.eurasia) #5 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold03" & pos > 4488488  & pos < 4494271)
dim(wd40.pos.sel.eurasia) #35 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold06" & pos > 289366  & pos < 291920)
dim(wd40.pos.sel.eurasia) #19 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 1895702  & pos < 1897213)
dim(wd40.pos.sel.eurasia) #7 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 2261246  & pos < 2262456)
dim(wd40.pos.sel.eurasia) #6 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.eurasia) #30 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 4883313  & pos < 4884898)
dim(wd40.pos.sel.eurasia) #16 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold12" & pos > 1300977  & pos < 1305868)
dim(wd40.pos.sel.eurasia) #0 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold10" & pos > 374220  & pos < 374936)
dim(wd40.pos.sel.eurasia) #0 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold10" & pos > 383820  & pos < 384862)
dim(wd40.pos.sel.eurasia) #0 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold01" & pos > 2130910  & pos < 2132724)
dim(wd40.pos.sel.eurasia) #9 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold01" & pos > 3309095  & pos < 3313763)
dim(wd40.pos.sel.eurasia) #25 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eursia <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 1882961  & pos < 1887017)
dim(wd40.pos.sel.eurasia) #25 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold07" & pos > 2495402  & pos < 2499338)
dim(wd40.pos.sel.eurasia) #29 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold02" & pos > 510389  & pos < 513008)
dim(wd40.pos.sel.eurasia) #15 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold02" & pos > 2185478  & pos < 2192354)
dim(wd40.pos.sel.eurasia) #33 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold02" & pos > 2794003  & pos < 2795777)
dim(wd40.pos.sel.eurasia) #12 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold02" & pos > 2822363  & pos < 2824510)
dim(wd40.pos.sel.eurasia) #14 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold04" & pos > 2352891  & pos < 2354877)
dim(wd40.pos.sel.eurasia) #3 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold04" & pos > 4515977  & pos < 4517169)
dim(wd40.pos.sel.eurasia) #0 SNPs
wd40.pos.sel.eurasia$pvalue <- 10^wd40.pos.sel.eurasia$when_mutation_has_freq2 
table(wd40.pos.sel.eurasia$pvalue < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

#euro
wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold03" & pos > 1180329  & pos < 1186588)
dim(wd40.pos.sel.euro) #39 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold03" & pos > 2277944  & pos < 2280931)
dim(wd40.pos.sel.euro) #15 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold03" & pos > 3356359  & pos < 3357989)
dim(wd40.pos.sel.euro) #5 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold03" & pos > 4488488  & pos < 4494271)
dim(wd40.pos.sel.euro) #35 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold06" & pos > 289366  & pos < 291920)
dim(wd40.pos.sel.euro) #19 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 1895702  & pos < 1897213)
dim(wd40.pos.sel.euro) #5 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 2261246  & pos < 2262456)
dim(wd40.pos.sel.euro) #6 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 3535366  & pos < 3540285)
dim(wd40.pos.sel.euro) #18 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 4883313  & pos < 4884898)
dim(wd40.pos.sel.euro) #16 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold12" & pos > 1300977  & pos < 1305868)
dim(wd40.pos.sel.euro) #16 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold10" & pos > 374220  & pos < 374936)
dim(wd40.pos.sel.euro) #0 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold10" & pos > 383820  & pos < 384862)
dim(wd40.pos.sel.euro) #0 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold01" & pos > 2130910  & pos < 2132724)
dim(wd40.pos.sel.euro) #9 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold01" & pos > 3309095  & pos < 3313763)
dim(wd40.pos.sel.euro) #12 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 1882961  & pos < 1887017)
dim(wd40.pos.sel.euro) #12 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold07" & pos > 2495402  & pos < 2499338)
dim(wd40.pos.sel.euro) #12 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold02" & pos > 510389  & pos < 513008)
dim(wd40.pos.sel.euro) #12 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold02" & pos > 2794003  & pos < 2795777)
dim(wd40.pos.sel.euro) #12 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold02" & pos > 2822363  & pos < 2824510)
dim(wd40.pos.sel.euro) #14 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold04" & pos > 2352891  & pos < 2354877)
dim(wd40.pos.sel.euro) #3 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

wd40.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold04" & pos > 4515977  & pos < 4517169)
dim(wd40.pos.sel.euro) #0 SNPs
wd40.pos.sel.euro$pvalue <- 10^wd40.pos.sel.euro$when_mutation_has_freq2 
table(wd40.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))

###check rcd
rcd1
dim(rcd1)
#deltafst all
rcd1.delta.fst.all <- filter(delta_fst.all, Scaffold=="Scaffold10" & Start > 896396 & End < 896396 + 11000)
rcd1.delta.fst.all #0
#deltafst noram
rcd1.delta.fst.noram <- filter(delta_fst.noram, Scaffold=="Scaffold10" & Start > 896396 & End < 896396 + 11000)
rcd1.delta.fst.noram #0
#pos.sel noramA
rcd1.pos.sel.noramA <- filter(pos.sel.noramA, rs_id=="Scaffold10" & pos > 896396  & pos < 898399)
dim(rcd1.pos.sel.noramA) #9 SNPs
rcd1.pos.sel.noramA$pvalue <- 10^rcd1.pos.sel.noramA$when_mutation_has_freq2 
rcd1.pos.sel.noramA #2 somewhat significant
rcd1.pos.sel.noramA$pvalue < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#pos sel noramB
rcd1.pos.sel.noramB <- filter(pos.sel.noramB, rs_id=="Scaffold10" & pos > 896396  & pos < 898399)
dim(rcd1.pos.sel.noramB) #9 SNPs
rcd1.pos.sel.noramB$pvalue <- 10^rcd1.pos.sel.noramB$when_mutation_has_freq2 
rcd1.pos.sel.noramB #0 significant
rcd1.pos.sel.noramB$pvalue < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#pos sel Eurasia
rcd1.pos.sel.eurasia <- filter(pos.sel.eurasia, rs_id=="Scaffold10" & pos > 896396  & pos < 898399)
dim(rcd1.pos.sel.eurasia) #0 SNPs
##pos sel euro
rcd1.pos.sel.euro <- filter(pos.sel.euro, rs_id=="Scaffold10" & pos > 896396  & pos < 898399)
dim(rcd1.pos.sel.euro) #11 SNPs
rcd1.pos.sel.euro$pvalue <- 10^rcd1.pos.sel.euro$when_mutation_has_freq2 
rcd1.pos.sel.euro #0 significant
rcd1.pos.sel.euro$pvalue < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05)
#

#gwas
rcd1.gwas <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold10", Pos > 896396  & Pos < 898399)
dim(rcd1.gwas) #0 SNPs
#scopa crosses noram
rcd1.scopa.noram <- filter(scopa.crosses.noram, Chromosome==10 & Position > 896396  & Position < 898399)
dim(rcd1.scopa.noram) #54 SNPs
rcd1.scopa.noram #1 very significant, several quite significant
#scopa crosses subset
rcd1.scopa.subset <- filter(scopa.crosses.subset, Chromosome==10 & Position > 896396  & Position < 898399)
dim(rcd1.scopa.subset) #54 SNPs
rcd1.scopa.subset #2 very significant, many quite significant
table(rcd1.scopa.subset$P.value < 1e-5)
###check caspases:
caspases
#redundant entries:
#1 and 2, 
#3 and 4, 
#5 6 7 8 and 9, 
#10 and 11,
#12, 13, and 14,
#15, 16 and 17
#18 and 19
#20 and 21, 22
#23 and 24

caspases <- caspases[c(1,3,5,10,12,15,18,20,23),]
caspases
#deltafst all:
caspases[1,]
blast1 <- filter(delta_fst.all, Scaffold=="Scaffold06" & Start > 1199504)
blast1 <- filter(blast1, End < 1199504 + 11000)
blast1 #1 hit
caspases[2,]
blast2 <- filter(delta_fst.all, Scaffold=="Scaffold05" & Start > 4750429 )
blast2 <- filter(blast2, End < 4750429  + 11000 )
blast2 #0 hits
caspases[3,]
blast3 <- filter(delta_fst.all, Scaffold=="Scaffold08" & Start > 4742004)
blast3 <- filter(blast3, End < 4742004 + 11000)
blast3 #0 hits
caspases[4,]
blast4 <- filter(delta_fst.all, Scaffold=="Scaffold08" & Start > 442633)
blast4 <- filter(blast4, End < 442633 + 11000)
blast4 #0 hits
caspases[5,]
blast5 <- filter(delta_fst.all, Scaffold=="Scaffold08" & Start > 2437877)
blast5 <- filter(blast5, End < 2437877 + 11000)
blast5 #0 rows
caspases[6,]
blast6 <- filter(delta_fst.all, Scaffold=="Scaffold08" & Start > 2460146)
blast6 <- filter(blast6, End < 2460146 + 11000)
blast6 #0 rows
caspases[7,]
blast7 <- filter(delta_fst.all, Scaffold=="Scaffold08" & Start > 2532548)
blast7 <- filter(blast7, End < 2532548 + 11000)
blast7 #0 rows
caspases[8,]
blast8 <- filter(delta_fst.all, Scaffold=="Scaffold12" & Start > 1306285)
blast8 <- filter(blast8, End < 1306285 + 11000)
blast8 #O rows
caspases[9,]
blast9 <- filter(delta_fst.all, Scaffold=="Scaffold10" & Start > 2042757)
blast9 <- filter(blast9, End < 2042757 + 11000)
blast9 #0 rows

#in deltafst noram:
caspases[1,]
blast1 <- filter(delta_fst.noram, Scaffold=="Scaffold06" & Start > 1199504)
blast1 <- filter(blast1, End < 1199504 + 11000)
blast1 #0 hit
caspases[2,]
blast2 <- filter(delta_fst.noram, Scaffold=="Scaffold05" & Start > 4750429 )
blast2 <- filter(blast2, End < 4750429  + 11000 )
blast2 #0 hits
caspases[3,]
blast3 <- filter(delta_fst.noram, Scaffold=="Scaffold08" & Start > 4742004)
blast3 <- filter(blast3, End < 4742004 + 11000)
blast3 #0 hits
caspases[4,]
blast4 <- filter(delta_fst.noram, Scaffold=="Scaffold08" & Start > 442633)
blast4 <- filter(blast4, End < 442633 + 11000)
blast4 #0 hits
caspases[5,]
blast5 <- filter(delta_fst.noram, Scaffold=="Scaffold08" & Start > 2437877)
blast5 <- filter(blast5, End < 2437877 + 11000)
blast5 #0 rows
caspases[6,]
blast6 <- filter(delta_fst.noram, Scaffold=="Scaffold08" & Start > 2460146)
blast6 <- filter(blast6, End < 2460146 + 11000)
blast6 #0 rows
caspases[7,]
blast7 <- filter(delta_fst.noram, Scaffold=="Scaffold08" & Start > 2532548)
blast7 <- filter(blast7, End < 2532548 + 11000)
blast7 #0 rows
caspases[8,]
blast8 <- filter(delta_fst.noram, Scaffold=="Scaffold12" & Start > 1306285)
blast8 <- filter(blast8, End < 1306285 + 11000)
blast8 #O rows
caspases[9,]
blast9 <- filter(delta_fst.noram, Scaffold=="Scaffold10" & Start > 2042757)
blast9 <- filter(blast9, End < 2042757 + 11000)
blast9 #0 rows

#positive selection noramA
str(pos.sel.noramA)
caspases[1,]
blast1 <- filter(pos.sel.noramA, rs_id=="Scaffold06" & pos > 1199504 & pos < 1207398)
dim(blast1) #39 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
sort(blast1$p.value.freq2) #4 somewhat significant
blast1$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
table(blast1$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
caspases[2,]
blast2 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 4750429 & pos < 4755890)
dim(blast2) #0 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2
caspases[3,]
blast3 <- filter(pos.sel.noramA, rs_id=="Scaffold05" & pos > 4742004 & pos < 4745188)
dim(blast3) #0 snps
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3
caspases[4,]
blast4 <- filter(pos.sel.noramA, rs_id=="Scaffold08" & pos > 442633 & pos < 445273)
dim(blast4) #5 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #1 somewhat significant
blast4$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
caspases[5,]
blast5 <- filter(pos.sel.noramA, rs_id=="Scaffold08" & pos > 2437877 & pos < 2444353)
dim(blast5) #29 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5 #2 somewhat significant
blast5$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
caspases[6,]
blast6 <- filter(pos.sel.noramA, rs_id=="Scaffold08" & pos > 2460146 & pos < 2474674)
dim(blast6) #68 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
sort(blast6$p.value.freq2) #many somewhat significant
blast6$p.value.freq2 < 1e-4
table(blast6$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
caspases[7,]
blast7 <- filter(pos.sel.noramA, rs_id=="Scaffold08" & pos > 2532548  & pos < 2537307)
dim(blast7) #12 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7 #3 somewhat significant
table(blast7$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05))
caspases[8,]
blast8 <- filter(pos.sel.noramA, rs_id=="Scaffold12" & pos > 1306285 & pos < 1309400)
dim(blast8) #8 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8
blast8$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
caspases[9,]
blast9 <- filter(pos.sel.noramA, rs_id=="Scaffold10" & pos > 2042757 & pos < 2057519)
dim(blast9) #38 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
sort(blast9$p.value.freq2) #several somewhat significant
blast9$p.value.freq2 < quantile(10^pos.sel.noramA.filtered$when_mutation_has_freq2, probs=0.05)
#Noram A:1,4,5,6,7,9 

#positive selection noramB
str(pos.sel.noramB)
caspases[1,]
blast1 <- filter(pos.sel.noramB, rs_id=="Scaffold06" & pos > 1199504 & pos < 1207398)
dim(blast1) #82 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
sort(blast1$p.value.freq2) #many somewhat significant
table(blast1$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))
caspases[2,]
blast2 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 4750429 & pos < 4755890)
dim(blast2) #0 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2
caspases[3,]
blast3 <- filter(pos.sel.noramB, rs_id=="Scaffold05" & pos > 4742004 & pos < 4745188)
dim(blast3) #0 snps
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3
caspases[4,]
blast4 <- filter(pos.sel.noramB, rs_id=="Scaffold08" & pos > 442633 & pos < 445273)
dim(blast4) #8 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #4 somewhat significant
blast4$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
caspases[5,]
blast5 <- filter(pos.sel.noramB, rs_id=="Scaffold08" & pos > 2437877 & pos < 2444353)
dim(blast5) #23 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5 #several somewhat significant
table(blast5$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))
caspases[6,]
blast6 <- filter(pos.sel.noramB, rs_id=="Scaffold08" & pos > 2460146 & pos < 2474674)
dim(blast6) #55 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
sort(blast6$p.value.freq2) #many somewhat significant
table(blast6$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05))
caspases[7,]
blast7 <- filter(pos.sel.noramB, rs_id=="Scaffold08" & pos > 2532548  & pos < 2537307)
dim(blast7) #12 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7 #1 somewhat significant
blast7$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
caspases[8,]
blast8 <- filter(pos.sel.noramB, rs_id=="Scaffold12" & pos > 1306285 & pos < 1309400)
dim(blast8) #18 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8 #6 somewhat significant
blast8$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
caspases[9,]
blast9 <- filter(pos.sel.noramB, rs_id=="Scaffold10" & pos > 2042757 & pos < 2057519)
dim(blast9) #36 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
sort(blast9$p.value.freq2) #several somewhat significant
blast9$p.value.freq2 < quantile(10^pos.sel.noramB.filtered$when_mutation_has_freq2, probs=0.05)
#Noram A:1,4,5,6,7,9
#Noram B: 1,4,5,6,7,8,9

#positive selection eurasia
str(pos.sel.eurasia)
caspases[1,]
blast1 <- filter(pos.sel.eurasia, rs_id=="Scaffold06" & pos > 1199504 & pos < 1207398)
dim(blast1) #89 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
sort(blast1$p.value.freq2) #many somewhat significant
table(blast1$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))

caspases[2,]
blast2 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 4750429 & pos < 4755890)
dim(blast2) #0 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2
caspases[3,]
blast3 <- filter(pos.sel.eurasia, rs_id=="Scaffold05" & pos > 4742004 & pos < 4745188)
dim(blast3) #0 snps
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3
caspases[4,]
blast4 <- filter(pos.sel.eurasia, rs_id=="Scaffold08" & pos > 442633 & pos < 445273)
dim(blast4) #26 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #1 somewhat significant
table(blast4$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
caspases[5,]
blast5 <- filter(pos.sel.eurasia, rs_id=="Scaffold08" & pos > 2437877 & pos < 2444353)
dim(blast5) #27 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5 #4 somewhat significant
table(blast5$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
caspases[6,]
blast6 <- filter(pos.sel.eurasia, rs_id=="Scaffold08" & pos > 2460146 & pos < 2474674)
dim(blast6) #63 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
sort(blast6$p.value.freq2) #many somewhat significant
table(blast6$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
caspases[7,]
blast7 <- filter(pos.sel.eurasia, rs_id=="Scaffold08" & pos > 2532548  & pos < 2537307)
dim(blast7) #11 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7 #2 somewhat significant
table(blast7$p.value.freq2 < quantile(10^pos.sel.eurasia.filtered$when_mutation_has_freq2, probs=0.05))
caspases[8,]
blast8 <- filter(pos.sel.eurasia, rs_id=="Scaffold12" & pos > 1306285 & pos < 1309400)
dim(blast8) #0 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8 
caspases[9,]
blast9 <- filter(pos.sel.eurasia, rs_id=="Scaffold10" & pos > 2042757 & pos < 2057519)
dim(blast9) #0 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
sort(blast9$p.value.freq2) 
#Noram A:1,4,5,6,7,9
#Noram B: 1,4,5,6,7,8,9
#Eurasia: 1,4,5,6,7

#positive selection Europe
str(pos.sel.euro)
caspases[1,]
blast1 <- filter(pos.sel.euro, rs_id=="Scaffold06" & pos > 1199504 & pos < 1207398)
dim(blast1) #22 SNPs
blast1$p.value.freq2 <- 10^blast1$when_mutation_has_freq2 
sort(blast1$p.value.freq2) #no significant
table(blast1$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[2,]
blast2 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 4750429 & pos < 4755890)
dim(blast2) #0 snps
blast2$p.value.freq2 <- 10^blast2$when_mutation_has_freq2 
blast2
caspases[3,]
blast3 <- filter(pos.sel.euro, rs_id=="Scaffold05" & pos > 4742004 & pos < 4745188)
dim(blast3) #0 snps
blast3$p.value.freq2 <- 10^blast3$when_mutation_has_freq2 
blast3
caspases[4,]
blast4 <- filter(pos.sel.euro, rs_id=="Scaffold08" & pos > 442633 & pos < 445273)
dim(blast4) #7 SNPs
blast4$p.value.freq2 <- 10^blast4$when_mutation_has_freq2 
blast4 #0 significant
table(blast4$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[5,]
blast5 <- filter(pos.sel.euro, rs_id=="Scaffold08" & pos > 2437877 & pos < 2444353)
dim(blast5) #10 SNPs
blast5$p.value.freq2 <- 10^blast5$when_mutation_has_freq2 
blast5 #0 significant
table(blast5$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[6,]
blast6 <- filter(pos.sel.euro, rs_id=="Scaffold08" & pos > 2460146 & pos < 2474674)
dim(blast6) #38 SNPs
blast6$p.value.freq2 <- 10^blast6$when_mutation_has_freq2 
sort(blast6$p.value.freq2) #no significant
table(blast6$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[7,]
blast7 <- filter(pos.sel.euro, rs_id=="Scaffold08" & pos > 2532548  & pos < 2537307)
dim(blast7) #7 SNPs
blast7$p.value.freq2 <- 10^blast7$when_mutation_has_freq2 
blast7 
table(blast7$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[8,]
blast8 <- filter(pos.sel.euro, rs_id=="Scaffold12" & pos > 1306285 & pos < 1309400)
dim(blast8) #8 SNPs
blast8$p.value.freq2 <- 10^blast8$when_mutation_has_freq2 
blast8 
table(blast8$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
caspases[9,]
blast9 <- filter(pos.sel.euro, rs_id=="Scaffold10" & pos > 2042757 & pos < 2057519)
dim(blast9) #22 SNPs
blast9$p.value.freq2 <- 10^blast9$when_mutation_has_freq2 
sort(blast9$p.value.freq2) 
table(blast9$p.value.freq2 < quantile(10^pos.sel.euro.filtered$when_mutation_has_freq2, probs=0.05))
#

#gwas
caspases
caspases[1,]
blast1 <- filter(gwas.all.inferred, gwas.all.inferred$Chr=="Scaffold06", Pos > 1199504 & Pos < 1207398)
dim(blast1) #8 SNPs
blast1 #all quite significant
caspases[2,]
blast2 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 4750429 & Pos < 4755890)
dim(blast2) #0 snps
caspases[3,]
blast3 <- filter(gwas.all.inferred, Chr=="Scaffold05" & Pos > 4742004 & Pos < 4745188)
dim(blast3) #0 SNPs
caspases[4,]
blast4 <- filter(gwas.all.inferred, Chr=="Scaffold08" & Pos > 442633 & Pos < 445273)
dim(blast4) #0 SNPs
caspases[5,]
blast5 <- filter(gwas.all.inferred, Chr=="Scaffold08" & Pos > 2437877 & Pos < 2444353)
dim(blast5) #0 SNPs
caspases[6,]
blast6 <- filter(gwas.all.inferred, Chr=="Scaffold08" & Pos > 2460146 & Pos < 2474674)
dim(blast6) #1 SNPs
blast6 #quite significant
caspases[7,]
blast7 <- filter(gwas.all.inferred, Chr=="Scaffold08" & Pos > 2532548  & Pos < 2537307)
dim(blast7) #0 SNPs
caspases[8,]
blast8 <- filter(gwas.all.inferred, Chr=="Scaffold12" & Pos > 1306285 & Pos < 1309400)
dim(blast8) #0 SNPs
caspases[9,]
blast9 <- filter(gwas.all.inferred, Chr=="Scaffold10" & Pos > 2042757 & Pos < 2057519 )
dim(blast9) #0 SNPs
blast9 #quite significant

#scopa noram
caspases[1,]
blast1 <- filter(scopa.crosses.noram, Chromosome== 6, Position > 1199504 & Position < 1207398 )
dim(blast1) #265 SNPs
blast1 #several with somewhat low p-value
caspases[2,]
blast2 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 4750429 & Position < 4755890)
dim(blast2) #0 snps
caspases[3,]
blast3 <- filter(scopa.crosses.noram, Chromosome==5 & Position > 4742004 & Position < 4745188)
dim(blast3) #0 SNPs
caspases[4,]
blast4 <- filter(scopa.crosses.noram, Chromosome==8 & Position > 442633 & Position < 445273)
dim(blast4) #137 SNPs
blast4 #1 somewhat low
caspases[5,]
blast5 <- filter(scopa.crosses.noram, Chromosome==8 & Position > 2437877 & Position < 2444353)
dim(blast5) #127 SNPs
blast5 #several quite low
caspases[6,]
blast6 <- filter(scopa.crosses.noram, Chromosome==8 & Position > 2460146 & Position < 2474674)
dim(blast6) #309 SNPs
blast6 #some quite low
caspases[7,]
blast7 <- filter(scopa.crosses.noram, Chromosome==8 & Position > 2532548 & Position < 2537307)
dim(blast7) #46 SNPs
blast7 #a few somewhat low
caspases[8,]
blast8 <- filter(scopa.crosses.noram, Chromosome==12 & Position > 1306285 & Position < 1309400)
dim(blast8) #69 SNPs
blast8 #one low
caspases[9,]
blast9 <- filter(scopa.crosses.noram, Chromosome==10 & Position > 2042757 & Position < 2057519)
dim(blast9) #153 SNPs
blast9 #a few low 

#scopa subset
caspases[1,]
blast1 <- filter(scopa.crosses.subset, Chromosome== 6, Position > 1199504 & Position < 1207398 )
dim(blast1) #265 SNPs
blast1 #many with low p-value
table(blast1$P.value < 1e-05) #33
caspases[2,]
blast2 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 4750429 & Position < 4755890)
dim(blast2) #0 snps
caspases[3,]
blast3 <- filter(scopa.crosses.subset, Chromosome==5 & Position > 4742004 & Position < 4745188)
dim(blast3) #0 SNPs
caspases[4,]
blast4 <- filter(scopa.crosses.subset, Chromosome==8 & Position > 442633 & Position < 445273)
dim(blast4) #137 SNPs
blast4 #many extremely low
table(blast4$P.value < 1e-05) #27
caspases[5,]
blast5 <- filter(scopa.crosses.subset, Chromosome==8 & Position > 2437877 & Position < 2444353)
dim(blast5) #127 SNPs
blast5 #many quite low
table(blast5$P.value < 1e-05) #14
caspases[6,]
blast6 <- filter(scopa.crosses.subset, Chromosome==8 & Position > 2460146 & Position < 2474674)
dim(blast6) #309 SNPs
blast6 #some quite low
table(blast6$P.value < 1e-05) 
caspases[7,]
blast7 <- filter(scopa.crosses.subset, Chromosome==8 & Position > 2532548 & Position < 2537307)
dim(blast7) #46 SNPs
blast7 #some quite low
table(blast7$P.value < 1e-05) 
caspases[8,]
blast8 <- filter(scopa.crosses.subset, Chromosome==12 & Position > 1306285 & Position < 1309400)
dim(blast8) #69 SNPs
blast8 #some very low
table(blast8$P.value < 1e-05) 
caspases[9,]
blast9 <- filter(scopa.crosses.subset, Chromosome==10 & Position > 2042757 & Position < 2057519)
dim(blast9) #153 SNPs
blast9 #some very low 
table(blast9$P.value < 1e-05) 
