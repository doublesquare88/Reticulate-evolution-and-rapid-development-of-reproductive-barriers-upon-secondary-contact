rm(list = ls())
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggvenn")
library("ggvenn")
install.packages("ggVennDiagram")
#library("ggVennDiagram")
.libPaths()

setwd("/Users/dabaosl/Dropbox/UiO/Phd/Data_analysis")
setwd("/Users/dabaosl/UiO Dropbox/Dabao Lu/UiO/Phd/Data_analysis")
pop <- read.table("pop_info_sub.txt")
#pop <- read.table("pop_info_superfine.txt")

levels(as.factor(pop$V2))

table(pop$V2)

# Prepare input files:
# Read in the file with sliding window estimates of FST, pi and dxy
setwd("/Users/dabaosl/Dropbox/UiO/Phd/Data_analysis/Genome_scans")
setwd("/Users/dabaosl/UiO Dropbox/Dabao Lu/UiO/Phd/Data_analysis/Genome_scans")

windowStats<-read.csv("no_scaffold9/Fst_Dxy_pi.csv",header=T)
windowStats<-read.csv("Fst_Dxy_pi_superfine.csv",header=T)
windowStats<-read.csv("Fst_Dxy_pi_superfine_10kb_5kb2.csv",header=T)
windowStats<-read.csv("Fst_Dxy_pi.csv", header=T)
windowStats<-read.csv("Fst_Dxy_pi_main_EuroSS_10kb_5kb.csv", header=T)

# Let's have a look at what columns are present
names(windowStats)
length(windowStats$scaffold)
head(windowStats)

# Get chrom ends for making positions additive
chrom<-read.table("chrEnds.txt",header=T)
chrom
chrom$add<-c(0,cumsum(chrom$end)[1:11])
chrom
windowStats$scaffold
chrom$scaffold

levels(as.factor(windowStats$scaffold))
# Make the positions of the divergence and diversity window estimates additive
windowStats$mid <-windowStats$mid+chrom[match(windowStats$scaffold,chrom$scaffold),3]
windowStats$mid

###########################################
#                                         #
#   get summary statistics                #
#                                         #
###########################################

names(windowStats)

mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro_N, na.rm = T)
mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro_S, na.rm = T)
mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm, na.rm = T)
mean(windowStats$Fst_Circumboreal_NorAm_North_American, na.rm = T)
mean(windowStats$Fst_Circumboreal_Asia_North_American, na.rm = T)
mean(windowStats$Fst_Circumboreal_Euro_S_North_American, na.rm = T)
mean(windowStats$Fst_Circumboreal_Euro_N_North_American, na.rm = T)
mean(windowStats$Fst_Circumboreal_Euro_S_Circumboreal_NorAm, na.rm = T)
mean(windowStats$Fst_Circumboreal_Euro_N_Circumboreal_NorAm, na.rm = T)
mean(windowStats$Fst_Circumboreal_Euro_S_Circumboreal_Euro_N, na.rm = T)

mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro_N, na.rm = T) #0.0819963
mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro_S, na.rm = T) #0.09624772
mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_NorAm, na.rm = T) #0.07805347
mean(windowStats$dxy_Circumboreal_NorAm_North_American, na.rm = T) #0.1012712
mean(windowStats$dxy_Circumboreal_Asia_North_American, na.rm = T) #0.1208102
mean(windowStats$dxy_Circumboreal_Euro_S_North_American, na.rm = T) #0.1098921
mean(windowStats$dxy_Circumboreal_Euro_N_North_American, na.rm = T) # 0.107001
mean(windowStats$dxy_Circumboreal_Euro_S_Circumboreal_NorAm, na.rm = T) #0.08566079
mean(windowStats$dxy_Circumboreal_Euro_N_Circumboreal_NorAm, na.rm = T) #0.066813
mean(windowStats$dxy_Circumboreal_Euro_S_Circumboreal_Euro_N, na.rm = T) #0.06230157

#with correct pop partition for fig2 in paper:
mean(windowStats$dxy_Circumboreal_Euro_N_Circumboreal_Asia, na.rm = T) #0.08018874
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_Asia, na.rm = T) #0.10282
mean(windowStats$dxy_Circumboreal_NorAm_Circumboreal_Asia, na.rm = T) #0.07585359
mean(windowStats$dxy_Circumboreal_NorAm_North_American, na.rm = T) #0.09604966
mean(windowStats$dxy_Circumboreal_Asia_North_American, na.rm = T) #0.1148088
mean(windowStats$dxy_Circumboreal_Euro_SS_North_American, na.rm = T) #0.1068726
mean(windowStats$dxy_Circumboreal_Euro_N_North_American, na.rm = T) #0.1015187
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_NorAm, na.rm = T) #0.0935133
mean(windowStats$dxy_Circumboreal_Euro_N_Circumboreal_NorAm, na.rm = T) #0.06505847
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_Euro_N, na.rm = T) #0.07651295
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_Euro_S, na.rm = T) #0.06108651
#add with noramC:
mean(windowStats$dxy_Circumboreal_Euro_SS_North_American_2, na.rm = T) #0.1032746
mean(windowStats$dxy_Circumboreal_Asia_North_American_2, na.rm = T) #0.1122938

mean(windowStats$dxy_Circumboreal_NorAm_North_American_2, na.rm = T) #0.09703957
mean(windowStats$dxy_Circumboreal_Euro_N_North_American_2, na.rm = T) #0.0978046
mean(windowStats$dxy_Circumboreal_NorAm_North_American_2, na.rm = T) #0.09703957

#get for asian sublineages only:
mean(windowStats$dxy_Asia_N_Asia_SW, na.rm = T) #0.06737796
mean(windowStats$dxy_Asia_NE_Asia_Sichuan, na.rm = T) #0.06745493
mean(windowStats$dxy_Asia_N_Asia_Sichuan, na.rm = T) #0.06371048
mean(windowStats$dxy_Asia_NE_Asia_N, na.rm = T) #0.06346031
mean(windowStats$dxy_Asia_NE_Asia_SW, na.rm = T) #0.05096864
mean(windowStats$dxy_Asia_SW_Asia_Sichuan, na.rm = T) #0.06675611
mean(windowStats$dxy_Circumboreal_Euro_N_North_American, na.rm = T) #0.1015187
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_NorAm, na.rm = T) #0.0935133
mean(windowStats$dxy_Circumboreal_Euro_N_Circumboreal_NorAm, na.rm = T) #0.06505847
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_Euro_N, na.rm = T) #0.07651295
mean(windowStats$dxy_Circumboreal_Euro_SS_Circumboreal_Euro_S, na.rm = T) #0.06108651


##########################################################################
# Plotting statistics along the genome:

#Set color palette:
help(palette)
palette.pals()
#palette("Dark 2")
palette("Paired")
#palette("Okabe-Ito")
#palette("Classic Tableau")
palette(c("grey40","grey60"))

#######################################
#                                     #
#   Superfine_pop_designation         #
#                                     #
#######################################

###################
#                 #  
#     Fst         #  
#                 #
###################

#pdf("Fst_windows.pdf", height=12, width=12)
par(mfrow=c(5,1),oma=c(3,0,0,0),mar=c(4,5,2,2))


par(mfrow=c(1,1))
### Asia  within ###
plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Asia_N",
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_N,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_SW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst EastAsia_NE x EastAsia_SW", 
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_SW,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Asia_Sichuan",
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Asia_N_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_N_Asia_Sichuan",
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Asia_SW_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_SW_Asia_Sichuan",
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Asia_N_Asia_SW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_N_Asia_SW",
     ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Asia_SW,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
#pdf("assess_mating.pdf")

names(windowStats)
###############################################################
#                                                             #
# Comparisons between groups that don´t seem to mate well     #
#                                                             #
###############################################################

##Delta
Delta_df <- as.data.frame(cbind(windowStats$scaffold,windowStats$start,windowStats$end, windowStats$mid, windowStats$sites))
colnames(Delta_df)[1] <-"Scaffold"
colnames(Delta_df)[2] <-"Start"
colnames(Delta_df)[3] <-"End"
colnames(Delta_df)[4] <-"Mid"
colnames(Delta_df)[5] <- "Sites"

palette(c("grey60", "firebrick2"))

#ASIA
par(mfrow=c(1,1))
Delta_df$D.Asia_SW.Sichuan_SW.NE <- windowStats$Fst_Asia_SW_Asia_Sichuan - windowStats$Fst_Asia_NE_Asia_SW
# identify the 95% and 99% percentile
quantile(Delta_df$D.Asia_SW.Sichuan_SW.NE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.Asia_SW.Sichuan_SW.NE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.Asia_SW.Sichuan_SW.NE = ifelse(Delta_df$D.Asia_SW.Sichuan_SW.NE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.Asia_SW.Sichuan_SW.NE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.Asia_SW.Sichuan_SW.NE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
      main="DeltaFst EastAsiaSW_EastAsiaSC - EastAsiaNE_EastAsiaSW",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.Asia_SW.Sichuan_SW.NE))
abline(h=mean(Delta_df$D.Asia_SW.Sichuan_SW.NE,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.Asia_SW.Sichuan_Sichuan.NE <- windowStats$Fst_Asia_SW_Asia_Sichuan - windowStats$Fst_Asia_NE_Asia_Sichuan
# identify the 95% and 99% percentile
quantile(Delta_df$D.Asia_SW.Sichuan_Sichuan.NE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.Asia_SW.Sichuan_SW.NE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.Asia_SW.Sichuan_Sichuan.NE = ifelse(Delta_df$D.Asia_SW.Sichuan_Sichuan.NE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.Asia_SW.Sichuan_Sichuan.NE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.Asia_SW.Sichuan_Sichuan.NE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst Asia_SW_Asia_Sichuan - Asia_Sichuan_Asia_NE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.Asia_SW.Sichuan_Sichuan.NE))
abline(h=mean(Delta_df$D.Asia_SW.Sichuan_Sichuan.NE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.Asia_SW.N_N.NE <- windowStats$Fst_Asia_N_Asia_SW - windowStats$Fst_Asia_NE_Asia_N
# identify the 95% and 99% percentile
quantile(Delta_df$D.Asia_SW.N_N.NE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.Asia_SW.N_N.NE, 0.975, na.rm = T)
# make an outlier column in the data.frame

par(mfrow=c(1,1))
Delta_df <- Delta_df %>% mutate(outlier.D.Asia_SW.N_N.NE = ifelse(Delta_df$D.Asia_SW.N_N.NE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.Asia_SW.N_N.NE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.Asia_SW.N_N.NE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst Asia_SW_Asia_N - Asia_N_Asia_NE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.Asia_SW.N_N.NE))
abline(h=mean(Delta_df$D.Asia_SW.N_N.NE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.Asia_SW.N_N.Sichuan <- windowStats$Fst_Asia_N_Asia_SW - windowStats$Fst_Asia_N_Asia_Sichuan
# identify the 95% and 99% percentile
quantile(Delta_df$D.Asia_SW.N_N.Sichuan, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.Asia_SW.N_N.Sichuan, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.Asia_SW.N_N.Sichuan = ifelse(Delta_df$D.Asia_SW.N_N.Sichuan > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.Asia_SW.N_N.Sichuan) %>% tally()
plot(Delta_df$Mid,Delta_df$D.Asia_SW.N_N.Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst Asia_SW_Asia_N - Asia_N_Asia_Sichuan",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.Asia_SW.N_N.Sichuan))
abline(h=mean(Delta_df$D.Asia_SW.N_N.Sichuan,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df_asia_within <- Delta_df

#North America
Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2 <- windowStats$Fst_Circumboreal_NorAm_E_North_American - windowStats$Fst_Circumboreal_NorAm_E_North_American_2
# identify the 95% and 99% percentile
quantile(Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.Noram_NoramE.CircNoramE_CircNoramE.Noram2 = ifelse(Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2 > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.Noram_NoramE.CircNoramE_CircNoramE.Noram2) %>% tally()
plot(Delta_df$Mid,Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst Noram_NoramE - CircNoramE_CircNoramE.Noram2",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.Noram_NoramE.CircNoramE_CircNoramE.Noram2))
abline(h=mean(Delta_df$D.Noram_NoramE.CircNoramE_CircNoramE.Noram2,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW <- windowStats$Fst_Circumboreal_NorAm_E_North_American - windowStats$Fst_Circumboreal_NorAm_W_North_American_S
# identify the 95% and 99% percentile
quantile(Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.NoramAE_NoramBE.NoramAS_NoramBW = ifelse(Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.NoramAE_NoramBE.NoramAS_NoramBW) %>% tally()
plot(Delta_df$Mid,Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW ,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="DeltaFst NoramAE_NoramBE - NoramAS_NoramBW",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.NoramAE_NoramBE.NoramAS_NoramBW))
abline(h=mean(Delta_df$D.NoramAE_NoramBE.NoramAS_NoramBW,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df_Noram_within <- Delta_df


#Asia X North America
Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE <- windowStats$Fst_Asia_SW_North_American - windowStats$Fst_Asia_NE_North_American
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.NoramE_AsiaNE.NoramE = ifelse(Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.NoramE_AsiaNE.NoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_NoramE - AsiaNE_NoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE))
abline(h=mean(Delta_df$D.AsiaSW.NoramE_AsiaNE.NoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE <- windowStats$Fst_Asia_SW_North_American - windowStats$Fst_Asia_N_North_American
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.NoramE_AsiaN.NoramE = ifelse(Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.NoramE_AsiaN.NoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_NoramE - AsiaN_NoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.NoramE_AsiaN.NoramE))
abline(h=mean(Delta_df$D.AsiaSW.NoramE_AsiaN.NoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE <- windowStats$Fst_Asia_Sichuan_North_American - windowStats$Fst_Asia_NE_North_American
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE = ifelse(Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSichuan_NoramE - AsiaNE_NoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE))
abline(h=mean(Delta_df$D.AsiaSichuan.NoramE_AsiaNE.NoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE <- windowStats$Fst_Asia_Sichuan_North_American - windowStats$Fst_Asia_N_North_American
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE = ifelse(Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSichuan_NoramE - AsiaN_NoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE))
abline(h=mean(Delta_df$D.AsiaSichuan.NoramE_AsiaN.NoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

#Asia X Circumboreal North AmericaE
Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE <- windowStats$Fst_Asia_SW_Circumboreal_NorAm_E - windowStats$Fst_Asia_NE_Circumboreal_NorAm_E
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE, 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE = ifelse(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_CircNoramE - AsiaNE_CircNoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE))
abline(h=mean(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE <- windowStats$Fst_Asia_SW_Circumboreal_NorAm_E - windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_E
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE = ifelse(Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_CircNoramE - AsiaSichuan_CircNoramE",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE))
abline(h=mean(D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

#Asia X Circumboreal North AmericaW
Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW <- windowStats$Fst_Asia_SW_Circumboreal_NorAm_W - windowStats$Fst_Asia_NE_Circumboreal_NorAm_W
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.CircNoramW_AsiaNE.CircNoramW = ifelse(Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.CircNoramW_AsiaNE.CircNoramW) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_CircNoramW - AsiaNE_CircNoramW",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaNE.CircNoramW))
abline(h=mean(Delta_df$D.AsiaSW.CircNoramW_AsiaNE.CircNoramW,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW <- windowStats$Fst_Asia_SW_Circumboreal_NorAm_W - windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_W
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW = ifelse(Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_CircNoramW - AsiaSicuan_CircNoramW",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW))
abline(h=mean(Delta_df$D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

#Asia EuroN
Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN <- windowStats$Fst_Asia_SW_Euro_N - windowStats$Fst_Asia_NE_Euro_N
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.EuroN_AsiaNE.EuroN = ifelse(Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.EuroN_AsiaNE.EuroN) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_EuroN - AsiaNE_EuroN",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN))
abline(h=mean(Delta_df$D.AsiaSW.EuroN_AsiaNE.EuroN,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN <- windowStats$Fst_Asia_SW_Euro_N - windowStats$Fst_Asia_Sichuan_Euro_N
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN = ifelse(Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_EuroN - AsiaSichuan_EuroN",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN))
abline(h=mean(Delta_df$D.AsiaSW.EuroN_AsiaSichuan.EuroN,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

#Asia Noram - Asia CircNoram
Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram <- windowStats$Fst_Asia_Sichuan_North_American - windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_E
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSC.Noram_AsiaSC.CircNoram = ifelse(Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSC_Noram - AsiaSC_CircNoram",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSC.Noram_AsiaSC.CircNoram))
abline(h=mean(Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df$D.AsiaSW.Noram_AsiaSW.CircNoram <- windowStats$Fst_Asia_SW_North_American - windowStats$Fst_Asia_SW_Circumboreal_NorAm_E
# identify the 95% and 99% percentile
quantile(Delta_df$D.AsiaSW.Noram_AsiaSW.CircNoram, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(Delta_df$D.AsiaSW.Noram_AsiaSW.CircNoram, 0.975, na.rm = T)
Delta_df <- Delta_df %>% mutate(outlier.D.AsiaSW.Noram_AsiaSW.CircNoram = ifelse(Delta_df$D.AsiaSC.Noram_AsiaSC.CircNoram > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.D.AsiaSW.Noram_AsiaSW.CircNoram) %>% tally()
plot(Delta_df$Mid,Delta_df$D.AsiaSW.Noram_AsiaSW.CircNoram,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Delta_Fst AsiaSW_Noram - AsiaSW_CircNoram",ylim=c(-1,1), col=as.factor(Delta_df$outlier.D.AsiaSW.Noram_AsiaSW.CircNoram))
abline(h=mean(Delta_df$D.AsiaSW.Noram_AsiaSW.CircNoram,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

Delta_df_between <- Delta_df

#Fst outlier Noram:
# identify the 95% and 99% percentile
quantile(windowStats$Fst_Circumboreal_NorAm_E_North_American, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(windowStats$Fst_Circumboreal_NorAm_E_North_American, 0.975, na.rm = T)

Delta_df <- Delta_df %>% mutate(outlier.fst_Noram.CircnoramE = ifelse(windowStats$Fst_Circumboreal_NorAm_E_North_American > my_threshold, "outlier", "background"))
Delta_df %>% group_by(outlier.fst_Noram.CircnoramE) %>% tally()
which(Delta_df$outlier.fst_Noram.CircnoramE=="outlier")

#############################
#                           #
# Explore Delta_df          #
#                           #
#############################


#all outliers:
str(Delta_df)
Delta_df_outlier_columns <- grep("outlier", Delta_df)
Delta_df_outlier_columns
Delta_df_outliers <- Delta_df[, c(1,2,3,4,5,6,Delta_df_outlier_columns)]
colnames(Delta_df_outliers)
dim(Delta_df_outliers)

Delta_df_all_outliers <- Delta_df %>% filter(outlier.D.Asia_SW.Sichuan_SW.NE =="outlier" | 
                           outlier.D.Asia_SW.Sichuan_Sichuan.NE =="outlier" | 
                           outlier.D.Asia_SW.N_N.NE =="outlier" |
                           outlier.D.Asia_SW.N_N.Sichuan =="outlier" |
                           outlier.D.Noram_NoramE.CircNoramE_CircNoramE.Noram2 =="outlier" |
                           outlier.D.AsiaSW.NoramE_AsiaNE.NoramE =="outlier" |
                           outlier.D.AsiaSW.NoramE_AsiaN.NoramE =="outlier" |
                           outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE =="outlier" |
                           outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE =="outlier" |
                           outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE =="outlier" |
                           outlier.D.AsiaSW.CircNoramW_AsiaNE.CircNoramW =="outlier" |
                           outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW =="outlier" |
                           outlier.D.AsiaSW.EuroN_AsiaNE.EuroN =="outlier" |
                           outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN =="outlier" |
                           outlier.D.AsiaSW.EuroN_AsiaNE.EuroN =="outlier" |
                           outlier.D.AsiaSC.Noram_AsiaSC.CircNoram =="outlier" |
                           outlier.D.AsiaSW.Noram_AsiaSW.CircNoram =="outlier")
 
#write output to file:                          
write.table(Delta_df_all_outliers, file="Delta_df_all_outliers.txt", sep = "\t", quote = FALSE, row.names = FALSE)
Delta_df_all_outliers_pos <- Delta_df_all_outliers[,1:3]
write.table(Delta_df_all_outliers_pos, file="Delta_df_all_outliers_pos.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#overall comparison of overlap between asia within, noram within, and between continents:

Delta_df_asia_within_outliers <- which(Delta_df_asia_within %>% filter(outlier.D.Asia_SW.Sichuan_Sichuan.NE=="outlier" | 
                                                                   outlier.D.Asia_SW.Sichuan_SW.NE=="outlier" |
                                                                   outlier.D.Asia_SW.N_N.NE=="outlier" |
                                                                   outlier.D.Asia_SW.N_N.Sichuan=="outlier")


my_groups <- list("asia"=c(which()), 
                  "noram"=c(which()),
                  "between"=c(which()))

Asia <- list("SW_SC - SW_NE"=c(which(Delta_df$outlier.D.Asia_SW.Sichuan_SW.NE=="outlier")),
             "SW_SC - SC_NE"=c(which(Delta_df$outlier.D.Asia_SW.Sichuan_Sichuan.NE=="outlier")),
             "SW_N - N_NE"=c(which(Delta_df$outlier.D.Asia_SW.N_N.NE=="outlier")),
             "SW_N - N_SC"=c(which(Delta_df$outlier.D.Asia_SW.N_N.Sichuan=="outlier")))
length(Asia)
ggvenn(Asia, show_percentage=FALSE,set_name_size = 4)
ggvenn(Asia, c("D.Asia_SW.Sichuan_SW.NE","D.Asia_SW.Sichuan_Sichuan.NE"))
ggVennDiagram(Asia, label_alpha = 0)



Delta_df_asia_within

Delta_df_Noram_within

Delta_df_between

#outliers Asia Euro N:
which(Delta_df$outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN=="outlier" & Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier")

#outliers Asia CircNoram:
which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier" & Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaNE.CircNoramW=="outlier")
which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE=="outlier" & Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW=="outlier")

which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier"
      & Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW=="outlier")

Delta_df[3165,]
#Scaffold05 35001 45000

Delta_df[4235,]
#Scaffold05 6140001 6150000

Delta_df[6260,]
#Scaffold09 2335001

#outliers Asia Noram
which(Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE=="outlier" 
      & Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE=="outlier")

which(Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE=="outlier" & 
        Delta_df$outlier.D.AsiaSW.NoramE_AsiaN.NoramE=="outlier" &
        Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE=="outlier" & 
        Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE=="outlier")
  #791 5907 5965 6572
Delta_df[791,]
#Scaffold02 310001 320000
Delta_df[5907,]
#Scaffold09 135001 145000
Delta_df[5965,]
#Scaffold09 500001 510000
Delta_df[6572,]
#Scaffold10 375001 385000

Asia.Noram <- list("AsiaSW_Noram - AsiaNE_Noram"=c(which(Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE=="outlier")),
                   "AsiaSW_Noram - AsiaN_Noram"=c(which(Delta_df$outlier.D.AsiaSW.NoramE_AsiaN.NoramE=="outlier")),
                  "AsiaSC_Noram - AsiaNE_Noram"=c(which(Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE=="outlier")),
                  "AsiaSC_Noram - AsiaN_Noram"=c(which(Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaN.NoramE=="outlier")))
ggvenn(Asia.Noram, show_percentage=FALSE,set_name_size = 4)

#outliers AsiaNoram - AsiaCirNoram
which(Delta_df$outlier.D.AsiaSW.Noram_AsiaSW.CircNoram=="outlier" &
         Delta_df$outlier.D.AsiaSC.Noram_AsiaSC.CircNoram=="outlier")
        
#outliers within Asia:
which(Delta_df$outlier.D.Asia_SW.Sichuan_SW.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.Sichuan_Sichuan.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.N_N.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.N_N.Sichuan=="outlier") 
Delta_df[4291,]
#Scaffold06 130001 140000
#protein of unknown function
Delta_df[4292,]
#Scaffold06 135001 145000

#is this also an outlier between noram and circnoram?
#no
Asia <- list("SW_SC - SW_NE"=c(which(Delta_df$outlier.D.Asia_SW.Sichuan_SW.NE=="outlier")),
             "SW_SC - SC_NE"=c(which(Delta_df$outlier.D.Asia_SW.Sichuan_Sichuan.NE=="outlier")),
             "SW_N - N_NE"=c(which(Delta_df$outlier.D.Asia_SW.N_N.NE=="outlier")),
             "SW_N - N_SC"=c(which(Delta_df$outlier.D.Asia_SW.N_N.Sichuan=="outlier")))
length(Asia)
ggvenn(Asia, show_percentage=FALSE,set_name_size = 4)
ggvenn(Asia, c("D.Asia_SW.Sichuan_SW.NE","D.Asia_SW.Sichuan_Sichuan.NE"))
ggVennDiagram(Asia, label_alpha = 0)

#outliers Asia EuroN & AsiaCircNoram:
which(Delta_df$outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN=="outlier" 
      & Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier" 
      & Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier"
      & Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW=="outlier")

Delta_df[6260,]
#Scaffold09 2335001 2345000
#4 genes:
  #1.Similar to MRD1: Multiple RNA-binding domain-containing protein 1 (Ustilago maydis (strain 521 / FGSC 9021)OX=237631)
  #2. unknown function
  #3. unknown function
  #4. Similar to SPT5: Transcription elongation factor SPT5 (Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923/ NRRL Y-1056) OX=284811
#leucine-rich repeat receptor-like serine/threonine-protein kinase 1272099	1276238

which(Delta_df$outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN=="outlier" 
      & Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier" 
      & Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE=="outlier")
#6260 6554 6555
Delta_df[6554,] #Scaffold10 175001 185000
Delta_df[6555,] #Scaffold10 180001 190000
#rcd1: 896396	898399


Asia.CircNoram.Euro <- list("AsiaSW_EuroN - AsiaSC_EuroN"=c(which(Delta_df$outlier.D.AsiaSW.EuroN_AsiaSichuan.EuroN=="outlier")),
                            "AsiaSW_EuroN - AsiaNE_EuroN"=c(which(Delta_df$outlier.D.AsiaSW.EuroN_AsiaNE.EuroN=="outlier")),
                            "AsiaSW_CircNoramE - AsiaNE_CircNoramE"=c(which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier")),
                            "AsiaSW_CirNoramE - AsiaSC_CircNoramE"=c(which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE=="outlier")))
ggvenn(Asia.CircNoram.Euro, show_percentage=FALSE,set_name_size = 3)

#outliers Asia Within  & Asia Noram:
which(Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE=="outlier" 
      & Delta_df$outlier.D.AsiaSichuan.NoramE_AsiaNE.NoramE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.Sichuan_SW.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.Sichuan_Sichuan.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.N_N.NE=="outlier" 
      & Delta_df$outlier.D.Asia_SW.N_N.Sichuan=="outlier") 
Delta_df[4292,]
#Scaffold06 135001 145000


#outliers across all groups
which(Delta_df$outlier.D.AsiaSW.CircNoramE_AsiaNE.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSichuan.CircNoramE_AsiaNE.CircNoramE=="outlier"
      & Delta_df$outlier.D.AsiaSW.NoramE_AsiaNE.NoramE=="outlier"
      & Delta_df$outlier.D.AsiaSC.NoramE_AsiaNE.NoramE=="outlier")
      


str(windowStats)
Delta_df <- as.data.frame(cbind(windowStats$scaffold,windowStats$start,windowStats$end, windowStats$mid, windowStats$sites,
              D.AsiaSW.CircNoramE_AsiaNE.CircNoramE,
              D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE,
              D.AsiaSW.CircNoramW_AsiaNE.CircNoramW,
              D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW,
              D.AsiaSW.EuroN_AsiaNE.EuroN,
              D.AsiaSW.EuroN_AsiaSichuan.EuroN
              ))

str(Delta_df)
colnames(Delta_df)[1] <-"Scaffold"
colnames(Delta_df)[2] <-"Start"
colnames(Delta_df)[3] <-"End"
colnames(Delta_df)[4] <-"Mid"
colnames(Delta_df)[5] <- "Sites"
str(Delta_df)

# identify the 95% and 99% percentile
quantile(as.numeric(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE), c(0.975, 0.995), na.rm = T)

# identify the 95% percentile
my_threshold <- quantile(as.numeric(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE), 0.975, na.rm = T)
# make an outlier column in the data.frame
Delta_df <- Delta_df %>% mutate(outlier = ifelse(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE > my_threshold, "outlier", "background"))
str(Delta_df)

Delta_df %>% group_by(outlier) %>% tally()

par(mfrow=c(1,1))
plot(Delta_df$Mid, Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="all",ylim=c(-1,1), col="red")
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
points(Delta_df$Mid, Delta_df$D.AsiaSW.CircNoramE_AsiaSichuan.CircNoramE, col="blue")
points(Delta_df$Mid,D.AsiaSW.CircNoramW_AsiaNE.CircNoramW, col="green")
points(Delta_df$Mid,D.AsiaSW.CircNoramW_AsiaSichuan.CircNoramW, col="orange")


which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578
which.max(Delta_df$D.AsiaSW.CircNoramE_AsiaNE.CircNoramE) #578

which.max(scf9$Fst_Asia_N_Asia_SW) #578
[587,] #start: 3625001 end: 3635000 mid: 39880262

## make venn diagrams:
#install.packages("ggvenn")
library("ggvenn")

#########################
#                       #
# Plots regular Fst     #
#                       #
#########################

palette(c("grey40","grey60"))

  #Asia SW Asia Sichuan
par(mfrow=c(2,2))

palette(c("#1F78B4","#A6CEE3"))
plot(windowStats$mid,windowStats$Fst_Asia_SW_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_EastAsia_SW_EastAsia_SC"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_Asia_Sichuan,na.rm=T)
#0.1406599

#Asia SW Asia N
plot(windowStats$mid,windowStats$Fst_Asia_N_Asia_SW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_N_Asia_SW"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Asia_SW,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_Asia_SW,na.rm=T)
#0.1564485

##North America
  #North America Circumboreal E North America
pdf("noramA_noramBE.pdf", height=4, width=12)
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_E_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Circumboreal_NorAm_E_North_American"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_E_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Circumboreal_NorAm_E_North_American,na.rm=T)
dev.off()
#0.4167556

#plot only scaffold5:
scaffold5 <- filter(windowStats, windowStats$scaffold == "Scaffold05")
plot(scaffold5$start,scaffold5$Fst_Circumboreal_NorAm_E_North_American,cex=0.5,pch=19,xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Circumboreal_NorAm_E_North_American"
     ,ylim=c(0,1), col="darkgrey")
abline(v=3213128, col="red") #half point
abline(v=3110000,col="blue") 
abline(v=3620000, col= "green")     
       
       
#North America Circumboreal W North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_W_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Circumboreal_NorAm_E_North_American"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_W_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Circumboreal_NorAm_W_North_American,na.rm=T)
#0.4621608

  #North America North America 2#
plot(windowStats$mid,windowStats$Fst_North_American_North_American_2,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_North_American_North_American_2",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_North_American_North_American_2,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_North_American_North_American_2,na.rm=T)
#0.08572177

##Asia North America
  #Asia SW North America
plot(windowStats$mid,windowStats$Fst_Asia_SW_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_SW_North_American"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_North_American,na.rm=T)
#0.4149724

  #Asia Sichuan North America
plot(windowStats$mid,windowStats$Fst_Asia_Sichuan_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_Sichuan_North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_Sichuan_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_Sichuan_North_American,na.rm=T)
#0.2910003

#Fst_Asia_SW_Circumboreal_NorAm_E
plot(windowStats$mid,windowStats$Fst_Asia_SW_Circumboreal_NorAm_E,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_SW_Circumboreal_NorAm_E"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_Circumboreal_NorAm_E,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_Circumboreal_NorAm_E,na.rm=T)
#0.2063156

#Fst_Asia_SW_North_American_2
plot(windowStats$mid,windowStats$Fst_Asia_SW_North_American_2,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_SW_North_American_2",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_North_American_2,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_North_American_2,na.rm=T)
#0.2155965

###############################################
#                                             #
#Comparisons between groups which do mate     #
#                                             #  
###############################################

  ##Asia
#Fst_Asia_NE_Asia_SW
plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_SW,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_EastAsia_NE_EastAsia_SW"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_SW,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Asia_SW,na.rm=T)
#0.102311

#Fst_Asia_NE_Asia_Sichuan
plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Asia_Sichuan"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Asia_Sichuan,na.rm=T)
#0.1807297

#Fst_Asia_N_Asia_Sichuan
plot(windowStats$mid,windowStats$Fst_Asia_N_Asia_Sichuan,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_N_Asia_Sichuan"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Asia_Sichuan,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_Asia_Sichuan,na.rm=T)
#0.1132778

#Fst_Asia_NE_Asia_N
palette("paired")
palette(c("#E31A1C","#FB9A99"))
plot(windowStats$mid,windowStats$Fst_Asia_NE_Asia_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Asia_N"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Asia_N,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Asia_N,na.rm=T)
#0.1638262

##North America
  #North America early diverging X Circumboreal North America E
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_E_North_American_2,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Circumboreal_NorAm_E_North_American_2",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_E_North_American_2,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Circumboreal_NorAm_E_North_American_2,na.rm=T)
#0.06302359

#Circumboreal North America W X Circumboreal North America E
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W,na.rm=T)
#0.1341001

##Asia North America
  #Asia_N x North_American
plot(windowStats$mid,windowStats$Fst_Asia_N_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_N_North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_North_American,na.rm=T)
#0.3132754

  #Asia_NE x North_American
plot(windowStats$mid,windowStats$Fst_Asia_NE_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_North_American"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_North_American,na.rm=T)
#0.4867606

  #Asia_NE x Circumboreal_NorAm_E
plot(windowStats$mid,windowStats$Fst_Asia_NE_Circumboreal_NorAm_E,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Circumboreal_NorAm_E"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Circumboreal_NorAm_E,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Circumboreal_NorAm_E,na.rm=T)
#0.2775764

  #Asia_N x Circumboreal_NorAm_E
plot(windowStats$mid,windowStats$Fst_Asia_N_Circumboreal_NorAm_E,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_N_Circumboreal_NorAm_E",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Circumboreal_NorAm_E,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_Circumboreal_NorAm_E,na.rm=T)
#0.0834743

  #Asia_Sichuan x Circumboreal_NorAm_E
plot(windowStats$mid,windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_E,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_Sichuan_Circumboreal_NorAm_E",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_E,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_E,na.rm=T)
#0.100946

###########################################
#                                         #
#             Comparisons with Europe     #                          
#                                         #
###########################################

palette("paired")
palette(c("#E31A1C","#FB9A99"))
## Noram x Europe N
#NorAm E x Euro N
plot(windowStats$mid,windowStats$Fst_Euro_N_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Euro_N_Circumboreal_NorAm_E"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Euro_N_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Euro_N_North_American,na.rm=T)
#0.4061016

  #Circumboreal NorAm E x Euro N
plot(windowStats$mid,windowStats$Fst_Euro_N_Circumboreal_NorAm_E,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Euro_N_Circumboreal_NorAm_E"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Euro_N_Circumboreal_NorAm_E,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Euro_N_Circumboreal_NorAm_E,na.rm=T)
#0.2551512

#NorAm E x Euro N
pdf("noramA_eurasia.pdf", height=4, width=12)
plot(windowStats$mid,windowStats$Fst_Euro_N_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Euro_N_Circumboreal_NorAm_E"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Euro_N_North_American,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Euro_N_North_American,na.rm=T)
dev.off()
#0.4061016

### Asia x Europe N
palette("paired")
palette()
palette(c("#1F78B4","#A6CEE3"))
  #Asia SW x Europe N 
plot(windowStats$mid,windowStats$Fst_Asia_SW_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_SW_Euro_N"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_Euro_N,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_Euro_N,na.rm=T)
#0.1801456

  #Asia NE x Europe N
plot(windowStats$mid,windowStats$Fst_Asia_NE_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     #main="Fst_Asia_NE_Euro_N"
     ,ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Euro_N,na.rm=T),col="grey2", lwd=2, lty=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Euro_N,na.rm=T)
#0.2477928
  #Asia N x Europe N
plot(windowStats$mid,windowStats$Fst_Asia_N_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_N_Euro_N",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Euro_N,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_Euro_N,na.rm=T)
#0.07182993
  #Asia Sichuan x Europe N
plot(windowStats$mid,windowStats$Fst_Asia_Sichuan_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_Sichuan_Euro_N",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_Sichuan_Euro_N,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_Sichuan_Euro_N,na.rm=T)
#0.08445428

### Asia x Europe S
#Asia SW x Europe S
plot(windowStats$mid,windowStats$Fst_Asia_SW_Euro_S,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_SW_Euro_S",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_SW_Euro_S,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_SW_Euro_S,na.rm=T)
#0.4114358
#Asia Sichuan x Europe S
plot(windowStats$mid,windowStats$Fst_Asia_Sichuan_Euro_S,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_Sichuan_Euro_S",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_Sichuan_Euro_S,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_Sichuan_Euro_S,na.rm=T)
#0.3362468
#Asia NE x Europe S
plot(windowStats$mid,windowStats$Fst_Asia_NE_Euro_S,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_NE_Euro_S",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_NE_Euro_S,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_NE_Euro_S,na.rm=T)
#0.4521234
#Asia N x Europe S
plot(windowStats$mid,windowStats$Fst_Asia_N_Euro_S,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst_Asia_N_Euro_S",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Asia_N_Euro_S,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)
mean(windowStats$Fst_Asia_N_Euro_S,na.rm=T)
#0.3491728

#Get rest of mean Fst values:
mean(windowStats$Fst_Asia_N_Euro_Mid,na.rm=T)
#0.1547135
mean(windowStats$Fst_Euro_N_North_American_S,na.rm=T)
#0.1728753
mean(windowStats$Fst_Asia_NE_Euro_Mid,na.rm=T)
#0.3296264
mean(windowStats$Fst_Asia_SW_Euro_Mid,na.rm=T)
#0.2668768
mean(windowStats$Fst_Asia_Sichuan_Euro_Mid,na.rm=T)
#0.1556292
mean(windowStats$Fst_Euro_Mid_Euro_N,na.rm=T)
#0.06488082
mean(windowStats$Fst_Euro_Mid_North_American,na.rm=T)
#0.3956119
mean(windowStats$Fst_Euro_Mid_North_American_S,na.rm=T)
#0.2668447
mean(windowStats$Fst_Euro_Mid_North_American_2,na.rm=T)
#0.07890068
mean(windowStats$Fst_Euro_Mid_Circumboreal_NorAm_E,na.rm=T)
#0.244127
mean(windowStats$Fst_Euro_Mid_Circumboreal_NorAm_W,na.rm=T)
#0.2857113
mean(windowStats$Fst_Euro_S_Circumboreal_NorAm_W,na.rm=T)
#0.4511593
mean(windowStats$Fst_Euro_S_Circumboreal_NorAm_E,na.rm=T)
#0.2883318
mean(windowStats$Fst_Euro_S_Euro_N,na.rm=T)
#0.1736305
mean(windowStats$Fst_Euro_S_Euro_Mid,na.rm=T)
#0.1034673
mean(windowStats$Fst_Euro_S_North_American,na.rm=T)
#0.4083817
mean(windowStats$Fst_Euro_S_North_American_2,na.rm=T)
#0.1990372
mean(windowStats$Fst_Euro_S_North_American_S,na.rm=T)
#0.4311234
mean(windowStats$Fst_Euro_N_North_American_S,na.rm=T)
#0.1728753
mean(windowStats$Fst_Euro_N_North_American_2,na.rm=T)
#0.04497038
mean(windowStats$Fst_Euro_N_Circumboreal_NorAm_W,na.rm=T)
#0.1873957
mean(windowStats$Fst_North_American_North_American_S,na.rm=T)
#0.0730569
mean(windowStats$Fst_Circumboreal_NorAm_W_North_American,na.rm=T)
#0.4621608
mean(windowStats$Fst_Circumboreal_NorAm_W_North_American_S,na.rm=T)
#0.4434679
mean(windowStats$Fst_Circumboreal_NorAm_W_North_American_2,na.rm=T)
#0.1953608
mean(windowStats$Fst_Circumboreal_NorAm_E_North_American_S,na.rm=T)
#0.2135186
mean(windowStats$Fst_North_American_2_North_American_S,na.rm=T)
#0.2269915
mean(windowStats$Fst_Asia_SW_North_American_S,na.rm=T)
#0.4378888
mean(windowStats$Fst_Asia_Sichuan_North_American_S,na.rm=T)
#0.4217726
mean(windowStats$Fst_Asia_NE_North_American_S,na.rm=T)
#0.4546826
mean(windowStats$Fst_Asia_N_North_American_S,na.rm=T)
#0.4252672

mean(windowStats$Fst_Asia_SW_Circumboreal_NorAm_W,na.rm=T)
#0.3508238
mean(windowStats$Fst_Asia_Sichuan_Circumboreal_NorAm_W,na.rm=T)
#0.2402303
mean(windowStats$Fst_Asia_NE_Circumboreal_NorAm_W,na.rm=T)
#0.406608
mean(windowStats$Fst_Asia_N_Circumboreal_NorAm_W,na.rm=T)
#0.2003453

mean(windowStats$Fst_Asia_SW_North_American_2,na.rm=T)
#0.2155965
mean(windowStats$Fst_Asia_Sichuan_North_American_2,na.rm=T)
#0.2220209
mean(windowStats$Fst_Asia_NE_North_American_2,na.rm=T)
#0.2227738
mean(windowStats$Fst_Asia_N_North_American_2,na.rm=T)
#0.2134425

#dev.off()
###################################
#                                 #
#   Investigate specific windows  #
#                                 #
###################################

## 10 kb
scf9 <- subset(windowStats, scaffold=="Scaffold09")
which.max(scf9$Fst_Asia_SW_Asia_Sichuan) #578
which.max(scf9$Fst_Asia_N_Asia_SW) #578
scf9[587,] #start: 3625001 end: 3635000 mid: 39880262

scf6 <- subset(windowStats, scaffold=="Scaffold06")
which.max(scf6$Fst_Asia_SW_Asia_Sichuan) #17
which.max(scf6$Fst_Asia_SW_North_American_2) #152
which.max(scf6$Fst_Asia_SW_North_American) #134
which.max(scf6$Fst_Asia_NE_Circumboreal_NorAm_E) #312
which.max(scf6$Fst_Asia_N_Euro_N) #176
which.max(scf6$Fst_Asia_SW_Euro_S) #176
which.max(scf6$Fst_Euro_S_Circumboreal_NorAm_E) #123


### 20kb

### Scaffold9 ###
scf9 <- subset(windowStats, scaffold=="Scaffold09")
str(scf9)


scf9[145,] #start: 3540001 #end:3560000 #sites: 258

#Asia
which.max(scf9$Fst_Asia_SW_Asia_Sichuan) #145
which.max(scf9$Fst_Asia_N_Asia_SW) #145
which.max(scf9$Fst_Asia_NE_Asia_Sichuan) #15
which.max(scf9$Fst_Asia_NE_Asia_SW) #145
which.max(scf9$Fst_Asia_NE_Asia_N) #106

#North America
which.max(scf9$Fst_North_American_North_American_2) #89
new <- scf9[-89,]
which.max(new$Fst_North_American_North_American_2) #149
which.max(scf9$Fst_Circumboreal_NorAm_E_North_American) #76
which.max(scf9$Fst_Circumboreal_NorAm_E_North_American_2) #89
which.max(scf9$Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W) #145

#North America x Asia
which.max(scf9$Fst_Asia_SW_North_American) #89
which.max(scf9$Fst_Asia_Sichuan_North_American) #89
which.max(scf9$Fst_Asia_N_North_American) #89
which.max(scf9$Fst_Asia_NE_North_American) #89
which.max(scf9$Fst_Asia_SW_Circumboreal_NorAm_E) #128
new <- scf9[-128,]
which.max(new$Fst_Asia_SW_Circumboreal_NorAm_E) #144
which.max(scf9$Fst_Asia_Sichuan_Circumboreal_NorAm_E) #145
which.max(scf9$Fst_Asia_NE_Circumboreal_NorAm_E) #128
which.max(scf9$Fst_Asia_N_Circumboreal_NorAm_E) #128
which.max(scf9$Fst_Asia_SW_North_American_2) #89

#Europe N x Asia
which.max(scf9$Fst_Asia_SW_Euro_N) #106
which.max(scf9$Fst_Asia_NE_Euro_N) #106
which.max(scf9$Fst_Asia_N_Euro_N) #107
which.max(scf9$Fst_Asia_Sichuan_Euro_N) #95

#Europe S  x Asia
which.max(scf9$Fst_Asia_SW_Euro_S) #113
which.max(scf9$Fst_Asia_NE_Euro_S) #77
which.max(scf9$Fst_Asia_N_Euro_S) #77
which.max(scf9$Fst_Asia_Sichuan_Euro_S) #77

### Scaffold06  ###
scf6 <- subset(windowStats, scaffold=="Scaffold06")
str(scf6)

scf6[38,] #start:1040001 #end:1060000 #sites: 166

#Asia
which.max(scf6$Fst_Asia_SW_Asia_Sichuan) #38
which.max(scf6$Fst_Asia_N_Asia_SW) #77
which.max(scf6$Fst_Asia_NE_Asia_Sichuan) #6
which.max(scf6$Fst_Asia_NE_Asia_SW) #85
which.max(scf6$Fst_Asia_NE_Asia_N) #77

#North America
which.max(scf6$Fst_Circumboreal_NorAm_E_Circumboreal_NorAm_W) #20
which.max(scf6$Fst_Circumboreal_NorAm_E_North_American) #13
new <- scf6[-13,]
which.max(new$Fst_Circumboreal_NorAm_E_North_American) #77

#North America x Asia 
which.max(scf6$Fst_Asia_SW_North_American) #38
which.max(scf6$Fst_Asia_NE_North_American) #78
new <- scf6[-78,]
which.max(new$Fst_Asia_NE_North_American) #38
which.max(scf6$Fst_Asia_N_North_American) #38
which.max(scf6$Fst_Asia_Sichuan_North_American) #78
new <- scf6[-78,]
which.max(new$Fst_Asia_Sichuan_North_American) #38

#Circumboreal North America x Asia 
which.max(scf6$Fst_Asia_SW_Circumboreal_NorAm_E) #77
which.max(scf6$Fst_Asia_SW_North_American_2) #38
which.max(scf6$Fst_Asia_Sichuan_Circumboreal_NorAm_E) #77
which.max(scf6$Fst_Asia_NE_Circumboreal_NorAm_E) #77

#Europe N x Asia
which.max(scf6$Fst_Asia_SW_Euro_N) #44
new <- scf6[-44,]
which.max(new$Fst_Asia_SW_Euro_N) #76
which.max(scf6$Fst_Asia_NE_Euro_N) #77
which.max(scf6$Fst_Asia_N_Euro_N) #38
which.max(scf6$Fst_Asia_Sichuan_Euro_N) #77

#Europe S  x Asia
which.max(scf6$Fst_Asia_SW_Euro_S) #38
which.max(scf6$Fst_Asia_NE_Euro_S) #44
which.max(scf6$Fst_Asia_N_Euro_S) #45
which.max(scf6$Fst_Asia_Sichuan_Euro_S) #37



### North America within  ###
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_E_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_NorAm_E & North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_E_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_W_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_NorAm_W & North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_W_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

### Asia North America  ###
plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm_W,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_NorAm_W & Asia",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm_W,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst North America & Asia",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

### Asia Europe ###
plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Asia & Circumboreal_Euro_N",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro_N,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,na.rm=T)
mean(windowStats$Fst_Circumboreal_NorAm_North_American,na.rm=T)

#############################################
#                                           #  
#     Main Pop designatoin                  #
#                                           #  
#############################################

#################
#               #
#       pi      #
#               #
#################

names(windowStats)
par(mfrow=c(4,1),oma=c(3,0,0,0),mar=c(4,5,2,2))

#North American
plot(windowStats$mid,windowStats$pi_North_American,cex=0.5,pch=19,xaxt="n",ylab="",xlab="",cex.main=2, cex.lab=2,
     main="pi_North_American",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$pi_North_American,na.rm=T),col="grey", lwd=2)
# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$pi_North_American,na.rm=T)

#North American 2
#plot(windowStats$mid,windowStats$pi_North_American_2,cex=0.5,pch=19,xaxt="n",
     #main="pi_North_American_2",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
#abline(h=mean(windowStats$pi_North_American_2,na.rm=T),col="grey", lwd=2)
# Add the LG names to the center of each LG
#axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

#mean(windowStats$pi_North_American_2,na.rm=T)

#Circumboreal_NorAm
plot(windowStats$mid,windowStats$pi_Circumboreal_NorAm,cex=0.5,pch=19,xaxt="n",ylab="",xlab="",cex.main=2, cex.lab=2,
     main="pi_Circumboreal_NorAm",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$pi_Circumboreal_NorAm,na.rm=T),col="grey", lwd=2)
# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$pi_Circumboreal_NorAm,na.rm=T)

#pi_Circumboreal_Asia
plot(windowStats$mid,windowStats$pi_Circumboreal_Asia,cex=0.5,pch=19,xaxt="n",ylab="",xlab="",cex.main=2, cex.lab=2,
     main="pi_Circumboreal_Asia",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$pi_Circumboreal_Asia,na.rm=T),col="grey", lwd=2)
# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$pi_Circumboreal_Asia,na.rm=T)

#pi_Circumboreal_Euro
plot(windowStats$mid,windowStats$pi_Circumboreal_Euro,cex=0.5,pch=19,xaxt="n",ylab="",xlab="",cex.main=2, cex.lab=2,
     main="pi_Circumboreal_Euro",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$pi_Circumboreal_Euro,na.rm=T),col="grey", lwd=2)
# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$pi_Circumboreal_Euro,na.rm=T)

#################
#               #
#       Fst     #
#               #
#################

names(windowStats)
ggplot(windowStats, aes(mid, Fst_Circumboreal_NorAm_North_American, color=scaffold)) + geom_line()

#pdf("Fst_windows.pdf", height=12, width=12)
par(mfrow=c(5,1),oma=c(3,0,0,0),mar=c(4,5,2,2))

# Plot Fst between Populations in North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_NorAm_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_NorAm & North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_NorAm_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

# Plot Fst between Asia and North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Asia & North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Asia_North_American,na.rm=T)

# Plot Fst between Europe and North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_Euro_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Euro & North_American",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Euro_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Euro_North_American,na.rm=T)

# Plot Fst between Asia and Circumboreal North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Asia & Circumboreal_NorAm",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_NorAm,na.rm=T)

# Plot Fst between Europe and Circumboreal North America
plot(windowStats$mid,windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Europe & Circumboreal_NorAm",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm,na.rm=T)

# Plot Fst between Asia and Europe
plot(windowStats$mid,windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Fst Circumboreal_Asia & Circumboreal_Euro",ylim=c(0,1), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,na.rm=T)

# Plot Fst between Populations in Europe
#plot(windowStats$mid,windowStats$Fst_Circumboreal_Euro_S_Circumboreal_Euro_N,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",
     #main="Fst Circumboreal_Euro_S & Circumboreal_Euro_N",ylim=c(0,1), col=as.factor(windowStats$scaffold))
#abline(h=mean(windowStats$Fst_Circumboreal_Euro_S_Circumboreal_Euro_N,na.rm=T),col="grey", lwd=2)
#axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)


#################
#               #
#       Dxy     #
#               #
#################

par(mfrow=c(6,1),oma=c(3,0,0,0),mar=c(4,5,2,2))

# Plot Dxy between Populations in North America
plot(windowStats$mid,windowStats$dxy_Circumboreal_NorAm_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_NorAm & North_American",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_NorAm_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_NorAm_North_American,na.rm=T)

# Plot Dxy between Asia and North America
plot(windowStats$mid,windowStats$dxy_Circumboreal_Asia_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_Asia & North_American",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_Asia_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_Asia_North_American,na.rm=T)

# Plot Fst between Europe and North America
plot(windowStats$mid,windowStats$dxy_Circumboreal_Euro_North_American,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_Euro & North_American",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_Euro_North_American,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_Euro_North_American,na.rm=T)

# Plot Fst between Asia and Circumboreal North America
plot(windowStats$mid,windowStats$dxy_Circumboreal_Asia_Circumboreal_NorAm,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_Asia & Circumboreal_NorAm",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_NorAm,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_NorAm,na.rm=T)

# Plot Fst between Europe and Circumboreal North America
plot(windowStats$mid,windowStats$dxy_Circumboreal_Euro_Circumboreal_NorAm,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_Europe & Circumboreal_NorAm",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_Euro_Circumboreal_NorAm,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_Euro_Circumboreal_NorAm,na.rm=T)

# Plot Fst between Asia and Europe
plot(windowStats$mid,windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,cex=0.5,pch=19,xaxt="n",xlab="",ylab="",cex.main=2, cex.lab=2,
     main="Dxy Circumboreal_Asia & Circumboreal_Euro",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
abline(h=mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,na.rm=T),col="grey", lwd=2)
axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

mean(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,na.rm=T)

# Plot Dxy between Populations in Europe
#plot(windowStats$mid,windowStats$dxy_Circumboreal_Euro_S_Circumboreal_Euro_N,cex=0.5,pch=19,xaxt="n",
     #ylab="dxy_Circumboreal_Euro_S_Circumboreal_Euro_N",ylim=c(0,0.5), col=as.factor(windowStats$scaffold))
#abline(h=mean(windowStats$dxy_Circumboreal_Euro_S_Circumboreal_Euro_N,na.rm=T),col="grey",lwd=2)
# Add the LG names to the center of each LG
#axis(1,at=chrom$add+chrom$end/2,tick = T,labels = 1:12)

# Are dxy and fst correlated?
par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(1,1,1,1))

plot(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,
     xlab="dxy",ylab="Fst",pch=19,cex=0.3)
# ignoring the outliers:
plot(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,
     xlab="dxy",ylab="Fst",pch=19,cex=0.3,xlim=c(0,0.01))

# Is dxy correlated with pi, indicating that it is mostly driven by variance in mutation and recombination rates?
plot(rowMeans(cbind(windowStats$pi_Circumboreal_Asia,windowStats$pi_Circumboreal_Euro),na.rm=T),
     windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,cex=0.3,
     xlab="pi",ylab="dxy")
abline(a=0,b=1,col="grey")

# If we correct for mutation rate differences with dxy to kivu
plot(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro/rowMeans(cbind(windowStats$dxy_Circumboreal_Asia_Circumboreal_Euro,
                                                    windowStats$dxy_Circumboreal_Euro_Circumboreal_NorAm),na.rm=T),
     windowStats$Fst_Circumboreal_Asia_Circumboreal_Eur,ylab="Fst",xlab="dxy normalized",pch=19,cex=0.3)


##############################################################################

# Supplementary


# What are the FST distributions between the two species pairs?
par(mfrow=c(2,2))
boxplot(windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro,windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm)
abline(h=0)

# Are shared high Fst windows in regions of low recombination?
sharedHighFst<-windowStats[windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro > 0.1 &
                             windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm > 0.1, ]
dim(sharedHighFst)

normalFST<-windowStats[windowStats$Fst_Circumboreal_Asia_Circumboreal_Euro<0.1 &
                         windowStats$Fst_Circumboreal_Euro_Circumboreal_NorAm<0.1,]

dim(normalFST)

# Add the ratio of dxy between Makobe to dxy to kivu
sharedHighFst$ratio <-sharedHighFst$dxy_Circumboreal_Asia_Circumboreal_Euro/sharedHighFst$dxy_Circumboreal_Euro_Circumboreal_NorAm

sharedHighFst$dxy_Circumboreal_Asia_Circumboreal_Euro
sharedHighFst$dxy_Circumboreal_Euro_Circumboreal_NorAm

length(sharedHighFst$ratio)
dim(sharedHighFst)

sharedHighFst<-normalFST[!is.infinite(sharedHighFst$ratio)|is.na(sharedHighFst$ratio),]
dim(sharedHighFst)
normalFST$ratio<-normalFST$dxy_Circumboreal_Asia_Circumboreal_Euro/normalFST$dxy_Circumboreal_Euro_Circumboreal_NorAm
length(normalFST$ratio)
normalFST<-normalFST[!is.infinite(normalFST$ratio)|is.na(normalFST$ratio),]

sharedHighFst
normalFST

# How many sites of shared high FST (>0.1) are there?
length(sharedHighFst$scaffold)

install.packages("vioplot")
library(vioplot)

# Compare stats at windows with high fst to windows with Fst<0.1
vioplot(sharedHighFst$dxy_Circumboreal_Asia_Circumboreal_Euro,normalFST$dxy_Circumboreal_Asia_Circumboreal_Euro)
vioplot(sharedHighFst$pi_Circumboreal_Asia,normalFST$pi_Circumboreal_Asia)
vioplot(sharedHighFst$pi_Circumboreal_Euro,normalFST$pi_Circumboreal_Euro)
vioplot(sharedHighFst$pi_Circumboreal_NorAm,normalFST$pi_Circumboreal_NorAm)

# Do shared high FST windows contain haplotypes predating the Victoria-Kivu split?
vioplot(sharedHighFst$ratio,normalFST$ratio,horizontal = T)

install.packages("beanplot")
require("beanplot")

par(mfrow=c(1,1),mar=c(1,3,1,1),mgp=c(1.6,0.5,0),cex.axis=0.9,cex=1.3,xaxs="i",yaxs="i")
bp<-beanplot(normalFST$ratio,sharedHighFst$ratio,log = "",
             side='both', border='NA',
             col=list('cornflowerblue','blue'),
             ylab='maximum relative divergence' ,what=c(0,1,0,0),xaxt="n",
             ylim=c(0,3),cex=2)
abline(h=1)
par(xpd=F)

# Function to plot the density of dots in a grid
colPlot <- function(dataset=stats,varx,vary,minx=0,miny=0,maxx=0.02,maxy=maxx,title="",xlab=varx,ylab=vary,corr=T){
  rcbpal<-c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026")
  data<-dataset
  rf<-colorRampPalette((rcbpal))
  r<-rf(12)
  xbin<-seq(minx,maxx,length.out = 100)
  ybin<-seq(miny,maxy,length.out = 100)
  # colx=grep(names(data),pattern=varx)
  # coly=grep(names(data),pattern=vary)
  colx=varx
  coly=vary
  freq<-as.data.frame(table(findInterval(x = data[,colx],vec = xbin,all.inside=T),findInterval(x = data[,coly],vec = ybin,all.inside = T)))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))
  freq2D<-matrix(0,nrow=length(xbin),ncol=length(ybin))
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  image(xbin,ybin,log10(freq2D),col=r,breaks = seq(minx,max(log10(freq2D)),length.out=length(r)+1),
        xlab=xlab,ylab=ylab,xlim=c(minx,maxx),ylim=c(miny,maxy))
  barval<-round(10^(seq(minx,max(log10(freq2D)),length.out = length(r)+1)))
  title(main=title)
  lm<-summary(lm(data[,coly]~data[,colx]))
  r2<-round(lm$r.squared,3)
  if(corr) legend("top",legend=bquote(r^2 == .(r2)),bty = "n",cex=1.2)
  abline(lm(data[,coly]~data[,colx]))
}

par(mfrow=c(2,2),mar=c(5,5,0,0))

# Are the FST values of the two species pairs correlated?
colPlot(dataset = windowStats,varx = "Fst_Circumboreal_Asia_Circumboreal_Euro",vary = "Fst_Circumboreal_Euro_Circumboreal_NorAm",maxx = 1,maxy=1)
colPlot(dataset = windowStats,varx = "dxy_Circumboreal_Asia_Circumboreal_Euro",vary = "dxy_Circumboreal_Euro_Circumboreal_NorAm",maxx = 0.03,maxy=0.03)

# Are dxy and pi correlated?
colPlot(dataset = windowStats,varx = "pi_Circumboreal_Asia",vary = "dxy_Circumboreal_Asia_Circumboreal_Euro",maxx = 0.03,maxy=0.03)
colPlot(dataset = windowStats,varx = "pi_Circumboreal_Euro",vary = "dxy_Circumboreal_Asia_Circumboreal_Euro",maxx = 0.03,maxy=0.03)


