# PAL1314_Exp13_heatmap.R

# Purpose: Generate heatmap for data from PAL1314 Exp13

# ****** Assumes all variables created by PAL1314_liposome_expts.R are already
# in user's workspace; this file can be found at
# https://github.com/jamesrco/LipidPhotoOxBox

# Created 10/24/16 by J.R.C.

# necessary libraries

library(gplots)

# pull data out of Exp_13_PC.samp.pmol.mL.norm.mean, Exp_13_FFA.neg.samp.norm.mean,
# Exp_13_LPC.samp.pmol.mL.neg.mean (from PAL1314_liposome_expts.R) to prepare a 
# matrix for heatmap plotting; need just the means, not the sd or se's

# note: also need the data frames Exp_13_PC.samp and Exp_13_neg.subset later on

Exp_13_PC.samp.pmol.mL.norm.HM = Exp_13_PC.samp.pmol.mL.norm.mean[,seq(1, 33, 3)]
Exp_13_FFA.neg.samp.pmol.mL.mean.HM = Exp_13_FFA.neg.samp.pmol.mL.mean[,seq(1, 33, 3)]
Exp_13_LPC.samp.pmol.mL.neg.mean.HM = Exp_13_LPC.samp.pmol.mL.neg.mean[,seq(1, 33, 3)]

# generate some good colnames and rownames for the heatmaps, consistent w/presentation
# in other figures and tables

colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("Dark_control_no_HB"),c("Dark control"),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("Dark_control_plus_HB"),c("Dark control +het. bact."),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("Quartz_no_HB"),c("+UVB, -het. bact."),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("Quartz_plus_HB"),c("+UVB, +het. bact."),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("EPA_no_HB"),c("-UVB, -het. bact."),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))

colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("_2013-12-14 09:30:00.mean"),c(", 0 h"),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("_2013-12-14 13:40:00.mean"),c(", + 4 h"),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))
colnames(Exp_13_PC.samp.pmol.mL.norm.HM) = 
  gsub(c("_2013-12-14 17:50:00.mean"),c(", + 8 h"),colnames(Exp_13_PC.samp.pmol.mL.norm.HM))

colnames(Exp_13_FFA.neg.samp.pmol.mL.mean.HM) = colnames(Exp_13_PC.samp.pmol.mL.norm.HM)
colnames(Exp_13_LPC.samp.pmol.mL.neg.mean.HM) = colnames(Exp_13_PC.samp.pmol.mL.norm.HM)

# swap bulk acyl C numbers for more specific notation

rownames(Exp_13_PC.samp.pmol.mL.norm.HM) =
  gsub(c("^PC 44:12"),c("PC 22:6, 22:6"),rownames(Exp_13_PC.samp.pmol.mL.norm.HM))
rownames(Exp_13_PC.samp.pmol.mL.norm.HM) =
  gsub(c("^PC 36:2"),c("PC 18:1, 18:1"),rownames(Exp_13_PC.samp.pmol.mL.norm.HM))
rownames(Exp_13_PC.samp.pmol.mL.norm.HM) =
  gsub(c("^PC 32:0"),c("PC 16:0, 16:0"),rownames(Exp_13_PC.samp.pmol.mL.norm.HM))
rownames(Exp_13_PC.samp.pmol.mL.norm.HM) =
  gsub(c("^PC 36:0"),c("PC 18:0, 18:0"),rownames(Exp_13_PC.samp.pmol.mL.norm.HM))

# append RT data to row names

rownames(Exp_13_PC.samp.pmol.mL.norm.HM) = apply(cbind(rownames(Exp_13_PC.samp.pmol.mL.norm.HM), as.character(round(Exp_13_PC.samp$peakgroup.rt/60,1))),1,paste,collapse = ", RT ")
rownames(Exp_13_PC.samp.pmol.mL.norm.HM) = apply(cbind(rownames(Exp_13_PC.samp.pmol.mL.norm.HM),rep("min.",13)),1,paste,collapse = " ")

# **** for FFA and LPC, taking a manual shortcut (bad bad bad)
rownames(Exp_13_FFA.neg.samp.pmol.mL.mean.HM) = apply(cbind(rownames(Exp_13_FFA.neg.samp.pmol.mL.mean.HM), as.character(round(Exp_13_neg.subset$peakgroup_rt/60,1)[c(2,3,5:9,12:18,20)])),1,paste,collapse = ", RT ")
rownames(Exp_13_FFA.neg.samp.pmol.mL.mean.HM) = apply(cbind(rownames(Exp_13_FFA.neg.samp.pmol.mL.mean.HM),rep("min.",15)),1,paste,collapse = " ")

rownames(Exp_13_LPC.samp.pmol.mL.neg.mean.HM) = apply(cbind(rownames(Exp_13_LPC.samp.pmol.mL.neg.mean.HM), as.character(round(Exp_13_neg.subset$peakgroup_rt/60,1)[c(4,10,11,19)])),1,paste,collapse = ", RT ")
rownames(Exp_13_LPC.samp.pmol.mL.neg.mean.HM) = apply(cbind(rownames(Exp_13_LPC.samp.pmol.mL.neg.mean.HM),rep("min.",4)),1,paste,collapse = " ")

# change to matrices

Exp_13_PC.samp.pmol.mL.norm.HM = as.matrix(Exp_13_PC.samp.pmol.mL.norm.HM)
Exp_13_FFA.neg.samp.pmol.mL.mean.HM = as.matrix(Exp_13_FFA.neg.samp.pmol.mL.mean.HM)
Exp_13_LPC.samp.pmol.mL.neg.mean.HM = as.matrix(Exp_13_LPC.samp.pmol.mL.neg.mean.HM)

# # fix NA's and small values
# 
# screenedpeaks_exptmeans.full_PC.byttp[screenedpeaks_exptmeans.full_PC.byttp<1000] <- 1000
# screenedpeaks_exptmeans.full_PC.byttp[1,is.na(screenedpeaks_exptmeans.full_PC.byttp[1,])] <- 10 # just need to do this for the first row

# calculate fold-change relative to initial

Exp_13_PC.samp.pmol.mL.norm.HM.foldchange <- sweep(Exp_13_PC.samp.pmol.mL.norm.HM, 1, Exp_13_PC.samp.pmol.mL.norm.HM[,1], "/")
Exp_13_FFA.neg.samp.pmol.mL.mean.HM.foldchange <- sweep(Exp_13_FFA.neg.samp.pmol.mL.mean.HM, 1, Exp_13_FFA.neg.samp.pmol.mL.mean.HM[,1], "/")
Exp_13_LPC.samp.pmol.mL.neg.mean.HM.foldchange <- sweep(Exp_13_LPC.samp.pmol.mL.neg.mean.HM, 1, Exp_13_LPC.samp.pmol.mL.neg.mean.HM[,1], "/")

# # return any initial values (first row) that were originally NA's back to NA's
# 
# screenedpeaks_exptmeans.full_PC.byttp[1,screenedpeaks_exptmeans.full_PC.byttp[1,]==10] <- NA 

# transform into log2 fold change

Exp_13_PC.samp.pmol.mL.norm.HM.foldchange.log2 <- log2(Exp_13_PC.samp.pmol.mL.norm.HM.foldchange)
Exp_13_FFA.neg.samp.pmol.mL.mean.HM.foldchange.log2 <- log2(Exp_13_FFA.neg.samp.pmol.mL.mean.HM.foldchange)
Exp_13_LPC.samp.pmol.mL.neg.mean.HM.foldchange.log2 <- log2(Exp_13_LPC.samp.pmol.mL.neg.mean.HM.foldchange)

# screenedpeaks_exptmeans.full_PC.byttp.rev <- screenedpeaks_exptmeans.full_PC.byttp[dim(screenedpeaks_exptmeans.full_PC.byttp):1,]

# build heatmaps

# first, combine all features into single matrix

Exp_13_features.HM.foldchange.log2 = rbind(Exp_13_PC.samp.pmol.mL.norm.HM.foldchange.log2,
                                           Exp_13_FFA.neg.samp.pmol.mL.mean.HM.foldchange.log2,
      Exp_13_LPC.samp.pmol.mL.neg.mean.HM.foldchange.log2)

# generate breaks and colors

breaks.neg = seq(min(Exp_13_features.HM.foldchange.log2, na.rm=T),0,length.out=50)
breaks.pos = seq(0,max(Exp_13_features.HM.foldchange.log2, na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-50] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

# full heatmap

heatmap.2(Exp_13_features.HM.foldchange.log2[,2:ncol(Exp_13_features.HM.foldchange.log2)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(12,12))

# maybe just PC 22:6 and derivatives

# create our subset
Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2[
  grepl(paste(c("PC 22:6",
                "FFA 22:6",
                "LPC 22:6"), collapse = "|"), rownames(Exp_13_features.HM.foldchange.log2)),]

# also, just want the -HB controls (not the +HB controls)
Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2.sub226[,-c(4:5)]

# lastly, get rid of a duplicate 22:6 FFA (ID did not stand up to manual examination of spectra)
Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2.sub226[-c(10),]

# create a matrix (yes, manually based on inspection of output from
# PAL1314_liposome_expts.R - I know, inefficient) of p-values for each block in the
# subsetted heatmap

Exp_13_features.HM.foldchange.log2.sub226.pvals = 
  matrix(nrow = nrow(Exp_13_features.HM.foldchange.log2.sub226),
         ncol = ncol(Exp_13_features.HM.foldchange.log2.sub226)-1)

# manually populate based on Tukey HSD test output from PAL1314_liposome_expts.R
# will use "1" as not significant indicator
Exp_13_features.HM.foldchange.log2.sub226.pvals[1,] =
  c(1,1,1,0.05,1,0.01,1,0.01) # PC 22:6, 22:6
Exp_13_features.HM.foldchange.log2.sub226.pvals[2,] =
  c(1,1,1,0.01,1,0.01,1,0.05) # PC 44:12 +2O
Exp_13_features.HM.foldchange.log2.sub226.pvals[3,] =
  c(1,1,1,0.05,1,0.01,1,1) # PC 44:12 +4O
Exp_13_features.HM.foldchange.log2.sub226.pvals[4,] =
  c(rep(1,8)) # PC 44:12 +1O
Exp_13_features.HM.foldchange.log2.sub226.pvals[5,] =
  c(rep(1,8)) # PC 44:12 +3O
Exp_13_features.HM.foldchange.log2.sub226.pvals[6,] =
  c(rep(1,8)) # PC 44:12 +3O
Exp_13_features.HM.foldchange.log2.sub226.pvals[7,] =
  c(rep(1,8)) # FFA 22:6
Exp_13_features.HM.foldchange.log2.sub226.pvals[8,] =
  c(1,1,1,1,1,0.0001,1,0.05) # FFA 22:6 +2O
Exp_13_features.HM.foldchange.log2.sub226.pvals[9,] =
  c(1,1,1,1,1,0.0001,1,1) # FFA 22:6 +3O
Exp_13_features.HM.foldchange.log2.sub226.pvals[10,] =
  c(1,1,1,1,1,0.01,1,1) # FFA 22:6 +10
Exp_13_features.HM.foldchange.log2.sub226.pvals[11,] =
  c(1,1,1,1,1,0.001,1,1) # LPC 22:6 +4O
Exp_13_features.HM.foldchange.log2.sub226.pvals[12,] =
  c(1,1,1,0.05,1,0.001,1,0.05) # LPC 22:6 +2O
Exp_13_features.HM.foldchange.log2.sub226.pvals[13,] =
  c(1,1,1,1,1,0.05,1,1) # LPC 22:6
Exp_13_features.HM.foldchange.log2.sub226.pvals[14,] =
  c(1,1,1,1,1,0.01,1,1) # LPC 22:6 +1O

# create a matrix of actual symbols to be plotted

Exp_13_features.HM.sigsymbols = Exp_13_features.HM.foldchange.log2.sub226.pvals
Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==1]=""
Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.05]="+"
Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.01]="*"
Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.0001]="***"

par(oma=c(0,0,0,0)) # set margins

pdf(file = "Exp13_heatmap_PC22-6plus.pdf",
    width = 7, height = 9.5, pointsize = 12,
    bg = "white")

#par(mar=c(5,5,1,1))

heatmap.2(Exp_13_features.HM.foldchange.log2.sub226[,2:ncol(Exp_13_features.HM.foldchange.log2.sub226)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(16,16),cellnote=Exp_13_features.HM.sigsymbols)

dev.off()

# an alternate plot without the -UVB treatment samples

par(oma=c(0,0,0,0)) # set margins

pdf(file = "Exp13_heatmap_PC22-6plus.pdf",
    width = 7, height = 9.5, pointsize = 12,
    bg = "white")

#par(mar=c(5,5,1,1))

heatmap.2(Exp_13_features.HM.foldchange.log2.sub226[,c(2:3,6:9)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(16,16),cellnote=Exp_13_features.HM.sigsymbols[,c(1:2,5:8)])

dev.off()

# now, a small heatmap of just the other unoxidized parent lipids in this same
# experiment, for comparison purposes

Exp_13_features.HM.foldchange.log2.subparents = Exp_13_features.HM.foldchange.log2[
  grepl(paste(c("PC 18:1, 18:1,",
                "PC 18:0, 18:0,",
                "PC 16:0, 16:0,"), collapse = "|"), rownames(Exp_13_features.HM.foldchange.log2)),]

# get rid of the +HB controls
Exp_13_features.HM.foldchange.log2.subparents = Exp_13_features.HM.foldchange.log2.subparents[,-c(4:5)]

# will need to order the treatments manually to fit with dendrogram clustering
# of bigger heatmap

Exp_13_features.HM.foldchange.log2.subparents = 
  Exp_13_features.HM.foldchange.log2.subparents[,c(2,3,8,6,4,5,7,9)]

# put species in a more logical order

Exp_13_features.HM.foldchange.log2.subparents = 
  Exp_13_features.HM.foldchange.log2.subparents[c(1,3,2),]

par(oma=c(0,0,0,0)) # set margins

pdf(file = "Exp13_heatmap_other_parents.pdf",
    width = 7, height = 10, pointsize = 12,
    bg = "white")

#par(mar=c(5,5,1,1))

heatmap.2(Exp_13_features.HM.foldchange.log2.subparents[,1:ncol(Exp_13_features.HM.foldchange.log2.subparents)],breaks=breaks.set,col=hm.colors,scale="none",Colv=FALSE,trace="none",Rowv=FALSE,dendrogram="none",na.color="grey",key=T,
          density.info=c("none"),margins=c(16,16))

dev.off()

