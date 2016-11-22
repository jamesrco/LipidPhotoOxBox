# PAL1314_liposome_expts.R

# Created 4/28/2016 by J.R.C.; initial commit to Github 9/23/16

# For workup and analysis of data from liposome experiments from PAL1314 field season

### workspace preparation ####

# load necessary libraries

library(LOBSTAHS)
library(chemCal)

# set wd

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")

### load some necessary metadata; define functions ###

# load in metadata

meta.raw = read.csv("Exactive data - PAL 1314 experimental samples.csv")

# some chemical data

DNPPE_mg_mL_1314 = 0.0565 # concentration DNPPE used in liquid/liquid extractions
# during  PAL1314 Antarctic work, mg/mL
DNPPE_MW = 875.081 # MW DNPPE, g/mol

# define a two-step function "splitpred" to compute pmol o.c. from raw values for PC and DNPPE, using one of two linear standard curves depending on magnitude of concentration (defined by cutoff value)

splitpred = function(x,linfit_low,linfit_hi,cutoff) {
  
  input = x
  
  if (input > cutoff) {
    
    predval = inverse.predict(linfit_hi, newdata = x)
    
  } else if (input <= cutoff) {
    
    predval = inverse.predict(linfit_low, newdata = x)
    
  }
  
  predval$Prediction
  
}

# define function allow "easy" retrieval of sample metadata based on sample ID in filename

getMetDat = function(fn,metadat.raw,whichdat) {
  
  # fn is sample descriptor (from verbose file name) from which sample ID will be extracted and then matched
  # whichdat is a list of metadata to be retrieved (corresponding to name of column in the metadat.raw table)
  
  match_ind = grepl(sub("^.*_","\\1",fn),metadat.raw$Orbi.sequence.ID) # get index of matching metadata

  this.metdat = metadat.raw[match_ind,whichdat] # extract called-for metadata
  
  this.metdat # return extracted data
  
}

### standards ###

### + mode standards  ###

# initial set, from 20160421

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/IPL_Standards_20160421_pos.Rdata") # load processed data
IPLstd_pos.raw = getLOBpeaklist(IPL_Standards_20160421_pos) # generate peaklist

# extract only the PC standards and DNPPE
PCstds = IPLstd_pos.raw[IPLstd_pos.raw$compound_name=="PC 32:0",]
PCstds = PCstds[14:26]
PCstds = PCstds[order(colnames(PCstds))]
DNPPEstds = IPLstd_pos.raw[IPLstd_pos.raw$compound_name=="DNPPE",]
DNPPEstds = DNPPEstds[14:26]
DNPPEstds = DNPPEstds[order(colnames(DNPPEstds))]

# separate 2 QC's from the standards
PCstds_QC = PCstds[c(1,length(PCstds))]
PCstds = PCstds[-c(1,length(PCstds))]
DNPPEstds_QC = DNPPEstds[c(1,length(DNPPEstds))]
DNPPEstds = DNPPEstds[-c(1,length(DNPPEstds))]

# fit standard curves using standard data

# necessary vectors 
PCoc = c(0.246,0.493,0.986,1.971,3.943,7.885,15.770,31.541,63.082,126.163,252.326) # level (pmol o.c.) of PC standards (assuming 20 uL injection, from HFF spreadsheet; MGDG at ~ 16k pmol)
DNPPEoc_20160331 = c(0.4553,0.9106,1.8212,3.6425,7.285,14.570,29.14025,58.2805,116.561,0.000,0.000) # level (pmol o.c.) of DNPPE standards (assuming 20 uL injection, based on 4k standard at 0.051 mg/mL DNPPE; this is DNPPE conc. from aliquots dated 3/31/16) 

# PC

# curve fitting & diagnostics

y = PCstds[1:9]
x = PCoc[1:9]
linfit_low.PC = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(PCoc,PCstds,pch="+")
points(PCoc[1:9],fitted(linfit_low.PC),col="red",pch="+")

# we will need some other fit for levels higher than ~ 100 pmol o.c. (not unexpected)
y = PCstds
x = PCoc
linfit_hi.PC = lm(as.numeric(y)~x)
points(PCoc,fitted(linfit_hi.PC),col="blue",pch="+")

# DNPPE

# curve fitting & diagnostics

y = DNPPEstds[1:7]
x = DNPPEoc_20160331[1:7]
linfit_low.DNPPE = lm(as.numeric(y)~x) # fit a first linear model for the first 7 standard levels
plot(DNPPEoc_20160331,DNPPEstds,pch="+")
points(DNPPEoc_20160331[1:7],fitted(linfit_low.DNPPE),col="red",pch="+")

y = DNPPEstds
x = DNPPEoc_20160331
linfit_hi.DNPPE = lm(as.numeric(y)~x) # fit other linear model for higher concentrations
points(DNPPEoc_20160331,fitted(linfit_hi.DNPPE),col="blue",pch="+")

# second set of standards, from 20161107 (needed for reanalysis of Exp_03a, and
# analysis of particulate environmental data)

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/6IPL_Standards_20161107_pos.RData") # load processed data
IPLstd_pos_20161107.raw = getLOBpeaklist(x6IPL_Standards_20161107_pos) # generate peaklist

# extract only the PC standards and DNPPE
PCstds.20161107 = IPLstd_pos_20161107.raw[IPLstd_pos_20161107.raw$compound_name=="PC 32:0",]
PCstds.20161107 = PCstds.20161107[14:25]
PCstds.20161107 = PCstds.20161107[order(colnames(PCstds.20161107))]
DNPPEstds.20161107 = IPLstd_pos_20161107.raw[IPLstd_pos_20161107.raw$compound_name=="DNPPE",]
DNPPEstds.20161107 = DNPPEstds.20161107[14:25]
DNPPEstds.20161107 = DNPPEstds.20161107[order(colnames(DNPPEstds.20161107))]

# separate a QC from the standards
PCstds_QC.20161107 = PCstds.20161107[c(length(PCstds.20161107))]
PCstds.20161107 = PCstds.20161107[-c(length(PCstds.20161107))]
DNPPEstds_QC.20161107 = DNPPEstds.20161107[c(length(DNPPEstds.20161107))]
DNPPEstds.20161107 = DNPPEstds.20161107[-c(length(DNPPEstds.20161107))]

# fit standard curves using standard data

# PC

# curve fitting & diagnostics

y = PCstds.20161107[1:9]
x = PCoc[1:9]
linfit_low.PC.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(PCoc,PCstds.20161107,pch="+")
points(PCoc[1:9],fitted(linfit_low.PC.20161107),col="red",pch="+")

# we will need some other fit for levels higher than ~ 100 pmol o.c. (not unexpected)
y = PCstds.20161107
x = PCoc
linfit_hi.PC.20161107 = lm(as.numeric(y)~x)
points(PCoc,fitted(linfit_hi.PC.20161107),col="blue",pch="+")

# DNPPE

# curve fitting & diagnostics

y = DNPPEstds.20161107[1:7]
x = DNPPEoc_20160331[1:7]
linfit_low.DNPPE.20161107 = lm(as.numeric(y)~x) # fit a first linear model for the first 7 standard levels
plot(DNPPEoc_20160331,DNPPEstds.20161107,pch="+")
points(DNPPEoc_20160331[1:7],fitted(linfit_low.DNPPE.20161107),col="red",pch="+")

y = DNPPEstds.20161107
x = DNPPEoc_20160331
linfit_hi.DNPPE.20161107 = lm(as.numeric(y)~x) # fit other linear model for higher concentrations
points(DNPPEoc_20160331,fitted(linfit_hi.DNPPE.20161107),col="blue",pch="+")

### - mode standards  ###

# initial set, from 20160421

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/IPL_Standards_20160421_neg.Rdata") # load processed data
IPLstd_neg.raw = getLOBpeaklist(IPL_Standards_20160421_neg) # generate peaklist

# extract only the PC standards and DNPPE
PCstds.neg = IPLstd_neg.raw[IPLstd_neg.raw$compound_name=="PC 32:0",]
PCstds.neg = PCstds.neg[14:26]
PCstds.neg = PCstds.neg[order(colnames(PCstds.neg))]
DNPPEstds.neg = IPLstd_neg.raw[IPLstd_neg.raw$compound_name=="DNPPE",]
DNPPEstds.neg = DNPPEstds.neg[14:26]
DNPPEstds.neg = DNPPEstds.neg[order(colnames(DNPPEstds.neg))]

# separate 2 QC's from the standards
PCstds_QC.neg = PCstds.neg[c(1,length(PCstds.neg))]
PCstds.neg = PCstds.neg[-c(1,length(PCstds.neg))]
DNPPEstds_QC.neg = DNPPEstds.neg[c(1,length(DNPPEstds.neg))]
DNPPEstds.neg = DNPPEstds.neg[-c(1,length(DNPPEstds.neg))]

# fit standard curves using standard data

# PC

# curve fitting & diagnostics

y = PCstds.neg[1:9]
x = PCoc[1:9]
linfit_low.PC.neg = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(PCoc,PCstds.neg,pch="+")
points(PCoc[1:9],fitted(linfit_low.PC.neg),col="red",pch="+")

# we will need some other fit for levels higher than ~ 100 pmol o.c. (not unexpected)
y = PCstds.neg
x = PCoc
linfit_hi.PC.neg = lm(as.numeric(y)~x)
points(PCoc,fitted(linfit_hi.PC.neg),col="blue",pch="+")

# DNPPE

# curve fitting & diagnostics

y = DNPPEstds.neg[1:7]
x = DNPPEoc_20160331[1:7]
linfit_low.DNPPE.neg = lm(as.numeric(y)~x) # fit a first linear model for the first 7 standard levels
plot(DNPPEoc_20160331,DNPPEstds.neg,pch="+")
points(DNPPEoc_20160331[1:7],fitted(linfit_low.DNPPE.neg),col="red",pch="+")

y = DNPPEstds.neg
x = DNPPEoc_20160331
linfit_hi.DNPPE.neg = lm(as.numeric(y)~x) # fit other linear model for higher concentrations
points(DNPPEoc_20160331,fitted(linfit_hi.DNPPE.neg),col="blue",pch="+")

# DHA standards, run on 20161109
# needed for analysis of Exp_13 data

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/DHA_Standards_20161109_neg.RData") # load processed data
DHA_std.raw = getLOBpeaklist(DHA_Standards_20161109_neg) # generate peaklist

# extract only the DHA standards
DHAstds.neg = DHA_std.raw[DHA_std.raw$compound_name=="FFA 22:6",]
DHAstds.neg = DHAstds.neg[13:23]
DHAstds.neg = DHAstds.neg[order(colnames(DHAstds.neg))]

# fit standard curves using standard data

# curve fitting & diagnostics

# calculate levels (pmol o.c.) of DHA standards used (assuming 20 uL injection)

# specify highest concentration standard used 
DHA_20161109_pmol_mL_highest_pmol_mL = 14323.2

# calculate amount in 20 uL injection
DHA_20161109_pmol_oc_highest = DHA_20161109_pmol_mL_highest_pmol_mL*.02

# preallocate vector
DHAoc = rep(NA,11)
DHAoc[1] = DHA_20161109_pmol_oc_highest

# calculate conc's in serial dilution
for (i in 2:length((DHAoc))) {
  DHAoc[i] = DHAoc[i-1]/2
}

DHAoc = rev(DHAoc)

y = DHAstds.neg
x = DHAoc
linfit.DHA.neg = lm(as.numeric(y)~x-1) # fit a linear model, force through origin
plot(DHAoc,DHAstds.neg,pch="+")
points(DHAoc,fitted(linfit.DHA.neg),col="red",pch="+")

### load in & perform initial analysis of experimental data ###

### positive ion mode ###

# Experiment 1

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_01_pos.RData")
Exp_01_pos.raw = getLOBpeaklist(Exp_01_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_01_pos.subset = Exp_01_pos.raw[(Exp_01_pos.raw$species==c("PC") & (
  (Exp_01_pos.raw$FA_total_no_C==32 &
     Exp_01_pos.raw$FA_total_no_DB==0 &
     Exp_01_pos.raw$degree_oxidation==0) |
    (Exp_01_pos.raw$FA_total_no_C==32 &
       Exp_01_pos.raw$FA_total_no_DB==2 &
       Exp_01_pos.raw$degree_oxidation>=0) |
    (Exp_01_pos.raw$FA_total_no_C==36 &
       Exp_01_pos.raw$FA_total_no_DB<=2 &
       Exp_01_pos.raw$degree_oxidation>=0))) |
    Exp_01_pos.raw$species==c("DNPPE")
  ,]

Exp_01_PC = Exp_01_pos.subset[grep("PC",Exp_01_pos.subset$compound_name),c(2,16:49)]
Exp_01_DNPPE = Exp_01_pos.subset[grep("DNPPE",Exp_01_pos.subset$compound_name),c(2,16:49)]

# extract QCs

Exp_01_PC.QC = Exp_01_PC[,c(1,grep("QC_",colnames(Exp_01_PC)))]
Exp_01_PC.QC = Exp_01_PC.QC[Exp_01_PC.QC$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_01_DNPPE.QC = Exp_01_DNPPE[,c(1,grep("QC_",colnames(Exp_01_DNPPE)))]

# eliminate QCs from list of data

Exp_01_PC.samp = Exp_01_PC[,-c(grep("QC_",colnames(Exp_01_PC)))]
Exp_01_DNPPE.samp = Exp_01_DNPPE[,-c(grep("QC_",colnames(Exp_01_DNPPE)))]

# calculate concentrations

# pmol o.c.

Exp_01_PC.samp.pmol_oc = apply(Exp_01_PC.samp[,2:ncol(Exp_01_PC.samp)],c(1,2),splitpred,linfit_low.PC,linfit_hi.PC,1e10)
Exp_01_DNPPE.samp.pmol_oc = apply(Exp_01_DNPPE.samp[,2:ncol(Exp_01_DNPPE.samp)],c(1,2),splitpred,linfit_low.DNPPE,linfit_hi.DNPPE,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_01_PC.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_01_PC.samp.pmol_oc),getMetDat,meta.raw,c(2:7,14)))
Exp_01_PC.metdat$Date.time.sample.collected = strptime(as.character(Exp_01_PC.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_01_PC.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_01_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Exp_01_DNPPE.samp.pmol_oc # recovery factor
Exp_01_PC.samp.pmol.total = sweep(Exp_01_PC.samp.pmol_oc, 2, Exp_01_DNPPE.samp.RF, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_01_PC.samp.pmol.mL = sweep(Exp_01_PC.samp.pmol.total, 2, Exp_01_PC.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# Experiment 2a

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_02a_pos_new.RData")
Exp_02a_pos.raw = getLOBpeaklist(Exp_02a_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_02a_pos.subset = Exp_02a_pos.raw[(Exp_02a_pos.raw$species==c("PC") & (
  (Exp_02a_pos.raw$FA_total_no_C==32 &
     Exp_02a_pos.raw$FA_total_no_DB==0 &
     Exp_02a_pos.raw$degree_oxidation==0) |
    (Exp_02a_pos.raw$FA_total_no_C==32 &
       Exp_02a_pos.raw$FA_total_no_DB==2 &
       Exp_02a_pos.raw$degree_oxidation>=0) |
    (Exp_02a_pos.raw$FA_total_no_C==36 &
       Exp_02a_pos.raw$FA_total_no_DB<=4 &
       Exp_02a_pos.raw$degree_oxidation>=0))) |
    Exp_02a_pos.raw$species==c("DNPPE")
  ,]

Exp_02a_PC = Exp_02a_pos.subset[grep("PC",Exp_02a_pos.subset$compound_name),c(2,16:51)]
Exp_02a_DNPPE = Exp_02a_pos.subset[grep("DNPPE",Exp_02a_pos.subset$compound_name),c(2,16:51)]

# extract QCs

Exp_02a_PC.QC = Exp_02a_PC[,c(1,grep("QC_",colnames(Exp_02a_PC)))]
Exp_02a_PC.QC = Exp_02a_PC.QC[Exp_02a_PC.QC$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_02a_DNPPE.QC = Exp_02a_DNPPE[,c(1,grep("QC_",colnames(Exp_02a_DNPPE)))]

# eliminate QCs from list of data

Exp_02a_PC.samp = Exp_02a_PC[,-c(grep("QC_",colnames(Exp_02a_PC)))]
Exp_02a_DNPPE.samp = Exp_02a_DNPPE[,-c(grep("QC_",colnames(Exp_02a_DNPPE)))]

# calculate concentrations

# pmol o.c.

Exp_02a_PC.samp.pmol_oc = apply(Exp_02a_PC.samp[,2:ncol(Exp_02a_PC.samp)],c(1,2),splitpred,linfit_low.PC,linfit_hi.PC,1e10)
Exp_02a_DNPPE.samp.pmol_oc = apply(Exp_02a_DNPPE.samp[,2:ncol(Exp_02a_DNPPE.samp)],c(1,2),splitpred,linfit_low.DNPPE,linfit_hi.DNPPE,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_02a_PC.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_02a_PC.samp.pmol_oc),getMetDat,meta.raw,c(2:7,14)))
Exp_02a_PC.metdat$Date.time.sample.collected = strptime(as.character(Exp_02a_PC.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_02a_PC.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_02a_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Exp_02a_DNPPE.samp.pmol_oc # recovery factor
Exp_02a_PC.samp.pmol.total = sweep(Exp_02a_PC.samp.pmol_oc, 2, Exp_02a_DNPPE.samp.RF, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_02a_PC.samp.pmol.mL = sweep(Exp_02a_PC.samp.pmol.total, 2, Exp_02a_PC.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# Experiment 3a, first run of samples 

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_03a_pos_new.RData")
Exp_03a_pos.raw = getLOBpeaklist(Exp_03a_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_03a_pos.subset = Exp_03a_pos.raw[(Exp_03a_pos.raw$species==c("PC") & (
  (Exp_03a_pos.raw$FA_total_no_C==32 &
     Exp_03a_pos.raw$FA_total_no_DB==0 &
     Exp_03a_pos.raw$degree_oxidation==0) |
    (Exp_03a_pos.raw$FA_total_no_C==36 &
       Exp_03a_pos.raw$FA_total_no_DB==0 &
       Exp_03a_pos.raw$degree_oxidation==0) |
    (Exp_03a_pos.raw$FA_total_no_C==44 &
       Exp_03a_pos.raw$FA_total_no_DB==0 &
       Exp_03a_pos.raw$degree_oxidation==0) |
    (Exp_03a_pos.raw$FA_total_no_C==44 &
       Exp_03a_pos.raw$FA_total_no_DB %in% c(10,11,12) &
       Exp_03a_pos.raw$degree_oxidation>=0))) |
    Exp_03a_pos.raw$species==c("DNPPE")
  ,]

Exp_03a_PC = Exp_03a_pos.subset[grep("PC",Exp_03a_pos.subset$compound_name),c(2,16:51)]
Exp_03a_DNPPE = Exp_03a_pos.subset[grep("DNPPE",Exp_03a_pos.subset$compound_name),c(2,16:51)]
Exp_03a_DNPPE = Exp_03a_DNPPE[1,] # the less stringent criteria for screening of the Exp 3 data allowed secondary peakgroups falsely ID'd as DNPPE & PC 32:0 to get through
Exp_03a_PC = Exp_03a_PC[-4,] # the less stringent criteria for screening of the Exp 3 data allowed secondary peakgroups falsely ID'd as DNPPE & PC 32:0 to get through

# extract QCs

Exp_03a_PC.QC = Exp_03a_PC[,c(1,grep("QC_",colnames(Exp_03a_PC)))]
Exp_03a_PC.QC = Exp_03a_PC.QC[Exp_03a_PC.QC$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_03a_DNPPE.QC = Exp_03a_DNPPE[,c(1,grep("QC_",colnames(Exp_03a_DNPPE)))]

# eliminate QCs from list of data

Exp_03a_PC.samp = Exp_03a_PC[,-c(grep("QC_",colnames(Exp_03a_PC)))]
Exp_03a_DNPPE.samp = Exp_03a_DNPPE[,-c(grep("QC_",colnames(Exp_03a_DNPPE)))]

# calculate concentrations

# pmol o.c.

Exp_03a_PC.samp.pmol_oc = apply(Exp_03a_PC.samp[,2:ncol(Exp_03a_PC.samp)],c(1,2),splitpred,linfit_low.PC,linfit_hi.PC,1e10)
Exp_03a_DNPPE.samp.pmol_oc = apply(Exp_03a_DNPPE.samp[,2:ncol(Exp_03a_DNPPE.samp)],c(1,2),splitpred,linfit_low.DNPPE,linfit_hi.DNPPE,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_03a_PC.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_03a_PC.samp.pmol_oc),getMetDat,meta.raw,c(2:7,14)))
Exp_03a_PC.metdat$Date.time.sample.collected = strptime(as.character(Exp_03a_PC.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_03a_PC.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_03a_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Exp_03a_DNPPE.samp.pmol_oc # recovery factor
Exp_03a_PC.samp.pmol.total = sweep(Exp_03a_PC.samp.pmol_oc, 2, Exp_03a_DNPPE.samp.RF, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_03a_PC.samp.pmol.mL = sweep(Exp_03a_PC.samp.pmol.total, 2, Exp_03a_PC.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# Experiment 3a, second run of samples 

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_03a_rerunNov16_pos.RData")
Exp_03a_pos.rerun.raw = getLOBpeaklist(Exp_03a_rerunNov16_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_03a_pos.rerun.subset = Exp_03a_pos.rerun.raw[(Exp_03a_pos.rerun.raw$species==c("PC") & (
  (Exp_03a_pos.rerun.raw$FA_total_no_C==32 &
     Exp_03a_pos.rerun.raw$FA_total_no_DB==0 &
     Exp_03a_pos.rerun.raw$degree_oxidation==0) |
    (Exp_03a_pos.rerun.raw$FA_total_no_C==36 &
       Exp_03a_pos.rerun.raw$FA_total_no_DB==0 &
       Exp_03a_pos.rerun.raw$degree_oxidation==0) |
    (Exp_03a_pos.rerun.raw$FA_total_no_C==44 &
       Exp_03a_pos.rerun.raw$FA_total_no_DB==0 &
       Exp_03a_pos.rerun.raw$degree_oxidation==0) |
    (Exp_03a_pos.rerun.raw$FA_total_no_C==44 &
       Exp_03a_pos.rerun.raw$FA_total_no_DB %in% c(10,11,12) &
       Exp_03a_pos.rerun.raw$degree_oxidation>=0))) |
    Exp_03a_pos.rerun.raw$species==c("DNPPE")
  ,]

# eliminate a group incorrectly ID'd as PC 22:6, 22:6
Exp_03a_pos.rerun.subset = Exp_03a_pos.rerun.subset[-2,]

Exp_03a_PC.rerun = Exp_03a_pos.rerun.subset[grep("PC",Exp_03a_pos.rerun.subset$compound_name),c(2,16:51)]
Exp_03a_DNPPE.rerun = Exp_03a_pos.rerun.subset[grep("DNPPE",Exp_03a_pos.rerun.subset$compound_name),c(2,16:51)]

# extract QCs

Exp_03a_PC.QC.rerun = Exp_03a_PC.rerun[,c(1,grep("QC_",colnames(Exp_03a_PC.rerun)))]
Exp_03a_PC.QC.rerun = Exp_03a_PC.QC.rerun[Exp_03a_PC.QC.rerun$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_03a_DNPPE.QC.rerun = Exp_03a_DNPPE.rerun[,c(1,grep("QC_",colnames(Exp_03a_DNPPE.rerun)))]

# eliminate QCs from list of data

Exp_03a_PC.samp.rerun = Exp_03a_PC.rerun[,-c(grep("QC_",colnames(Exp_03a_PC.rerun)))]
Exp_03a_DNPPE.samp.rerun = Exp_03a_DNPPE.rerun[,-c(grep("QC_",colnames(Exp_03a_DNPPE.rerun)))]

# calculate concentrations

# pmol o.c.

Exp_03a_PC.samp.pmol_oc.rerun = apply(Exp_03a_PC.samp.rerun[,2:ncol(Exp_03a_PC.samp.rerun)],c(1,2),splitpred,linfit_low.PC.20161107,linfit_hi.PC.20161107,1e10)
Exp_03a_DNPPE.samp.pmol_oc.rerun = apply(Exp_03a_DNPPE.samp.rerun[,2:ncol(Exp_03a_DNPPE.samp.rerun)],c(1,2),splitpred,linfit_low.DNPPE.20161107,linfit_hi.DNPPE.20161107,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_03a_PC.metdat.rerun = do.call(rbind.data.frame,lapply(colnames(Exp_03a_PC.samp.pmol_oc.rerun),getMetDat,meta.raw,c(2:7,14)))
Exp_03a_PC.metdat.rerun$Date.time.sample.collected = strptime(as.character(Exp_03a_PC.metdat.rerun$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp.rerun = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_03a_PC.metdat.rerun$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_03a_DNPPE.samp.RF.rerun = DNPPE_pmol_added_per_samp.rerun/Exp_03a_DNPPE.samp.pmol_oc.rerun # recovery factor
Exp_03a_PC.samp.pmol.total.rerun = sweep(Exp_03a_PC.samp.pmol_oc.rerun, 2, Exp_03a_DNPPE.samp.RF.rerun, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_03a_PC.samp.pmol.mL.rerun = sweep(Exp_03a_PC.samp.pmol.total.rerun, 2, Exp_03a_PC.metdat.rerun$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# Experiment 12

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_12_pos.RData")
Exp_12_pos.raw = getLOBpeaklist(Exp_12_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_12_pos.subset = Exp_12_pos.raw[(Exp_12_pos.raw$species==c("PC") & (
  (Exp_12_pos.raw$FA_total_no_C==32 &
     Exp_12_pos.raw$FA_total_no_DB==0 &
     Exp_12_pos.raw$degree_oxidation==0) |
    (Exp_12_pos.raw$FA_total_no_C==32 &
       Exp_12_pos.raw$FA_total_no_DB==2 &
       Exp_12_pos.raw$degree_oxidation>=0) |
    (Exp_12_pos.raw$FA_total_no_C==36 &
       Exp_12_pos.raw$FA_total_no_DB<=4 &
       Exp_12_pos.raw$degree_oxidation>=0))) |
    Exp_12_pos.raw$species==c("DNPPE")
  ,]

Exp_12_PC = Exp_12_pos.subset[grep("PC",Exp_12_pos.subset$compound_name),c(2,18:53)]
Exp_12_DNPPE = Exp_12_pos.subset[grep("DNPPE",Exp_12_pos.subset$compound_name),c(2,18:53)]

# extract QCs

Exp_12_PC.QC = Exp_12_PC[,c(1,grep("QC_",colnames(Exp_12_PC)))]
Exp_12_PC.QC = Exp_12_PC.QC[Exp_12_PC.QC$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_12_DNPPE.QC = Exp_12_DNPPE[,c(1,grep("QC_",colnames(Exp_12_DNPPE)))]

# eliminate QCs from list of data

Exp_12_PC.samp = Exp_12_PC[,-c(grep("QC_",colnames(Exp_12_PC)))]
Exp_12_DNPPE.samp = Exp_12_DNPPE[,-c(grep("QC_",colnames(Exp_12_DNPPE)))]

# calculate concentrations

# pmol o.c.

Exp_12_PC.samp.pmol_oc = apply(Exp_12_PC.samp[,2:ncol(Exp_12_PC.samp)],c(1,2),splitpred,linfit_low.PC,linfit_hi.PC,1e10)
Exp_12_DNPPE.samp.pmol_oc = apply(Exp_12_DNPPE.samp[,2:ncol(Exp_12_DNPPE.samp)],c(1,2),splitpred,linfit_low.DNPPE,linfit_hi.DNPPE,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_12_PC.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_12_PC.samp.pmol_oc),getMetDat,meta.raw,c(2:7,14)))
Exp_12_PC.metdat$Date.time.sample.collected = strptime(as.character(Exp_12_PC.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_12_PC.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_12_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Exp_12_DNPPE.samp.pmol_oc # recovery factor
Exp_12_PC.samp.pmol.total = sweep(Exp_12_PC.samp.pmol_oc, 2, Exp_12_DNPPE.samp.RF, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_12_PC.samp.pmol.mL = sweep(Exp_12_PC.samp.pmol.total, 2, Exp_12_PC.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# Experiment 13

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_13_pos.RData")
Exp_13_pos.raw = getLOBpeaklist(Exp_13_pos) # generate peaklist

# extract only possible PC species from liposomes added to this experiment, and DNPPE

Exp_13_pos.subset = Exp_13_pos.raw[(Exp_13_pos.raw$species==c("PC") & (
  (Exp_13_pos.raw$FA_total_no_C==32 &
     Exp_13_pos.raw$FA_total_no_DB==0 &
     Exp_13_pos.raw$degree_oxidation==0) |
    (Exp_13_pos.raw$FA_total_no_C==36 &
       Exp_13_pos.raw$FA_total_no_DB<=2 &
       Exp_13_pos.raw$degree_oxidation>=0) |
    (Exp_13_pos.raw$FA_total_no_C==44 &
       Exp_13_pos.raw$FA_total_no_DB==0 &
       Exp_13_pos.raw$degree_oxidation==0) |
    (Exp_13_pos.raw$FA_total_no_C==44 &
       Exp_13_pos.raw$FA_total_no_DB %in% c(10,11,12) &
       Exp_13_pos.raw$degree_oxidation>=0))) |
    Exp_13_pos.raw$species==c("DNPPE")
  ,]

Exp_13_PC = Exp_13_pos.subset[grep("PC",Exp_13_pos.subset$compound_name),c(2,18:ncol(Exp_13_pos.subset))]

# append some more sample metadata
Exp_13_PC = cbind(Exp_13_PC,Exp_13_pos.subset[grep("PC",Exp_13_pos.subset$compound_name),c(1:17)])

Exp_13_DNPPE = Exp_13_pos.subset[grep("DNPPE",Exp_13_pos.subset$compound_name),c(2,18:54)]

# extract QCs

Exp_13_PC.QC = Exp_13_PC[,c(1,grep("QC_",colnames(Exp_13_PC)))]
Exp_13_PC.QC = Exp_13_PC.QC[Exp_13_PC.QC$compound_name=="PC 32:0",] # subset to just PC 16:0/16:0
Exp_13_DNPPE.QC = Exp_13_DNPPE[,c(1,grep("QC_",colnames(Exp_13_DNPPE)))]

# eliminate QCs from list of data

Exp_13_PC.samp = Exp_13_PC[,-c(grep("QC_",colnames(Exp_13_PC)))]
Exp_13_DNPPE.samp = Exp_13_DNPPE[,-c(grep("QC_",colnames(Exp_13_DNPPE)))]

# calculate concentrations

# pmol o.c.

Exp_13_PC.samp.pmol_oc = apply(Exp_13_PC.samp[,2:34],c(1,2),splitpred,linfit_low.PC,linfit_hi.PC,1e10)
Exp_13_DNPPE.samp.pmol_oc = apply(Exp_13_DNPPE.samp[,2:ncol(Exp_13_DNPPE.samp)],c(1,2),splitpred,linfit_low.DNPPE,linfit_hi.DNPPE,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_13_PC.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_13_PC.samp.pmol_oc),getMetDat,meta.raw,c(2:7,14)))
Exp_13_PC.metdat$Date.time.sample.collected = strptime(as.character(Exp_13_PC.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale pmol o.c. to initial pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_13_PC.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_13_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Exp_13_DNPPE.samp.pmol_oc # recovery factor
Exp_13_PC.samp.pmol.total = sweep(Exp_13_PC.samp.pmol_oc, 2, Exp_13_DNPPE.samp.RF, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_13_PC.samp.pmol.mL = sweep(Exp_13_PC.samp.pmol.total, 2, Exp_13_PC.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

### adjust experimental data for variance across time using QCs ###

# create aggregate objects containing QC data from the various Orbitrap sample runs (and the standard curve run)

all.QCs.PC = c(Exp_01_PC.QC[2:length(Exp_01_PC.QC)],
               Exp_02a_PC.QC[2:length(Exp_02a_PC.QC)],
               Exp_03a_PC.QC[2:length(Exp_03a_PC.QC)],
               Exp_12_PC.QC[2:length(Exp_12_PC.QC)],
               Exp_13_PC.QC[2:length(Exp_13_PC.QC)],
               PCstds_QC)
all.QCs.PC = unlist(all.QCs.PC[!duplicated(names(all.QCs.PC))]) # eliminate duplicates
all.QCs.PC = all.QCs.PC[sort(names(all.QCs.PC),index.return=TRUE)$ix] # reorder by Orbi or QE sample ID (corresponds to time)

all.QCs.DNPPE = c(Exp_01_DNPPE.QC[2:length(Exp_01_DNPPE.QC)],
                  Exp_02a_DNPPE.QC[2:length(Exp_02a_DNPPE.QC)],
                  Exp_03a_DNPPE.QC[2:length(Exp_03a_DNPPE.QC)],
                  Exp_12_DNPPE.QC[2:length(Exp_12_DNPPE.QC)],
                  Exp_13_DNPPE.QC[2:length(Exp_13_DNPPE.QC)],
                  DNPPEstds_QC
                  )
all.QCs.DNPPE = unlist(all.QCs.DNPPE[!duplicated(names(all.QCs.DNPPE))]) # eliminate duplicates
all.QCs.DNPPE = all.QCs.DNPPE[sort(names(all.QCs.DNPPE),index.return=TRUE)$ix] # reorder by Orbi or QE sample ID (corresponds to time)

# some plots

plot(all.QCs.DNPPE)
plot(all.QCs.PC)
plot(all.QCs.DNPPE,all.QCs.PC)

# calculate some QC correction factors assuming the two QCs run at the time of standards represent "1"

QC.PC.cf = all.QCs.PC/mean(all.QCs.PC[c("QE001276","QE001264")])
QC.DNPPE.cf = all.QCs.DNPPE/mean(all.QCs.DNPPE[c("QE001276","QE001264")])
QC.cf = apply(cbind(QC.PC.cf,QC.DNPPE.cf),1,mean) # a combined normalization factor

# apply the correction factors to data appropriately; while we're at it, update row names to reflect compound names

# Experiment 1

Exp_01_PC.samp.pmol.mL.norm = Exp_01_PC.samp.pmol.mL
rownames(Exp_01_PC.samp.pmol.mL.norm) = Exp_01_PC.samp[,1]

for (i in 1:ncol(Exp_01_PC.samp.pmol.mL.norm)) {
  
  # extract cf's for nearest QCs to this sample, based on what's in the metadata table
  these.QC.IDs = unlist(strsplit(as.character(Exp_01_PC.metdat$Relevant.QC[i]),";"))
  these.cfs = QC.PC.cf[sub("^.*QC_","\\1",names(QC.PC.cf)) %in% these.QC.IDs]
  cf = mean(these.cfs)
  
  # normalize appropriately
  Exp_01_PC.samp.pmol.mL.norm[,i] = Exp_01_PC.samp.pmol.mL[,i]*cf
  
}

# Experiment 2a

Exp_02a_PC.samp.pmol.mL.norm = Exp_02a_PC.samp.pmol.mL
rownames(Exp_02a_PC.samp.pmol.mL.norm) = Exp_02a_PC.samp[,1]

for (i in 1:ncol(Exp_02a_PC.samp.pmol.mL.norm)) {
  
  # extract cf's for nearest QCs to this sample, based on what's in the metadata table
  these.QC.IDs = unlist(strsplit(as.character(Exp_02a_PC.metdat$Relevant.QC[i]),";"))
  these.cfs = QC.PC.cf[sub("^.*QC_","\\1",names(QC.PC.cf)) %in% these.QC.IDs]
  cf = mean(these.cfs)
  
  # normalize appropriately
  Exp_02a_PC.samp.pmol.mL.norm[,i] = Exp_02a_PC.samp.pmol.mL[,i]*cf
  
}

# Experiment 3a, initial run

Exp_03a_PC.samp.pmol.mL.norm = Exp_03a_PC.samp.pmol.mL
rownames(Exp_03a_PC.samp.pmol.mL.norm) = Exp_03a_PC.samp[,1]

for (i in 1:ncol(Exp_03a_PC.samp.pmol.mL.norm)) {
  
  # extract cf's for nearest QCs to this sample, based on what's in the metadata table
  these.QC.IDs = unlist(strsplit(as.character(Exp_03a_PC.metdat$Relevant.QC[i]),";"))
  these.cfs = QC.PC.cf[sub("^.*QC_","\\1",names(QC.PC.cf)) %in% these.QC.IDs]
  cf = mean(these.cfs)
  
  # normalize appropriately
  Exp_03a_PC.samp.pmol.mL.norm[,i] = Exp_03a_PC.samp.pmol.mL[,i]*cf
  
}

# Experiment 3a, rerun

# don't need to make this adjustment since standard were run right with the samples

Exp_03a_PC.samp.pmol.mL.norm.rerun = Exp_03a_PC.samp.pmol.mL.rerun
rownames(Exp_03a_PC.samp.pmol.mL.norm.rerun) = Exp_03a_PC.rerun[,1]

# Experiment 12

Exp_12_PC.samp.pmol.mL.norm = Exp_12_PC.samp.pmol.mL
rownames(Exp_12_PC.samp.pmol.mL.norm) = Exp_12_PC.samp[,1]

for (i in 1:ncol(Exp_12_PC.samp.pmol.mL.norm)) {
  
  # extract cf's for nearest QCs to this sample, based on what's in the metadata table
  these.QC.IDs = unlist(strsplit(as.character(Exp_12_PC.metdat$Relevant.QC[i]),";"))
  these.cfs = QC.PC.cf[sub("^.*QC_","\\1",names(QC.PC.cf)) %in% these.QC.IDs]
  cf = mean(these.cfs)
  
  # normalize appropriately
  Exp_12_PC.samp.pmol.mL.norm[,i] = Exp_12_PC.samp.pmol.mL[,i]*cf
  
}

# Experiment 13

Exp_13_PC.samp.pmol.mL.norm = Exp_13_PC.samp.pmol.mL
rownames(Exp_13_PC.samp.pmol.mL.norm) = Exp_13_PC.samp[,1]

for (i in 1:ncol(Exp_13_PC.samp.pmol.mL.norm)) {
  
  # extract cf's for nearest QCs to this sample, based on what's in the metadata table
  these.QC.IDs = unlist(strsplit(as.character(Exp_13_PC.metdat$Relevant.QC[i]),";"))
  these.cfs = QC.PC.cf[sub("^.*QC_","\\1",names(QC.PC.cf)) %in% these.QC.IDs]
  cf = mean(these.cfs)
  
  # normalize appropriately
  Exp_13_PC.samp.pmol.mL.norm[,i] = Exp_13_PC.samp.pmol.mL[,i]*cf
  
}

### now, can actually do something with the data ###

# Experiment 1

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_01_PC.samp.pmol.mL.norm)) {
  
  for (j in 1:length(unique(Exp_01_PC.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_01_PC.metdat[Exp_01_PC.metdat$Treatment.ID==unique(Exp_01_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_01_PC.samp.pmol.mL.norm[i,Exp_01_PC.metdat$Treatment.ID==unique(Exp_01_PC.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_01_PC.samp.pmol.mL.norm[i,]),max(Exp_01_PC.samp.pmol.mL.norm[i,]))
           )
      }  else {
        points(Exp_01_PC.metdat[Exp_01_PC.metdat$Treatment.ID==unique(Exp_01_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
               Exp_01_PC.samp.pmol.mL.norm[i,Exp_01_PC.metdat$Treatment.ID==unique(Exp_01_PC.metdat$Treatment.ID)[j]],pch=j,col=j)
      }
    
  }
  
  title(main = paste0("Exp 1: ",rownames(Exp_01_PC.samp.pmol.mL.norm)[i]))
  legend(x="topright",
         legend=unique(Exp_01_PC.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Experiment 02a

# need to fix some data for samples that were switched (see PAL1314 lab notebook, p. 22)

Exp_02a_PC.metdat[colnames(Exp_02a_PC.samp.pmol.mL.norm) %in% c("EPA_no_HB_20131009_1900_Orbi_2479",
                                                                "EPA_no_HB_20131009_1900_QE000907",
                                                                "EPA_no_HB_20131009_1900_QE000908"),
                  c("Treatment.ID")] = "Dark_control_no_HB"

Exp_02a_PC.metdat[colnames(Exp_02a_PC.samp.pmol.mL.norm) %in% c("Dark_control_no_HB_20131009_1900_Orbi_2474",
                                                                "Dark_control_no_HB_20131009_1900_Orbi_2475",
                                                                "Dark_control_no_HB_20131009_1900_Orbi_2476"),
                  c("Treatment.ID")] = "EPA_no_HB"

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_02a_PC.samp.pmol.mL.norm)) {
  
  for (j in 1:length(unique(Exp_02a_PC.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_02a_PC.metdat[Exp_02a_PC.metdat$Treatment.ID==unique(Exp_02a_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_02a_PC.samp.pmol.mL.norm[i,Exp_02a_PC.metdat$Treatment.ID==unique(Exp_02a_PC.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_02a_PC.samp.pmol.mL.norm[i,]),max(Exp_02a_PC.samp.pmol.mL.norm[i,]))
           )
    }  else {
      points(Exp_02a_PC.metdat[Exp_02a_PC.metdat$Treatment.ID==unique(Exp_02a_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_02a_PC.samp.pmol.mL.norm[i,Exp_02a_PC.metdat$Treatment.ID==unique(Exp_02a_PC.metdat$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 2a: ",rownames(Exp_02a_PC.samp.pmol.mL.norm)[i]))
  legend(x="topright",
         legend=unique(Exp_02a_PC.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Experiment 3a, initial run

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_03a_PC.samp.pmol.mL.norm)) {
  
  for (j in 1:length(unique(Exp_03a_PC.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_03a_PC.metdat[Exp_03a_PC.metdat$Treatment.ID==unique(Exp_03a_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_03a_PC.samp.pmol.mL.norm[i,Exp_03a_PC.metdat$Treatment.ID==unique(Exp_03a_PC.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_03a_PC.samp.pmol.mL.norm[i,]),max(Exp_03a_PC.samp.pmol.mL.norm[i,]))
      )
    }  else {
      points(Exp_03a_PC.metdat[Exp_03a_PC.metdat$Treatment.ID==unique(Exp_03a_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_03a_PC.samp.pmol.mL.norm[i,Exp_03a_PC.metdat$Treatment.ID==unique(Exp_03a_PC.metdat$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 3a: ",rownames(Exp_03a_PC.samp.pmol.mL.norm)[i]))
  legend(x="topright",
         legend=unique(Exp_03a_PC.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Experiment 3a, rerun

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_03a_PC.samp.pmol.mL.norm.rerun)) {
  
  for (j in 1:length(unique(Exp_03a_PC.metdat.rerun$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_03a_PC.metdat.rerun[Exp_03a_PC.metdat.rerun$Treatment.ID==unique(Exp_03a_PC.metdat.rerun$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_03a_PC.samp.pmol.mL.norm.rerun[i,Exp_03a_PC.metdat.rerun$Treatment.ID==unique(Exp_03a_PC.metdat.rerun$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_03a_PC.samp.pmol.mL.norm.rerun[i,]),max(Exp_03a_PC.samp.pmol.mL.norm.rerun[i,]))
      )
    }  else {
      points(Exp_03a_PC.metdat.rerun[Exp_03a_PC.metdat.rerun$Treatment.ID==unique(Exp_03a_PC.metdat.rerun$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_03a_PC.samp.pmol.mL.norm.rerun[i,Exp_03a_PC.metdat.rerun$Treatment.ID==unique(Exp_03a_PC.metdat.rerun$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 3a: ",rownames(Exp_03a_PC.samp.pmol.mL.norm.rerun)[i]))
  legend(x="topright",
         legend=unique(Exp_03a_PC.metdat.rerun$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Experiment 12

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_12_PC.samp.pmol.mL.norm)) {
  
  for (j in 1:length(unique(Exp_12_PC.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_12_PC.metdat[Exp_12_PC.metdat$Treatment.ID==unique(Exp_12_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_12_PC.samp.pmol.mL.norm[i,Exp_12_PC.metdat$Treatment.ID==unique(Exp_12_PC.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_12_PC.samp.pmol.mL.norm[i,], na.rm = TRUE),max(Exp_12_PC.samp.pmol.mL.norm[i,], na.rm = TRUE))
      )
    }  else {
      points(Exp_12_PC.metdat[Exp_12_PC.metdat$Treatment.ID==unique(Exp_12_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_12_PC.samp.pmol.mL.norm[i,Exp_12_PC.metdat$Treatment.ID==unique(Exp_12_PC.metdat$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 12: ",rownames(Exp_12_PC.samp.pmol.mL.norm)[i]))
  legend(x="topright",
         legend=unique(Exp_12_PC.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Experiment 13

# exploratory plots, all data by compound ID

for (i in 1:nrow(Exp_13_PC.samp.pmol.mL.norm)) {
  
  for (j in 1:length(unique(Exp_13_PC.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_13_PC.metdat[Exp_13_PC.metdat$Treatment.ID==unique(Exp_13_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_13_PC.samp.pmol.mL.norm[i,Exp_13_PC.metdat$Treatment.ID==unique(Exp_13_PC.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="Concentration (pmol/mL)",
           xlab="Time",
           ylim=c(min(Exp_13_PC.samp.pmol.mL.norm[i,]),max(Exp_13_PC.samp.pmol.mL.norm[i,]))
      )
    }  else {
      points(Exp_13_PC.metdat[Exp_13_PC.metdat$Treatment.ID==unique(Exp_13_PC.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_13_PC.samp.pmol.mL.norm[i,Exp_13_PC.metdat$Treatment.ID==unique(Exp_13_PC.metdat$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 13: ",rownames(Exp_13_PC.samp.pmol.mL.norm)[i]))
  legend(x="topright",
         legend=unique(Exp_13_PC.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# Perform multiple-way ANOVAs between t=0 and t=final timepoints for each experiment,
# to determine whether we observed any statistically significant removal of added moieties

# Experiment 1

# print name of experiment
print("Experiment 1")

# subset to only first and last timepoint

Exp_01_PC.fl = Exp_01_PC.samp.pmol.mL.norm[,(Exp_01_PC.metdat$Date.time.sample.collected %in% unique(Exp_01_PC.metdat$Date.time.sample.collected)[c(1,3,6)])]
Exp_01_PC.metdat.fl = Exp_01_PC.metdat[(Exp_01_PC.metdat$Date.time.sample.collected %in% unique(Exp_01_PC.metdat$Date.time.sample.collected)[c(1,3,6)]),]
Exp_01_PC.metdat.fl$ttp.ID = paste0(Exp_01_PC.metdat.fl$Treatment.ID,"_",Exp_01_PC.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_01_PC.fl)) { # subset by moiety
  
  print(rownames(Exp_01_PC.fl)[i]) # print name of this moiety
  
  Exp_01_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_01_PC.fl[i,]),Exp_01_PC.metdat.fl$ttp.ID))
  Exp_01_PC.fl.subs$V1 = as.numeric(as.character(Exp_01_PC.fl.subs$V1))
  colnames(Exp_01_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_01_PC.fl.subs)
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov)
  print(tukey)
}

# Experiment 2a

# print name of experiment
print("Experiment 2a")

# subset to only first and last timepoint

Exp_02a_PC.fl = Exp_02a_PC.samp.pmol.mL.norm[,(Exp_02a_PC.metdat$Date.time.sample.collected %in% unique(Exp_02a_PC.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_02a_PC.metdat.fl = Exp_02a_PC.metdat[(Exp_02a_PC.metdat$Date.time.sample.collected %in% unique(Exp_02a_PC.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_02a_PC.metdat.fl$ttp.ID = paste0(Exp_02a_PC.metdat.fl$Treatment.ID,"_",Exp_02a_PC.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_02a_PC.fl)) { # subset by moiety
  
  print(rownames(Exp_02a_PC.fl)[i]) # print name of this moiety
  
  Exp_02a_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_02a_PC.fl[i,]),Exp_02a_PC.metdat.fl$ttp.ID))
  Exp_02a_PC.fl.subs$V1 = as.numeric(as.character(Exp_02a_PC.fl.subs$V1))
  colnames(Exp_02a_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_02a_PC.fl.subs)
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov)
  print(tukey)
  
}

# Experiment 3a, initial run

# print name of experiment
print("Experiment 3a, initial run")

# subset to only first and last timepoint

Exp_03a_PC.fl = Exp_03a_PC.samp.pmol.mL.norm[,(Exp_03a_PC.metdat$Date.time.sample.collected %in% unique(Exp_03a_PC.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_03a_PC.metdat.fl = Exp_03a_PC.metdat[(Exp_03a_PC.metdat$Date.time.sample.collected %in% unique(Exp_03a_PC.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_03a_PC.metdat.fl$ttp.ID = paste0(Exp_03a_PC.metdat.fl$Treatment.ID,"_",Exp_03a_PC.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_03a_PC.fl)) { # subset by moiety
  
  print(rownames(Exp_03a_PC.fl)[i]) # print name of this moiety
  
  Exp_03a_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_03a_PC.fl[i,]),Exp_03a_PC.metdat.fl$ttp.ID))
  Exp_03a_PC.fl.subs$V1 = as.numeric(as.character(Exp_03a_PC.fl.subs$V1))
  colnames(Exp_03a_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
 Exp_03a_PC.fl.subs= Exp_03a_PC.fl.subs[-2,]
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_03a_PC.fl.subs)
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov)
  print(tukey)
  
}

# Experiment 3a, rerun

# print name of experiment
print("Experiment 3a, rerun")

# subset to only first and last timepoint

Exp_03a_PC.fl.rerun = Exp_03a_PC.samp.pmol.mL.norm.rerun[,(Exp_03a_PC.metdat.rerun$Date.time.sample.collected %in% unique(Exp_03a_PC.metdat.rerun$Date.time.sample.collected)[c(1,3)])]
Exp_03a_PC.metdat.fl.rerun = Exp_03a_PC.metdat.rerun[(Exp_03a_PC.metdat.rerun$Date.time.sample.collected %in% unique(Exp_03a_PC.metdat.rerun$Date.time.sample.collected)[c(1,3)]),]
Exp_03a_PC.metdat.fl.rerun$ttp.ID = paste0(Exp_03a_PC.metdat.fl.rerun$Treatment.ID,"_",Exp_03a_PC.metdat.fl.rerun$Date.time.sample.collected) # create a single treatment-timepoint ID

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_03a_PC.fl.rerun)) { # subset by moiety
  
  print(rownames(Exp_03a_PC.fl.rerun)[i]) # print name of this moiety
  
  Exp_03a_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_03a_PC.fl.rerun[i,]),Exp_03a_PC.metdat.fl.rerun$ttp.ID))
  Exp_03a_PC.fl.subs$V1 = as.numeric(as.character(Exp_03a_PC.fl.subs$V1))
  colnames(Exp_03a_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  Exp_03a_PC.fl.subs= Exp_03a_PC.fl.subs[-2,]
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_03a_PC.fl.subs)
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov)
  print(tukey)
  
}

# Experiment 12

# print name of experiment
print("Experiment 12")

# subset to only first and last timepoint

Exp_12_PC.fl = Exp_12_PC.samp.pmol.mL.norm[,(Exp_12_PC.metdat$Date.time.sample.collected %in% unique(Exp_12_PC.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_12_PC.metdat.fl = Exp_12_PC.metdat[(Exp_12_PC.metdat$Date.time.sample.collected %in% unique(Exp_12_PC.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_12_PC.metdat.fl$ttp.ID = paste0(Exp_12_PC.metdat.fl$Treatment.ID,"_",Exp_12_PC.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_12_PC.fl)) { # subset by moiety
  
  print(rownames(Exp_12_PC.fl)[i]) # print name of this moiety
  
  Exp_12_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_12_PC.fl[i,]),Exp_12_PC.metdat.fl$ttp.ID))
  Exp_12_PC.fl.subs$V1 = as.numeric(as.character(Exp_12_PC.fl.subs$V1))
  colnames(Exp_12_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_12_PC.fl.subs)
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov)
  print(tukey)
  
}

# Experiment 13

# print name of experiment
print("Experiment 13")

# subset to only first and last timepoint

Exp_13_PC.fl = Exp_13_PC.samp.pmol.mL.norm[,(Exp_13_PC.metdat$Date.time.sample.collected %in% unique(Exp_13_PC.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_13_PC.metdat.fl = Exp_13_PC.metdat[(Exp_13_PC.metdat$Date.time.sample.collected %in% unique(Exp_13_PC.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_13_PC.metdat.fl$ttp.ID = paste0(Exp_13_PC.metdat.fl$Treatment.ID,"_",Exp_13_PC.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# combine the +HB and -HB controls, since have no apparent difference
Exp_13_PC.metdat.fl$ttp.ID[Exp_13_PC.metdat.fl$ttp.ID %in% c("Dark_control_no_HB_2013-12-14 09:30:00",
                                                 "Dark_control_plus_HB_2013-12-14 09:30:00")] = c("Dark_control_2013-12-14 09:30:00")

# eliminate a bad sample
Exp_13_PC.fl = Exp_13_PC.fl[,-c(6)]
Exp_13_PC.metdat.fl = Exp_13_PC.metdat.fl[-c(6),]

# create subsets for a given moiety and perform ANOVA, then Tukey HSD 

for (i in 1:nrow(Exp_13_PC.fl)) { # subset by moiety
  
  print(rownames(Exp_13_PC.fl)[i]) # print name of this moiety
  
  Exp_13_PC.fl.subs = as.data.frame(cbind(as.numeric(Exp_13_PC.fl[i,]),Exp_13_PC.metdat.fl$ttp.ID))
  
  Exp_13_PC.fl.subs$V1 = as.numeric(as.character(Exp_13_PC.fl.subs$V1))
  colnames(Exp_13_PC.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  PC.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_13_PC.fl.subs)
  
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov, conf.level = 0.95)
  print(tukey)
  
  # calculate & display mean differences ± SD for selected treatment pairs
  
  xbar_init = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  sd_init = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  se_init = sd_init/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"]))

    xbar_final = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("+ UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("- UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("- UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"]))

    delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB, + HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("+ UVB, + HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_no_HB_2013-12-14 17:50:00"]))

    delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("dark control, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("dark control, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_plus_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_plus_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_PC.fl.subs$Conc_pmol_mL[Exp_13_PC.fl.subs$Treatment=="Dark_control_plus_HB_2013-12-14 17:50:00"]))

    delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("dark control, + HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n\n")
  cat("dark control, + HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n\n")

  }

# add'l plots for Exp't 13

Exp_13_PC.metdat$ttp.ID = paste0(Exp_13_PC.metdat$Treatment.ID,"_",Exp_13_PC.metdat$Date.time.sample.collected) # create a single treatment-timepoint ID

# combine the +HB and -HB controls, since have no apparent difference
Exp_13_PC.metdat$ttp.ID[Exp_13_PC.metdat$ttp.ID %in% c("Dark_control_no_HB_2013-12-14 09:30:00",
                                                             "Dark_control_plus_HB_2013-12-14 09:30:00")] = c("Dark_control_2013-12-14 09:30:00")

# eliminate the bad sample
Exp_13_PC.metdat = Exp_13_PC.metdat[-c(9),]
Exp_13_PC.samp.pmol.mL.norm = Exp_13_PC.samp.pmol.mL.norm[,-c(9)]

# calculate mean & SD for each set of replicates

# preallocate array

Exp_13_PC.unique.ttps = unique(Exp_13_PC.metdat$ttp.ID)

Exp_13_PC.samp.pmol.mL.norm.mean = matrix(data = NA, 
                                          nrow = nrow(Exp_13_PC.samp.pmol.mL.norm),
                                          ncol = 3*length(Exp_13_PC.unique.ttps)
)

rownames(Exp_13_PC.samp.pmol.mL.norm.mean) = rownames(Exp_13_PC.samp.pmol.mL.norm)
colnames(Exp_13_PC.samp.pmol.mL.norm.mean) = rep("",3*length(Exp_13_PC.unique.ttps))

# calculate stats, populate array

for (i in 1:length(Exp_13_PC.unique.ttps)) {
  
  current.data = Exp_13_PC.samp.pmol.mL.norm[,Exp_13_PC.metdat$ttp.ID==Exp_13_PC.unique.ttps[i]]
  
  mean.current = apply(current.data,1,mean)
  sd.current = apply(current.data,1,sd)
  se.current = sd.current/sqrt(ncol(current.data))

  # insert into our array
  
  Exp_13_PC.samp.pmol.mL.norm.mean[,3*i-2] = mean.current
  Exp_13_PC.samp.pmol.mL.norm.mean[,3*i-1] = sd.current
  Exp_13_PC.samp.pmol.mL.norm.mean[,3*i] = se.current
  
  # update column labels
  
  colnames(Exp_13_PC.samp.pmol.mL.norm.mean)[(3*i-2)] =
    paste0(Exp_13_PC.unique.ttps[i],".mean")
  colnames(Exp_13_PC.samp.pmol.mL.norm.mean)[(3*i-1)] =
    paste0(Exp_13_PC.unique.ttps[i],".sd")
  colnames(Exp_13_PC.samp.pmol.mL.norm.mean)[(3*i)] =
    paste0(Exp_13_PC.unique.ttps[i],".se")
  
  
}


# bar plots

# PC 22:6

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Exp13_PC22-6_barplot.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey","darkgrey")
treatments <- c("Dark","-UVB, -het. bact.","+UVB, -het. bact.","+UVB, +het. bact.")

bardata.mean <- matrix(NA,4,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.mean",
                                                                                                                                                        "Dark_control_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                        "Dark_control_no_HB_2013-12-14 17:50:00.mean")]
bardata.mean[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                            "Quartz_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                            "Quartz_no_HB_2013-12-14 17:50:00.mean")]
bardata.mean[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                            "Quartz_plus_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                            "Quartz_plus_HB_2013-12-14 17:50:00.mean")]
bardata.mean[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                            "EPA_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                            "EPA_no_HB_2013-12-14 17:50:00.mean")]

bardata.se <- matrix(NA,4,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.se",
                                                                                                                                                            "Dark_control_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                            "Dark_control_no_HB_2013-12-14 17:50:00.se")]
bardata.se[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.se",
                                                                                                                                                               "Quartz_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                               "Quartz_no_HB_2013-12-14 17:50:00.se")]
bardata.se[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.se",
                                                                                                                                                               "Quartz_plus_HB_2013-12-14 13:40:00.se",
                                                                                                                                                               "Quartz_plus_HB_2013-12-14 17:50:00.se")]
bardata.se[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.se",
                                                                                                                                                               "EPA_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                               "EPA_no_HB_2013-12-14 17:50:00.se")]

bardata.sd <- matrix(NA,4,3)
bardata.sd[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.sd",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 17:50:00.sd")]
bardata.sd[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 17:50:00.sd")]
bardata.sd[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 17:50:00.sd")]
bardata.sd[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "EPA_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "EPA_no_HB_2013-12-14 17:50:00.sd")]


plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 22:6, 22:6 PC (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topright",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()

# PC 18:1

par(oma = c(1, 2, 4, 1))
barcolors <- c("black","darkgrey","lightgrey")
treatments <- c("Dark","-UVB","+UVB")

bardata.mean <- matrix(NA,3,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(1,4,7)]
bardata.mean[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(NA,19,22)]
bardata.mean[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(NA,25,28)]

bardata.se <- matrix(NA,3,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(3,6,9)]
bardata.se[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(NA,21,24)]
bardata.se[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:2",c(NA,27,30)]

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 18:1, 18:1 PC (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topright",treatments,fill=barcolors,xpd = TRUE, inset = c(0,-.45))

# PC 18:0

par(oma = c(1, 2, 4, 1))
barcolors <- c("black","darkgrey","lightgrey")
treatments <- c("Dark","-UVB","+UVB")

bardata.mean <- matrix(NA,3,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(1,4,7)]
bardata.mean[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(NA,19,22)]
bardata.mean[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(NA,25,28)]

bardata.se <- matrix(NA,3,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(3,6,9)]
bardata.se[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(NA,21,24)]
bardata.se[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 36:0",c(NA,27,30)]

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 18:0, 18:0 PC (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topright",treatments,fill=barcolors,xpd = TRUE, inset = c(0,-.45))

# PC 16:0

par(oma = c(1, 2, 4, 1))
barcolors <- c("black","darkgrey","lightgrey")
treatments <- c("Dark","-UVB","+UVB")

bardata.mean <- matrix(NA,3,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(1,4,7)]
bardata.mean[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(NA,19,22)]
bardata.mean[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(NA,25,28)]

bardata.se <- matrix(NA,3,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(3,6,9)]
bardata.se[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(NA,21,24)]
bardata.se[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 32:0",c(NA,27,30)]

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 16:0, 16:0 PC (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topright",treatments,fill=barcolors,xpd = TRUE, inset = c(0,-.45))

# PC 22:6 +2O

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Exp13_PC22-6+2O_barplot.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey","darkgrey")
treatments <- c("Dark","-UVB, -het. bact.","+UVB, -het. bact.","+UVB, +het. bact.")

bardata.mean <- matrix(NA,4,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.mean",
                                                                                                                                                            "Dark_control_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                            "Dark_control_no_HB_2013-12-14 17:50:00.mean")]
bardata.mean[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                               "Quartz_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                               "Quartz_no_HB_2013-12-14 17:50:00.mean")]
bardata.mean[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                               "Quartz_plus_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                               "Quartz_plus_HB_2013-12-14 17:50:00.mean")]
bardata.mean[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.mean",
                                                                                                                                                               "EPA_no_HB_2013-12-14 13:40:00.mean",
                                                                                                                                                               "EPA_no_HB_2013-12-14 17:50:00.mean")]

bardata.se <- matrix(NA,4,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.se",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 17:50:00.se")]
bardata.se[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.se",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 17:50:00.se")]
bardata.se[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.se",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 13:40:00.se",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 17:50:00.se")]
bardata.se[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.se",
                                                                                                                                                             "EPA_no_HB_2013-12-14 13:40:00.se",
                                                                                                                                                             "EPA_no_HB_2013-12-14 17:50:00.se")]

bardata.sd <- matrix(NA,4,3)
bardata.sd[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Dark_control_2013-12-14 09:30:00.sd",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                          "Dark_control_no_HB_2013-12-14 17:50:00.sd")]
bardata.sd[3,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_no_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "Quartz_no_HB_2013-12-14 17:50:00.sd")]
bardata.sd[4,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("Quartz_plus_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "Quartz_plus_HB_2013-12-14 17:50:00.sd")]
bardata.sd[2,2:3]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +2O",colnames(Exp_13_PC.samp.pmol.mL.norm.mean) %in% c("EPA_no_HB_2013-12-14 09:30:00.sd",
                                                                                                                                                             "EPA_no_HB_2013-12-14 13:40:00.sd",
                                                                                                                                                             "EPA_no_HB_2013-12-14 17:50:00.sd")]


plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 22:6, 22:6 PC +2O (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topleft",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()

# PC 22:6 +4O

par(oma = c(1, 2, 4, 1))
barcolors <- c("black","darkgrey","lightgrey")
treatments <- c("Dark","-UVB","+UVB")

bardata.mean <- matrix(NA,3,3)
bardata.mean[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(1,4,7)]
bardata.mean[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(NA,19,22)]
bardata.mean[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(NA,25,28)]

bardata.se <- matrix(NA,3,3)
bardata.se[1,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(3,6,9)]
bardata.se[2,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(NA,21,24)]
bardata.se[3,]<-Exp_13_PC.samp.pmol.mL.norm.mean[rownames(Exp_13_PC.samp.pmol.mL.norm.mean)=="PC 44:12 +4O",c(NA,27,30)]

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("Concentration 22:6, 22:6 PC +4O (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topright",treatments,fill=barcolors,xpd = TRUE, inset = c(0,-.45))

#### negative ion mode ###

# Experiment 13

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/Exp_13_neg_12ppm.RData")

# note that this data was screened in LOBSTAHS using a 12 ppm match tolerance
# lock mass was not working properly in - mode during data acquisiton and the
# DNPPE peak is ~ 11.5 ppm off

Exp_13_neg.raw = getLOBpeaklist(Exp_13_neg) # generate peaklist

# extract only FFA, DNPPE, LPC

Exp_13_neg.subset = Exp_13_neg.raw[(Exp_13_neg.raw$lipid_class==c("FFA") & (
  (Exp_13_neg.raw$FA_total_no_C %in% c(16) & Exp_13_neg.raw$FA_total_no_DB==0) |
    (Exp_13_neg.raw$FA_total_no_C %in% c(18) & Exp_13_neg.raw$FA_total_no_DB<=1) |
    (Exp_13_neg.raw$FA_total_no_C %in% c(22) & Exp_13_neg.raw$FA_total_no_DB %in% c(5,6))
  )) |
    Exp_13_neg.raw$species==c("DNPPE") |
    (Exp_13_neg.raw$lipid_class==c("IP_MAG") & (
      (Exp_13_neg.raw$FA_total_no_C %in% c(22) & Exp_13_neg.raw$FA_total_no_DB==6)))
  ,]

Exp_13_FFA.neg = Exp_13_neg.subset[grep("FFA",Exp_13_neg.subset$compound_name),c(2,18:54)]
Exp_13_LPC.neg = Exp_13_neg.subset[grep("LPC",Exp_13_neg.subset$compound_name),c(2,18:54)]
Exp_13_DNPPE.neg = Exp_13_neg.subset[grep("DNPPE",Exp_13_neg.subset$compound_name),c(2,18:54)]

# extract QCs

Exp_13_DNPPE.neg.QC = Exp_13_DNPPE.neg[,c(1,grep("QC_",colnames(Exp_13_DNPPE.neg)))]

# eliminate QCs from list of data

Exp_13_FFA.neg.samp = Exp_13_FFA.neg[,-c(grep("QC_",colnames(Exp_13_FFA.neg)))]
Exp_13_LPC.neg.samp = Exp_13_LPC.neg[,-c(grep("QC_",colnames(Exp_13_LPC.neg)))]
Exp_13_DNPPE.neg.samp = Exp_13_DNPPE.neg[,-c(grep("QC_",colnames(Exp_13_DNPPE.neg)))]

# eliminate bad sample from negative mode data as well
Exp_13_LPC.neg.samp = Exp_13_LPC.neg.samp[,-c(10)]
Exp_13_FFA.neg.samp = Exp_13_FFA.neg.samp[,-c(10)]
Exp_13_DNPPE.neg.samp = Exp_13_DNPPE.neg.samp[,-c(10)]

# pmol o.c.

Exp_13_FFA.neg.samp.pmol_oc = apply(Exp_13_FFA.neg.samp[,2:ncol(Exp_13_FFA.neg.samp)],c(1,2),splitpred,linfit.DHA.neg,linfit.DHA.neg,1e10)
Exp_13_LPC.neg.samp.pmol_oc = apply(Exp_13_LPC.neg.samp[,2:ncol(Exp_13_LPC.neg.samp)],c(1,2),splitpred,linfit_low.PC.neg,linfit_hi.PC.neg,1e10)
Exp_13_DNPPE.samp.neg.pmol_oc = apply(Exp_13_DNPPE.neg.samp[,2:ncol(Exp_13_DNPPE.neg.samp)],c(1,2),splitpred,linfit_low.DNPPE.neg,linfit_hi.DNPPE.neg,1.25e9)

# retrieve necessary metadata
# do.call(rbind.data.frame,...) syntax necessary to prevent coercion of factors to integers (as would happen with unlist())

Exp_13_FFA.metdat = do.call(rbind.data.frame,lapply(colnames(Exp_13_FFA.neg.samp),getMetDat,meta.raw,c(2:7,14)))
Exp_13_FFA.metdat$Date.time.sample.collected = strptime(as.character(Exp_13_FFA.metdat$Date.time.sample.collected),"%m/%d/%y %H:%M")

# scale peak areas using DNPPE (recovery standard added at time of extraction)

DNPPE_pmol_added_per_samp = DNPPE_mg_mL*(1/DNPPE_MW)*(10^9)*(1/10^3)*Exp_13_FFA.metdat$Vol.DNP.PE..uL. # quantity of DNPPE (pmol) added per sample these experiments (DNPPE was added to vials containing 40 or 45 mL sample just prior to liquid/liquid extraction); should be 30 uL for almost all samples, except a few as noted
Exp_13_DNPPE.samp.RF.neg = DNPPE_pmol_added_per_samp/Exp_13_DNPPE.samp.neg.pmol_oc # recovery factor

Exp_13_LPC.samp.pmol.total.neg = sweep(Exp_13_LPC.neg.samp.pmol_oc, 2, Exp_13_DNPPE.samp.RF.neg, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_13_LPC.samp.pmol.mL.neg = sweep(Exp_13_LPC.samp.pmol.total.neg, 2, Exp_13_FFA.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

Exp_13_FFA.neg.samp.total.neg = sweep(Exp_13_FFA.neg.samp.pmol_oc, 2, Exp_13_DNPPE.samp.RF.neg, "*") # apply RF to samples, calculate total # pmol each species in given sample
Exp_13_FFA.neg.samp.pmol.mL = sweep(Exp_13_FFA.neg.samp.total.neg, 2, Exp_13_FFA.metdat$Vol.sample.extracted.or.filtered..mL., "/")  # calculate pmol/mL, using correct volumes

# generate matrices w/summary stats

# for the FFAs

Exp_13_FFA.neg.samp.pmol.mL.mean = matrix(data = NA, 
                                          nrow = nrow(Exp_13_FFA.neg.samp.pmol.mL),
                                          ncol = 3*length(Exp_13_PC.unique.ttps)
)

rownames(Exp_13_FFA.neg.samp.pmol.mL.mean) = Exp_13_FFA.neg.samp$compound_name
colnames(Exp_13_FFA.neg.samp.pmol.mL.mean) = rep("",3*length(Exp_13_PC.unique.ttps))

# calculate stats, populate array

for (i in 1:length(Exp_13_PC.unique.ttps)) {
  
  current.data = Exp_13_FFA.neg.samp.pmol.mL[,Exp_13_PC.metdat$ttp.ID==Exp_13_PC.unique.ttps[i]]
  
  mean.current = apply(current.data,1,mean)
  sd.current = apply(current.data,1,sd)
  se.current = sd.current/sqrt(ncol(current.data))
  
  # insert into our array
  
  Exp_13_FFA.neg.samp.pmol.mL.mean[,3*i-2] = mean.current
  Exp_13_FFA.neg.samp.pmol.mL.mean[,3*i-1] = sd.current
  Exp_13_FFA.neg.samp.pmol.mL.mean[,3*i] = se.current
  
  # update column labels
  
  colnames(Exp_13_FFA.neg.samp.pmol.mL.mean)[(3*i-2)] =
    paste0(Exp_13_PC.unique.ttps[i],".mean")
  colnames(Exp_13_FFA.neg.samp.pmol.mL.mean)[(3*i-1)] =
    paste0(Exp_13_PC.unique.ttps[i],".sd")
  colnames(Exp_13_FFA.neg.samp.pmol.mL.mean)[(3*i)] =
    paste0(Exp_13_PC.unique.ttps[i],".se")
  
  
}

# for the LPCs

Exp_13_LPC.samp.pmol.mL.neg.mean = matrix(data = NA, 
                                       nrow = nrow(Exp_13_LPC.samp.pmol.mL.neg),
                                       ncol = 3*length(Exp_13_PC.unique.ttps)
)

rownames(Exp_13_LPC.samp.pmol.mL.neg.mean) = Exp_13_LPC.neg$compound_name
colnames(Exp_13_LPC.samp.pmol.mL.neg.mean) = rep("",3*length(Exp_13_PC.unique.ttps))

# calculate stats, populate array

for (i in 1:length(Exp_13_PC.unique.ttps)) {
  
  current.data = Exp_13_LPC.samp.pmol.mL.neg[,Exp_13_PC.metdat$ttp.ID==Exp_13_PC.unique.ttps[i]]
  
  mean.current = apply(current.data,1,mean)
  sd.current = apply(current.data,1,sd)
  se.current = sd.current/sqrt(ncol(current.data))
  
  # insert into our array
  
  Exp_13_LPC.samp.pmol.mL.neg.mean[,3*i-2] = mean.current
  Exp_13_LPC.samp.pmol.mL.neg.mean[,3*i-1] = sd.current
  Exp_13_LPC.samp.pmol.mL.neg.mean[,3*i] = se.current
  
  # update column labels
  
  colnames(Exp_13_LPC.samp.pmol.mL.neg.mean)[(3*i-2)] =
    paste0(Exp_13_PC.unique.ttps[i],".mean")
  colnames(Exp_13_LPC.samp.pmol.mL.neg.mean)[(3*i-1)] =
    paste0(Exp_13_PC.unique.ttps[i],".sd")
  colnames(Exp_13_LPC.samp.pmol.mL.neg.mean)[(3*i)] =
    paste0(Exp_13_PC.unique.ttps[i],".se")
  
  
}

# exploratory plots, all FFA data by compound ID

for (i in 1:nrow(Exp_13_FFA.neg.samp.pmol.mL)) {
  
  for (j in 1:length(unique(Exp_13_FFA.metdat$Treatment.ID))) {
    
    if (j==1) {
      plot(Exp_13_FFA.metdat[Exp_13_FFA.metdat$Treatment.ID==unique(Exp_13_FFA.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
           Exp_13_FFA.neg.samp.pmol.mL[i,Exp_13_FFA.metdat$Treatment.ID==unique(Exp_13_FFA.metdat$Treatment.ID)[j]],pch=j,col=j,
           ylab="pmol/mL",
           xlab="Time",
           ylim=c(min(Exp_13_FFA.neg.samp.pmol.mL[i,]),max(Exp_13_FFA.neg.samp.pmol.mL[i,]))
      )
    }  else {
      points(Exp_13_FFA.metdat[Exp_13_FFA.metdat$Treatment.ID==unique(Exp_13_FFA.metdat$Treatment.ID)[j],c("Date.time.sample.collected")],
             Exp_13_FFA.neg.samp.pmol.mL[i,Exp_13_FFA.metdat$Treatment.ID==unique(Exp_13_FFA.metdat$Treatment.ID)[j]],pch=j,col=j)
    }
    
  }
  
  title(main = paste0("Exp 13: ",Exp_13_FFA.neg.samp[i,1]))
  legend(x="topright",
         legend=unique(Exp_13_FFA.metdat$Treatment.ID),
         pch=c(1:j),
         col=c(1:j),
         bty="n",
         cex=0.7,
         pt.cex=0.7)
  
}

# stats

# print name of experiment
print("Experiment 13 - Negative Ion Mode")

# subset to only first and last timepoint

Exp_13_FFA.neg.fl = Exp_13_FFA.neg.samp.pmol.mL[,(Exp_13_FFA.metdat$Date.time.sample.collected %in% unique(Exp_13_FFA.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_13_FFA.neg.metdat.fl = Exp_13_FFA.metdat[(Exp_13_FFA.metdat$Date.time.sample.collected %in% unique(Exp_13_FFA.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_13_FFA.neg.metdat.fl$ttp.ID = paste0(Exp_13_FFA.neg.metdat.fl$Treatment.ID,"_",Exp_13_FFA.neg.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# combine the +HB and -HB controls, since have no apparent difference
Exp_13_FFA.neg.metdat.fl$ttp.ID[Exp_13_FFA.neg.metdat.fl$ttp.ID %in% c("Dark_control_no_HB_2013-12-14 09:30:00",
                                                       "Dark_control_plus_HB_2013-12-14 09:30:00")] = c("Dark_control_2013-12-14 09:30:00")

Exp_13_LPC.neg.fl = Exp_13_LPC.samp.pmol.mL.neg[,(Exp_13_FFA.metdat$Date.time.sample.collected %in% unique(Exp_13_FFA.metdat$Date.time.sample.collected)[c(1,3)])]
Exp_13_LPC.neg.metdat.fl = Exp_13_FFA.metdat[(Exp_13_FFA.metdat$Date.time.sample.collected %in% unique(Exp_13_FFA.metdat$Date.time.sample.collected)[c(1,3)]),]
Exp_13_LPC.neg.metdat.fl$ttp.ID = paste0(Exp_13_FFA.neg.metdat.fl$Treatment.ID,"_",Exp_13_FFA.neg.metdat.fl$Date.time.sample.collected) # create a single treatment-timepoint ID

# combine the +HB and -HB controls, since have no apparent difference
Exp_13_LPC.neg.metdat.fl$ttp.ID[Exp_13_LPC.neg.metdat.fl$ttp.ID %in% c("Dark_control_no_HB_2013-12-14 09:30:00",
                                                                       "Dark_control_plus_HB_2013-12-14 09:30:00")] = c("Dark_control_2013-12-14 09:30:00")

# create subsets for a given moiety and perform ANOVA, then Tukey HSD

for (i in 1:nrow(Exp_13_FFA.neg.fl)) { # subset by moiety
  
  print(Exp_13_FFA.neg.samp[i,1]) # print name of this moiety
  
  Exp_13_FFA.neg.fl.subs = as.data.frame(cbind(as.numeric(Exp_13_FFA.neg.fl[i,]),Exp_13_FFA.neg.metdat.fl$ttp.ID))
  Exp_13_FFA.neg.fl.subs$V1 = as.numeric(as.character(Exp_13_FFA.neg.fl.subs$V1))
  colnames(Exp_13_FFA.neg.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  FFA.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_13_FFA.neg.fl.subs)
  print(anova(FFA.mod))
  FFA.aov = aov(FFA.mod)
  tukey = TukeyHSD(FFA.aov, conf.level = 0.95, "Treatment")
  print(tukey)
  
  
  # calculate & display mean differences ± SD for selected treatment pairs
  
  xbar_init = mean(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  sd_init = sd(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  se_init = sd_init/sqrt(length(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"]))
  
  xbar_final = mean(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("+ UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")
    
  xbar_final = mean(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("- UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("- UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_FFA.neg.fl.subs$Conc_pmol_mL[Exp_13_FFA.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB, + HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n\n")
  cat("+ UVB, + HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n\n")

  }


for (i in 1:nrow(Exp_13_LPC.neg.fl)) { # subset by moiety
  
  print(Exp_13_LPC.neg.samp[i,1]) # print name of this moiety
  
  Exp_13_LPC.neg.fl.subs = as.data.frame(cbind(as.numeric(Exp_13_LPC.neg.fl[i,]),Exp_13_LPC.neg.metdat.fl$ttp.ID))
  Exp_13_LPC.neg.fl.subs$V1 = as.numeric(as.character(Exp_13_LPC.neg.fl.subs$V1))
  colnames(Exp_13_LPC.neg.fl.subs) = c("Conc_pmol_mL","Treatment")
  
  FFA.mod = lm(Conc_pmol_mL ~ Treatment, data = Exp_13_LPC.neg.fl.subs)
  print(anova(FFA.mod))
  FFA.aov = aov(FFA.mod)
  tukey = TukeyHSD(FFA.aov, conf.level = 0.95, "Treatment")
  print(tukey)
  
  # calculate & display mean differences ± SD for selected treatment pairs
  
  xbar_init = mean(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  sd_init = sd(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  se_init = sd_init/sqrt(length(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"]))

    xbar_final = mean(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("+ UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"]))

    delta = xbar_final-xbar_init
    uncert = sqrt(se_final^2 + se_init^2)
    
  cat("- UVB, - HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")
  cat("- UVB, - HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n")

    xbar_final = mean(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  se_final = sd_final/sqrt(length(Exp_13_LPC.neg.fl.subs$Conc_pmol_mL[Exp_13_LPC.neg.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"]))

    delta = xbar_final-xbar_init
    uncert = sqrt(se_final^2 + se_init^2)
    
  cat("+ UVB, + HB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n\n")
  cat("+ UVB, + HB rate, pmol per mL per hr, mean ± uncertainty: ",delta/8.2," ± ",uncert/8.2,"\n\n")
  
}

##### Exp 13 - subplot with MDA (malondialdehyde) assay results

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")
MDA_assay_Exp13_14Dec13 = read.csv("data/nice/MDA_TBARS_assay_results/MDA_assay_Exp13_14Dec13.csv")

MDA_assay_Exp13_14Dec13$ttp.ID = paste0(MDA_assay_Exp13_14Dec13$Treatment,"_",MDA_assay_Exp13_14Dec13$Timepoint_GMT) # create a single treatment-timepoint ID

# calculate mean & SD for each set of replicates

# preallocate array

Exp_13_MDA.unique.ttps = unique(MDA_assay_Exp13_14Dec13$ttp.ID)

Exp_13_MDA.umol_L.mean = matrix(data = NA, 
                                          nrow = 1,
                                          ncol = 3*length(Exp_13_MDA.unique.ttps)
)

colnames(Exp_13_MDA.umol_L.mean) = rep("",3*length(Exp_13_MDA.unique.ttps))

# calculate stats, populate array

for (i in 1:length(Exp_13_MDA.unique.ttps)) {
  
  current.data = MDA_assay_Exp13_14Dec13[MDA_assay_Exp13_14Dec13$ttp.ID==Exp_13_MDA.unique.ttps[i],3]
  
    mean.current = mean(current.data)
    sd.current = sd(current.data)
    se.current = sd.current/sqrt(length(current.data))
    
    # insert into our array
    
    Exp_13_MDA.umol_L.mean[,3*i-2] = mean.current
    Exp_13_MDA.umol_L.mean[,3*i-1] = sd.current
    Exp_13_MDA.umol_L.mean[,3*i] = se.current
    
    # update column labels
    
    colnames(Exp_13_MDA.umol_L.mean)[(3*i-2)] =
      paste0(Exp_13_MDA.unique.ttps[i],".mean")
    colnames(Exp_13_MDA.umol_L.mean)[(3*i-1)] =
      paste0(Exp_13_MDA.unique.ttps[i],".sd")
    colnames(Exp_13_MDA.umol_L.mean)[(3*i)] =
      paste0(Exp_13_MDA.unique.ttps[i],".se")
    
    
  }
  
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Exp13_MDA_assay_barplot.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey","darkgrey")
treatments <- c("Dark","-UVB, -het. bact.","+UVB, -het. bact.","+UVB, +het. bact.")

bardata.mean <- matrix(NA,4,3)
bardata.mean[1,]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Control_no_hetbact_12/14/13 9:30.mean",
                                                                                                                                                                "Control_no_hetbact_12/14/13 13:40.mean",
                                                                                                                                                                "Control_no_hetbact_12/14/13 17:50.mean")]
bardata.mean[3,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Quartz_no_hetbact_12/14/13 19:30.mean",
                                                                                                                                                                   "Quartz_no_hetbact_12/14/13 13:40.mean",
                                                                                                                                                                   "Quartz_no_hetbact_12/14/13 17:50.mean")]
bardata.mean[4,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Quartz_plus_hetbact_12/14/13 9:30.mean",
                                                                                                                                                                   "Quartz_plus_hetbact_12/14/13 13:40.mean",
                                                                                                                                                                   "Quartz_plus_hetbact_12/14/13 17:50.mean")]
bardata.mean[2,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("EPA_no_hetbact_12/14/13 9:30.mean",
                                                                                                                                                                   "EPA_no_hetbact_12/14/13 13:40.mean",
                                                                                                                                                                   "EPA_no_hetbact_12/14/13 17:50.mean")]

bardata.se <- matrix(NA,4,3)
bardata.se[1,]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Control_no_hetbact_12/14/13 9:30.se",
                                                                                 "Control_no_hetbact_12/14/13 13:40.se",
                                                                                 "Control_no_hetbact_12/14/13 17:50.se")]
bardata.se[3,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Quartz_no_hetbact_12/14/13 19:30.se",
                                                                                    "Quartz_no_hetbact_12/14/13 13:40.se",
                                                                                    "Quartz_no_hetbact_12/14/13 17:50.se")]
bardata.se[4,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("Quartz_plus_hetbact_12/14/13 9:30.se",
                                                                                    "Quartz_plus_hetbact_12/14/13 13:40.se",
                                                                                    "Quartz_plus_hetbact_12/14/13 17:50.se")]
bardata.se[2,2:3]<-Exp_13_MDA.umol_L.mean[colnames(Exp_13_MDA.umol_L.mean) %in% c("EPA_no_hetbact_12/14/13 9:30.se",
                                                                                    "EPA_no_hetbact_12/14/13 13:40.se",
                                                                                    "EPA_no_hetbact_12/14/13 17:50.se")]

# scale up to pmol/mL
bardata.mean = bardata.mean*1000
bardata.se = bardata.se*1000

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,plotTop), ylab=paste("MDA (pmol/mL)"), xlab="Timepoint", names.arg=c("Initial", "+ 4 h", "+ 8 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topleft",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()
 
##### some basic calculations of oxidation state based on calculations, etc.

# Experiment 13

# preallocate a matrix for our calculations

Exp_13.ox.sums = as.data.frame(matrix(ncol = ncol(Exp_13_PC.samp.pmol.mL.norm),
                                            nrow = 6))
colnames(Exp_13.ox.sums) = colnames(Exp_13_PC.samp.pmol.mL.norm)
rownames(Exp_13.ox.sums) = c("pmol_mL_ox","pmol_mL_unox","ox_state_weighted",
                                   "total_PC_pmol_mL","frac_ox",
                             "frac_ox_weighted")

# cycle through, perform calculations

for (i in 1:ncol(Exp_13.ox.sums)) {
  
  Exp_13.ox.sums[c("pmol_mL_ox"),i] =
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation>=1,i])
  
  Exp_13.ox.sums[c("pmol_mL_unox"),i] =
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation==0,i])
  
  Exp_13.ox.sums[c("total_PC_pmol_mL"),i] =
    sum(Exp_13_PC.samp.pmol.mL.norm[,i])
  
  Exp_13.ox.sums[c("ox_state_weighted"),i] =
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation==1,i])*1 +
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation==2,i])*2 +
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation==3,i])*3 +
    sum(Exp_13_PC.samp.pmol.mL.norm[Exp_13_PC.samp$degree_oxidation==4,i])*4

  Exp_13.ox.sums[c("frac_ox"),i] =
    Exp_13.ox.sums[c("pmol_mL_ox"),i]/Exp_13.ox.sums[c("total_PC_pmol_mL"),i]
  
  Exp_13.ox.sums[c("frac_ox_weighted"),i] =
    Exp_13.ox.sums[c("ox_state_weighted"),i]/Exp_13.ox.sums[c("total_PC_pmol_mL"),i]
  
  
}

# now, evaluate statistical significance

# subset to only first and last timepoint

Exp_13.ox.sums.fl = Exp_13.ox.sums[,(Exp_13_PC.metdat$Date.time.sample.collected %in% unique(Exp_13_PC.metdat$Date.time.sample.collected)[c(1,3)])]

for (i in 1:nrow(Exp_13.ox.sums.fl)) { # subset by moiety
  
  print(rownames(Exp_13.ox.sums.fl)[i]) # print name of this moiety
  
  Exp_13.ox.sums.fl.subs = as.data.frame(cbind(as.numeric(Exp_13.ox.sums.fl[i,]),Exp_13_PC.metdat.fl$ttp.ID))
  Exp_13.ox.sums.fl.subs$V1 = as.numeric(as.character(Exp_13.ox.sums.fl.subs$V1))
  colnames(Exp_13.ox.sums.fl.subs) = c("Value","Treatment")
  
  PC.mod = lm(Value ~ Treatment, data = Exp_13.ox.sums.fl.subs)
  
  print(anova(PC.mod))
  PC.aov = aov(PC.mod)
  tukey = TukeyHSD(PC.aov, conf.level = 0.95)
  print(tukey)
  
  # calculate & display mean differences ± SD for selected treatment pairs
  
  xbar_init = mean(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  sd_init = sd(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Dark_control_2013-12-14 09:30:00"])
  
  xbar_final = mean(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Quartz_no_HB_2013-12-14 17:50:00"])
  
  delta = xbar_final-xbar_init
  uncert = sqrt(sd_final^2 + sd_init^2)
  
  cat("+ UVB, - HB vs. initial, mean ± SD: ",delta," ± ",uncert,"\n")
  
  xbar_final = mean(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="EPA_no_HB_2013-12-14 17:50:00"])
  
  delta = xbar_final-xbar_init
  uncert = sqrt(sd_final^2 + sd_init^2)
  
  cat("- UVB, - HB vs. initial, mean ± SD: ",delta," ± ",uncert,"\n")
  
  xbar_final = mean(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  sd_final = sd(Exp_13.ox.sums.fl.subs$Value[Exp_13.ox.sums.fl.subs$Treatment=="Quartz_plus_HB_2013-12-14 17:50:00"])
  
  delta = xbar_final-xbar_init
  uncert = sqrt(sd_final^2 + sd_init^2)
  
  cat("+ UVB, + HB vs. initial, mean ± SD: ",delta," ± ",uncert,"\n\n")
  
}

##### some basic calculations for Exp_03a

Exp_03a_PC.metdat.rerun$ttp.ID = paste0(Exp_03a_PC.metdat.rerun$Treatment.ID,"_",Exp_03a_PC.metdat.rerun$Date.time.sample.collected) # create a single treatment-timepoint ID

# calculate mean & SD for each set of replicates

# preallocate array

Exp_03a_PC.rerun.unique.ttps = unique(Exp_03a_PC.metdat.rerun$ttp.ID)

Exp_03a_PC.samp.pmol.mL.norm.rerun.mean = matrix(data = NA, 
                                          nrow = nrow(Exp_03a_PC.samp.pmol.mL.norm.rerun),
                                          ncol = 3*length(Exp_03a_PC.rerun.unique.ttps)
)

rownames(Exp_03a_PC.samp.pmol.mL.norm.rerun.mean) = rownames(Exp_03a_PC.samp.pmol.mL.norm.rerun)
colnames(Exp_03a_PC.samp.pmol.mL.norm.rerun.mean) = rep("",3*length(Exp_03a_PC.rerun.unique.ttps))

# calculate stats, populate array

for (i in 1:length(Exp_03a_PC.rerun.unique.ttps)) {
  
  current.data = Exp_03a_PC.samp.pmol.mL.norm.rerun[,Exp_03a_PC.metdat.rerun$ttp.ID==Exp_03a_PC.rerun.unique.ttps[i]]
  
  mean.current = apply(current.data,1,mean)
  sd.current = apply(current.data,1,sd)
  se.current = sd.current/sqrt(ncol(current.data))
  
  # insert into our array
  
  Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[,3*i-2] = mean.current
  Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[,3*i-1] = sd.current
  Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[,3*i] = se.current
  
  # update column labels
  
  colnames(Exp_03a_PC.samp.pmol.mL.norm.rerun.mean)[(3*i-2)] =
    paste0(Exp_03a_PC.rerun.unique.ttps[i],".mean")
  colnames(Exp_03a_PC.samp.pmol.mL.norm.rerun.mean)[(3*i-1)] =
    paste0(Exp_03a_PC.rerun.unique.ttps[i],".sd")
  colnames(Exp_03a_PC.samp.pmol.mL.norm.rerun.mean)[(3*i)] =
    paste0(Exp_03a_PC.rerun.unique.ttps[i],".se")
  
  
}
  