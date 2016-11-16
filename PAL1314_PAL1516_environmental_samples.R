# PAL1314_PAL1516_environmental_samples.R

# Purpose: Analysis of environmental samples from the PAL1314 and PAL1516 field
# seasons

# Created 11/13/2016 J.R.C.

### workspace preparation ####

# load necessary libraries

library(LOBSTAHS)
library(chemCal)

# set wd

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")

### load some necessary metadata; define functions ###

# some chemical data

DNPPE_mg_mL_1314 = 0.0565 # concentration DNPPE added in liquid/liquid extractions
# during PAL1314 Antarctic work, mg/mL
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

### standards from 20161107  ###

# QE003063-QE003073

# standards from 20161107, for PAL1314 & LMG1401 particulate data
# these didn't include any betaine standards

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/6IPL_Standards_20161107_pos.RData") # load processed standard data
IPLstd_pos_20161107.raw = getLOBpeaklist(x6IPL_Standards_20161107_pos) # generate peaklist

# extract standards for each species, and DNPPE
Std_peakareas.20161107 = IPLstd_pos_20161107.raw[
  IPLstd_pos_20161107.raw$compound_name %in% c("PG 32:0","PE 32:0","PC 32:0",
                                               "MGDG 36:0","SQDG 34:3","DGDG 36:4",
                                               "DNPPE"),]

rownames(Std_peakareas.20161107) = Std_peakareas.20161107$compound_name
Std_peakareas.20161107 = Std_peakareas.20161107[,14:25]
Std_peakareas.20161107 = Std_peakareas.20161107[order(colnames(Std_peakareas.20161107))]

# separate a QC from the standards
Std_peakareas.20161107_QC = Std_peakareas.20161107[,c(ncol(Std_peakareas.20161107))]
Std_peakareas.20161107 = Std_peakareas.20161107[,-c(ncol(Std_peakareas.20161107))]

# define quantities on column for standards (in pmol), assuming 20 uL injection,
# per HFF spreadsheet for current VML standards: 4/12/16, with DNPPE from 3/31/16

# will iterate so I don't make a mistake

# assumes user followed VML protocol for serial dilution of standards; molar
# quantities in the vectors below should be listed from lowest to highest

# DNPPE: 4k standard at 0.051 mg/mL
# MGDG: 16k pmol/mL in highest concentration standard

# create df, populate first row (and rows 2-3, for DNPPE)
Stds_20161107_oc = as.data.frame(matrix(NA,11,8))
colnames(Stds_20161107_oc) = c("pmol_mL_MGDG","pmol_oc_PG","pmol_oc_PE",
                               "pmol_oc_PC","pmol_oc_MGDG","pmol_oc_SQDG",
                               "pmol_oc_DGDG","pmol_oc_DNPPE")
Stds_20161107_oc[1,] = c(16000,250.130,252.326,252.326,318.839,501.626,254.179,0.000)
Stds_20161107_oc$pmol_oc_DNPPE[2:3] = c(0,116.56)
  
# fill out the rest of the matrix
for (i in 2:nrow(Stds_20161107_oc)) {
  
  Stds_20161107_oc[i,1:7] = Stds_20161107_oc[i-1,1:7]/2
  
  if (i>3) {
    
    Stds_20161107_oc[i,8] = Stds_20161107_oc[i-1,8]/2
    
  }
  
}

# fit standard curves using standard data, define breakpoints in cases where we will need two prediction ranges
# because the MS response was different 

# PG

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PG 32:0",1:9]
x = rev(Stds_20161107_oc$pmol_oc_PG)[1:9]
linfit_low.PG.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_PG),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PG 32:0",],
     pch="+",
     ylab = "Peak area, PG 32:0",
     xlab = "pmol o.c., PG 32:0")
points(rev(Stds_20161107_oc$pmol_oc_PG)[1:9],fitted(linfit_low.PG.20161107),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PG 32:0",c(9,11)]
x = rev(Stds_20161107_oc$pmol_oc_PG)[c(9,11)]
linfit_hi.PG.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_PG)[c(9,11)],fitted(linfit_hi.PG.20161107),col="blue",pch="+")

PG_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PG 32:0",9]

# PE

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PE 32:0",1:9]
x = rev(Stds_20161107_oc$pmol_oc_PE)[1:9]
linfit_low.PE.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_PE),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PE 32:0",],
     pch="+",
     ylab = "Peak area, PE 32:0",
     xlab = "pmol o.c., PE 32:0")
points(rev(Stds_20161107_oc$pmol_oc_PE)[1:9],fitted(linfit_low.PE.20161107),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PE 32:0",c(9,11)]
x = rev(Stds_20161107_oc$pmol_oc_PE)[c(9,11)]
linfit_hi.PE.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_PE)[c(9,11)],fitted(linfit_hi.PE.20161107),col="blue",pch="+")

PE_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PE 32:0",9]

# PC

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PC 32:0",1:10]
x = rev(Stds_20161107_oc$pmol_oc_PC)[1:10]
linfit_low.PC.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 10 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_PC),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PC 32:0",],
     pch="+",
     ylab = "Peak area, PC 32:0",
     xlab = "pmol o.c., PC 32:0")
points(rev(Stds_20161107_oc$pmol_oc_PC)[1:10],fitted(linfit_low.PC.20161107),col="red",pch="+")

# we will need some other fit for levels higher than ~ 125 pmol o.c.

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PC 32:0",10:11]
x = rev(Stds_20161107_oc$pmol_oc_PC)[10:11]
linfit_hi.PC.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_PC)[10:11],fitted(linfit_hi.PC.20161107),col="blue",pch="+")

PC_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="PC 32:0",10]

# MGDG

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="MGDG 36:0",1:9]
x = rev(Stds_20161107_oc$pmol_oc_MGDG)[1:9]
linfit_low.MGDG.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 10 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_MGDG),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="MGDG 36:0",],
     pch="+",
     ylab = "Peak area, MGDG 36:0",
     xlab = "pmol o.c., MGDG 36:0")
points(rev(Stds_20161107_oc$pmol_oc_MGDG)[1:9],fitted(linfit_low.MGDG.20161107),col="red",pch="+")

# we will need some other fit for levels higher than ~ 155 pmol o.c.

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="MGDG 36:0",10:11]
x = rev(Stds_20161107_oc$pmol_oc_MGDG)[10:11]
linfit_hi.MGDG.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_MGDG)[10:11],fitted(linfit_hi.MGDG.20161107),col="blue",pch="+")

MGDG_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="MGDG 36:0",10]

# SQDG

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="SQDG 34:3",1:9]
x = rev(Stds_20161107_oc$pmol_oc_SQDG)[1:9]
linfit_low.SQDG.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_SQDG),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="SQDG 34:3",],
     pch="+",
     ylab = "Peak area, SQDG 34:3",
     xlab = "pmol o.c., SQDG 34:3")
points(rev(Stds_20161107_oc$pmol_oc_SQDG)[1:9],fitted(linfit_low.SQDG.20161107),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="SQDG 34:3",c(9,11)]
x = rev(Stds_20161107_oc$pmol_oc_SQDG)[c(9,11)]
linfit_hi.SQDG.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_SQDG)[c(9,11)],fitted(linfit_hi.SQDG.20161107),col="blue",pch="+")

SQDG_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="SQDG 34:3",9]

# DGDG

# curve fitting & diagnostics
# just need one curve here, as long as we omit the second-highest level

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DGDG 36:4",c(1:9,11)]
x = rev(Stds_20161107_oc$pmol_oc_DGDG)[c(1:9,11)]
linfit_low.DGDG.20161107 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_DGDG),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DGDG 36:4",],
     pch="+",
     ylab = "Peak area, DGDG 36:4",
     xlab = "pmol o.c., DGDG 36:4")
points(rev(Stds_20161107_oc$pmol_oc_DGDG)[c(1:9,11)],fitted(linfit_low.DGDG.20161107),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

linfit_hi.DGDG.20161107 = linfit_low.DGDG.20161107
DGDG_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DGDG 36:4",9]

# DNPPE
# no DNPPE added at two highest standard levels, per VML lab SOP

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",1:9]
x = 1/(4 + (rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:9])*.025)
polyfit_low.DNPPE.20161107 = lm(as.numeric(y)~x) # fit a first-degee inverse polynomial model
plot(rev(Stds_20161107_oc$pmol_oc_DNPPE),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",],
     pch="+",
     ylab = "Peak area, DNPPE",
     xlab = "pmol o.c., DNPPE",
     ylim = c(0,3.5e+09))
points(rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:9],fitted(polyfit_low.DNPPE.20161107),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

polyfit_hi.DNPPE.20161107 = polyfit_low.DNPPE.20161107
DNPPE_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",8]

### standards from 20160421  ###

# QE001265-QE001276

# maybe not needed for environmental samples
# these were really run for the main batch of the liposome experiment data

### standards from 20161005  ###

# QE002850-QE002859

# these were the closest set of standards (temporally) to PAL1314 and LMG1401 
# dissolved phase environmental samples, the Marchetti Antarctic diatom extracts,
# the KM1605 UV-ox experiment, and all the PAL1516 particulate samples, including
# all SPE prefilters

# these standards contain DGTS/DGTA, but the DGTS/DGTA has a crazy low response factor
# also don't know exact concentration of DGTS/DGTA... must check w Helen

# Kevin Becker also ran TAG standards around the same time, will use them later
# TAG standards are QE002840-QE002847

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/6IPL_plus_DGTS_Standards_20161005_pos.RData") # load processed standard data
IPLstd_pos_20161005.raw = getLOBpeaklist(IPL_plus_DGTS_Standards_20161005) # generate peaklist

# extract standards for each species, and DNPPE
Std_peakareas.20161005 = IPLstd_pos_20161005.raw[
  IPLstd_pos_20161005.raw$compound_name %in% c("PG 32:0","PE 32:0","PC 32:0",
                                               "MGDG 36:0","SQDG 34:3","DGDG 36:4",
                                               "DGTS_DGTA 32:0","DNPPE"),]

rownames(Std_peakareas.20161005) = Std_peakareas.20161005$compound_name
Std_peakareas.20161005 = Std_peakareas.20161005[,13:23]
Std_peakareas.20161005 = Std_peakareas.20161005[order(colnames(Std_peakareas.20161005))]

# separate a QC from the standards
Std_peakareas.20161005_QC = Std_peakareas.20161005[,c(ncol(Std_peakareas.20161005))]
Std_peakareas.20161005 = Std_peakareas.20161005[,-c(ncol(Std_peakareas.20161005))]

# define quantities on column for standards (in pmol), assuming 20 uL injection,
# per HFF spreadsheet for current VML standards: 4/12/16, with DNPPE from 3/31/16

# will iterate so I don't make a mistake

# assumes user followed VML protocol for serial dilution of standards; molar
# quantities in the vectors below should be listed from lowest to highest

# these standards only went up to 8k pmol/mL MGDG

# DNPPE: 4k standard at 0.051 mg/mL
# MGDG: 8k pmol/mL in highest concentration standard

# create df, populate first row (and row 2, for DNPPE)
Stds_20161005_oc = as.data.frame(matrix(NA,10,9))
colnames(Stds_20161005_oc) = c("pmol_mL_MGDG","pmol_oc_PG","pmol_oc_PE",
                               "pmol_oc_PC","pmol_oc_MGDG","pmol_oc_SQDG",
                               "pmol_oc_DGDG","pmol_oc_DNPPE","pmol_oc_DGTS_DGTA")
Stds_20161005_oc[1,] = c(8000,125.0650,126.1630,126.1630,159.4195,250.8130,127.0895,0.000,125)
# 125 for DGTS/A is a placeholder right now; need to confirm w/Helen

Stds_20161005_oc$pmol_oc_DNPPE[2] = c(116.56)

# fill out the rest of the matrix
for (i in 2:nrow(Stds_20161005_oc)) {
  
  Stds_20161005_oc[i,c(1:7,9)] = Stds_20161005_oc[i-1,c(1:7,9)]/2
  
  if (i>2) {
    
    Stds_20161005_oc[i,8] = Stds_20161005_oc[i-1,8]/2
    
  }
  
}

# fit standard curves using standard data, define breakpoints in cases where we
# will need two prediction ranges because the MS response was different 

# PG

# curve fitting & diagnostics

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PG 32:0",c(1:7)]
#x = rev(Stds_20161005_oc$pmol_oc_PG)[c(1:8,10)]/(20+0.1*rev(Stds_20161005_oc$pmol_oc_PG)[c(1:8,10)])
#x = (rev(Stds_20161005_oc$pmol_oc_PG)[c(1:8,10)]^0.8)
#x = 8*(1-exp(-0.008*rev(Stds_20161005_oc$pmol_oc_PG)[c(1:8,10)])^1)
x = rev(Stds_20161005_oc$pmol_oc_PG)[c(1:7)]

linfit_low.PG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_PG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PG 32:0",],
     pch="+",
     ylab = "Peak area, PG 32:0",
     xlab = "pmol o.c., PG 32:0")
points(rev(Stds_20161005_oc$pmol_oc_PG)[c(1:7)],fitted(linfit_low.PG.20161005),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PG 32:0",c(7,8,10)]
x = rev(Stds_20161005_oc$pmol_oc_PG)[c(7,8,10)]
linfit_hi.PG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_PG)[c(7,8,10)],fitted(linfit_hi.PG.20161005),col="blue",pch="+")

PG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PG 32:0",7]

# PE

# curve fitting & diagnostics

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PE 32:0",c(1:8)]
x = rev(Stds_20161005_oc$pmol_oc_PE)[c(1:8)]
linfit_low.PE.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_PE),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PE 32:0",],
     pch="+",
     ylab = "Peak area, PE 32:0",
     xlab = "pmol o.c., PE 32:0")
points(rev(Stds_20161005_oc$pmol_oc_PE)[c(1:8)],fitted(linfit_low.PE.20161005),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PE 32:0",c(8,10)]
x = rev(Stds_20161005_oc$pmol_oc_PE)[c(8,10)]
linfit_hi.PE.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_PE)[c(8,10)],fitted(linfit_hi.PE.20161005),col="blue",pch="+")

PE_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PE 32:0",8]

# PC

# curve fitting & diagnostics

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",c(1:8,10)]
x = rev(Stds_20161005_oc$pmol_oc_PC)[c(1:8,10)]
linfit_low.PC.20161005 = lm(as.numeric(y)~x) # fit a linear model for the entire range,
# while skipping the 9th standard because something was wonky with it
plot(rev(Stds_20161005_oc$pmol_oc_PC),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",],
     pch="+",
     ylab = "Peak area, PC 32:0",
     xlab = "pmol o.c., PC 32:0")
points(rev(Stds_20161005_oc$pmol_oc_PC)[c(1:8,10)],fitted(linfit_low.PC.20161005),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

linfit_hi.PC.20161005 = linfit_low.PC.20161005
PC_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",9]

# MGDG

# curve fitting & diagnostics

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="MGDG 36:0",1:8]
x = rev(Stds_20161005_oc$pmol_oc_MGDG)[1:8]
linfit_low.MGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 10 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_MGDG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="MGDG 36:0",],
     pch="+",
     ylab = "Peak area, MGDG 36:0",
     xlab = "pmol o.c., MGDG 36:0")
points(rev(Stds_20161005_oc$pmol_oc_MGDG)[1:8],fitted(linfit_low.MGDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="MGDG 36:0",c(8:10)]
x = rev(Stds_20161005_oc$pmol_oc_MGDG)[8:10]
linfit_hi.MGDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_MGDG)[8:10],fitted(linfit_hi.MGDG.20161005),col="blue",pch="+")

MGDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="MGDG 36:0",8]

# SQDG

# curve fitting & diagnostics
# something appears to be very weird with the SQDG in these standards
# will generate one linear curve while ommitting the 8th and 9th points

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(1:7,10)]
x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)]
linfit_low.SQDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_SQDG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",],
     pch="+",
     ylab = "Peak area, SQDG 34:3",
     xlab = "pmol o.c., SQDG 34:3")
points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)],fitted(linfit_low.SQDG.20161005),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

linfit_hi.SQDG.20161005 = linfit_low.SQDG.20161005
SQDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",9]

# DGDG

# curve fitting & diagnostics
# just need one curve here, as long as we omit the second-highest level

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",c(1:8,10)]
x = 1/(4 + (rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8,10)])*.01)
hyperfit_low.DGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGDG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",],
     pch="+",
     ylab = "Peak area, DGDG 36:4",
     xlab = "pmol o.c., DGDG 36:4")
points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8,10)],fitted(hyperfit_low.DGDG.20161005),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

hyperfit_hi.DGDG.20161005 = hyperfit_low.DGDG.20161005
DGDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",9]

# DNPPE
# no DNPPE added at highest standard level, per VML lab SOP

# curve fitting & diagnostics
# will proceed with single linear fit after omitting 8th and 9th points

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DNPPE",c(1:7,9)]
x = rev(Stds_20161005_oc$pmol_oc_DNPPE)[c(1:7,9)]
linfit_low.DNPPE.20161005 = lm(as.numeric(y)~x) # fit a first-degee inverse polynomial model
plot(rev(Stds_20161005_oc$pmol_oc_DNPPE),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DNPPE",],
     pch="+",
     ylab = "Peak area, DNPPE",
     xlab = "pmol o.c., DNPPE",
     ylim = c(0,12e+08))
points(rev(Stds_20161005_oc$pmol_oc_DNPPE)[c(1:7,9)],fitted(linfit_low.DNPPE.20161005),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

linfit_hi.DNPPE.20161005= linfit_low.DNPPE.20161005
DNPPE_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DNPPE",8]

# DGTS & A

# curve fitting & diagnostics
# will fit a linear curve to first 8 points; hard to tell what's going on at
# greater concentrations

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",c(1:8)]
x = rev(Stds_20161005_oc$pmol_oc_DGTS_DGTA)[c(1:8)]
linfit_low.DGTS_DGTA.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGTS_DGTA),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",],
     pch="+",
     ylab = "Peak area, DGTS_DGTA 32:0",
     xlab = "pmol o.c., DGTS_DGTA 32:0")
points(rev(Stds_20161005_oc$pmol_oc_DGTS_DGTA)[c(1:8)],fitted(linfit_low.DGTS_DGTA.20161005),col="red",pch="+")

# will define the high-range curve for now as "NA," so we don't accidentally
# invoke it outside of its known linear range 

linfit_hi.DGTS_DGTA.20161005 = NA
DGTS_DGTA_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",8]

### pull in data ####

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/PAL1314_LMG1401_particulate_enviro_samples.RData")
PAL1314_LMG1401_0.2um_enviro_pos = getLOBpeaklist(PAL1314_LMG1401_particulate_enviro_samples_pos) # generate peaklist

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/PAL1314_LMG1401_dissolved_phase_enviro_samples,Marchetti_diatom_cultures_particulate.20161029.RData")
PAL1314_LMG1401_dissolved_pos = getLOBpeaklist(Pooled_environmental_samples_20161029) # generate peaklist

# split this one into two data frames, since it really contains two data sets
Marchetti_diatom_cultures_pos = PAL1314_LMG1401_dissolved_pos[,15:21]
PAL1314_LMG1401_dissolved_pos = PAL1314_LMG1401_dissolved_pos[,-c(15:21)]

