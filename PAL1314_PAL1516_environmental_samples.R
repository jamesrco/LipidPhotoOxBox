# PAL1314_PAL1516_environmental_samples.R

# Purpose: Analysis of environmental samples from the PAL1314 and PAL1516 field
# seasons

# Created 11/13/2016 J.R.C.

### workspace preparation ####

# load necessary libraries

library(LOBSTAHS)
library(chemCal)
library(stats)

# set wd

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")

### load some necessary metadata; define functions ###

# some chemical data

DNPPE_mg_mL_1314 = 0.0565 # concentration DNPPE added in liquid/liquid extractions
# during PAL1314 Antarctic work, mg/mL
DNPPE_mg_mL_BD_extracts_2016 = 0.051 # concentration of the DNPPE added during B&D extractions
# in Van Mooy Lab in 2016
DNPPE_MW = 875.081 # MW DNPPE, g/mol

# define a two-step function "splitpred" to compute pmol o.c. from raw values for PC and DNPPE, using one of two linear standard curves depending on magnitude of concentration (defined by cutoff value)

splitpred = function(x,linfit_low,linfit_hi,cutoff) {
  
  input = x
  
  if (!is.na(x)) {
    
    if (input > cutoff) {
      
      predval = inverse.predict(linfit_hi, newdata = x)
      
    } else if (input <= cutoff) {
      
      predval = inverse.predict(linfit_low, newdata = x)
      
    }
    
    predval$Prediction
    
  } else {
    
    NA
    
  }
  
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

### standards from 20161107  ####

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

# # DNPPE
# # no DNPPE added at two highest standard levels, per VML lab SOP
# 
# # curve fitting & diagnostics
# 
# y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",1:9]
# x = 1/(4 + (rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:9])*.025)
# polyfit_low.DNPPE.20161107 = lm(as.numeric(y)~x) # fit a first-degee inverse polynomial model
# plot(rev(Stds_20161107_oc$pmol_oc_DNPPE),
#      Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",],
#      pch="+",
#      ylab = "Peak area, DNPPE",
#      xlab = "pmol o.c., DNPPE",
#      ylim = c(0,3.5e+09))
# points(rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:9],fitted(polyfit_low.DNPPE.20161107),col="red",pch="+")
# 
# # can define the high-range curve to be the same as the low-range curve, in this case
# 
# polyfit_hi.DNPPE.20161107 = polyfit_low.DNPPE.20161107
# DNPPE_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",8]

# DNPPE
# no DNPPE added at two highest standard levels, per VML lab SOP

# curve fitting & diagnostics

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",1:8]
x = rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:8]
linfit_low.DNPPE.20161107 = lm(as.numeric(y)~x) # fit a linear model to first 8 points
plot(rev(Stds_20161107_oc$pmol_oc_DNPPE),
     Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",],
     pch="+",
     ylab = "Peak area, DNPPE",
     xlab = "pmol o.c., DNPPE",
     ylim = c(0,3.5e+09))
points(rev(Stds_20161107_oc$pmol_oc_DNPPE)[1:8],fitted(linfit_low.DNPPE.20161107),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

y = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",c(8,9)]
x = rev(Stds_20161107_oc$pmol_oc_DNPPE)[c(8,9)]
linfit_hi.DNPPE.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_DNPPE)[c(8,9)],fitted(linfit_hi.DNPPE.20161107),col="blue",pch="+")

DNPPE_std_breakpoint.20161107 = Std_peakareas.20161107[rownames(Std_peakareas.20161107)=="DNPPE",8]

### standards from 20160421  ####

# QE001265-QE001276

# maybe not needed for environmental samples
# these were really run for the main batch of the liposome experiment data

### standards from 20161005  ####

# QE002850-QE002859

# these were the closest set of standards (temporally) to PAL1314 and LMG1401 
# dissolved phase environmental samples, the Marchetti Antarctic diatom extracts,
# the KM1605 UV-ox experiment, and all the PAL1516 particulate samples, including
# all SPE prefilters

# these standards contain DGTS 16:0, 16:0 (from Avanti), but the DGTS has a
# crazy low response factor
# per HFF, DGTS was added to 4000 pmol/mL MGDG standard to achieve 4000 pmol/mL
# concentration 

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

# create df, populate first row (and row 2, for DNPPE and DGTS)
Stds_20161005_oc = as.data.frame(matrix(NA,10,9))
colnames(Stds_20161005_oc) = c("pmol_mL_MGDG","pmol_oc_PG","pmol_oc_PE",
                               "pmol_oc_PC","pmol_oc_MGDG","pmol_oc_SQDG",
                               "pmol_oc_DGDG","pmol_oc_DNPPE","pmol_oc_DGTS")
Stds_20161005_oc[1,] = c(8000,125.0650,126.1630,126.1630,159.4195,250.8130,127.0895,0.000,0.000)

Stds_20161005_oc$pmol_oc_DNPPE[2] = c(116.56)
Stds_20161005_oc$pmol_oc_DGTS[2] = c(79.7097500)

# fill out the rest of the matrix
for (i in 2:nrow(Stds_20161005_oc)) {
  
  Stds_20161005_oc[i,c(1:7)] = Stds_20161005_oc[i-1,c(1:7)]/2
  
  if (i>2) {
    
    Stds_20161005_oc[i,8:9] = Stds_20161005_oc[i-1,8:9]/2
    
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

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",c(1:8)]
x = rev(Stds_20161005_oc$pmol_oc_PC)[c(1:8)]
linfit_low.PC.20161005 = lm(as.numeric(y)~x) # fit a linear model for the entire range,
# while skipping the 9th standard because something was wonky with it
plot(rev(Stds_20161005_oc$pmol_oc_PC),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",],
     pch="+",
     ylab = "Peak area, PC 32:0",
     xlab = "pmol o.c., PC 32:0")
points(rev(Stds_20161005_oc$pmol_oc_PC)[c(1:8)],fitted(linfit_low.PC.20161005),col="red",pch="+")

# high range

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",c(8,10)]
x = rev(Stds_20161005_oc$pmol_oc_PC)[c(8,10)]
linfit_hi.PC.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_PC)[c(8,10)],fitted(linfit_hi.PC.20161005),col="blue",pch="+")

PC_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="PC 32:0",8]

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
# will generate one curve while ommitting the 8th and 9th points

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(1:7)]
# x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)]
x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7)]
linfit_low.SQDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 7 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_SQDG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",],
     pch="+",
     ylab = "Peak area, SQDG 34:3",
     xlab = "pmol o.c., SQDG 34:3")
points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7)],fitted(linfit_low.SQDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(7,10)]
x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(7,10)]
linfit_hi.SQDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(7,10)],fitted(linfit_hi.SQDG.20161005),col="blue",pch="+")

SQDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",7]

# # SQDG
# 
# # curve fitting & diagnostics
# # something appears to be very weird with the SQDG in these standards
# # will generate one curve while ommitting the 8th and 9th points
# 
# y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(1:7,10)]
# # x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)]
# x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)]*1/(4 + (rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)])*.01)
# hyperfit_low.SQDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
# plot(rev(Stds_20161005_oc$pmol_oc_SQDG),
#      Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",],
#      pch="+",
#      ylab = "Peak area, SQDG 34:3",
#      xlab = "pmol o.c., SQDG 34:3")
# points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)],fitted(hyperfit_low.SQDG.20161005),col="red",pch="+")
# 
# # can define the high-range curve to be the same as the low-range curve, in this case
# 
# hyperfit_hi.SQDG.20161005 = hyperfit_low.SQDG.20161005
# SQDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",9]

# xhat=c((4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(1)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(1)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(2)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(3)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(3)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(3)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(4)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(4)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(5)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(5)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(6)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(6)])[c("Prediction")])*0.01),
#        (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(7)])[c("Prediction")]))/
#          (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(7)])[c("Prediction")])*0.01),
# (4*as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(10)])[c("Prediction")]))/
#   (1-as.numeric(inverse.predict(hyperfit_low.SQDG.20161005,Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="SQDG 34:3",c(10)])[c("Prediction")])*0.01))

# points(xhat,fitted(hyperfit_low.SQDG.20161005),col="green",pch="o")

# # DGDG
# 
# # curve fitting & diagnostics
# # just need one curve here, as long as we omit the second-highest level
# 
# y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",c(1:8,10)]
# x = rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8,10)]*1/(4 + (rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8,10)])*.01)
# hyperfit_low.DGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 9 standard levels
# plot(rev(Stds_20161005_oc$pmol_oc_DGDG),
#      Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",],
#      pch="+",
#      ylab = "Peak area, DGDG 36:4",
#      xlab = "pmol o.c., DGDG 36:4")
# points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8,10)],fitted(hyperfit_low.DGDG.20161005),col="red",pch="+")
# 
# # can define the high-range curve to be the same as the low-range curve, in this case
# 
# hyperfit_hi.DGDG.20161005 = hyperfit_low.DGDG.20161005
# DGDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",9]

# DGDG

# curve fitting & diagnostics

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",c(1:8)]
x = rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8)]
linfit_low.DGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 7 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGDG),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",],
     pch="+",
     ylab = "Peak area, DGDG 36:4",
     xlab = "pmol o.c., DGDG 36:4")
points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8)],fitted(linfit_low.DGDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",c(8,10)]
x = rev(Stds_20161005_oc$pmol_oc_DGDG)[c(8,10)]
linfit_hi.DGDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(8,10)],fitted(linfit_hi.DGDG.20161005),col="blue",pch="+")

DGDG_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGDG 36:4",8]

# DNPPE
# no DNPPE added at highest standard level, per VML lab SOP

# curve fitting & diagnostics
# will proceed with single linear fit after omitting 8th and 9th points

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DNPPE",c(1:7,9)]
x = rev(Stds_20161005_oc$pmol_oc_DNPPE)[c(1:7,9)]
linfit_low.DNPPE.20161005 = lm(as.numeric(y)~x) # fit linear model
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

# DGTS

# curve fitting & diagnostics

y = as.numeric(Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",c(1:8)])
x = rev(Stds_20161005_oc$pmol_oc_DGTS)[c(1:8)]

linfit_low.DGTS_DGTA.20161005 = lm(as.numeric(y)~x) # fit a model for the first 8 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGTS),
     Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",],
     pch="+",
     ylab = "Peak area, DGTS_DGTA 32:0",
     xlab = "pmol o.c., DGTS_DGTA 32:0")
points(rev(Stds_20161005_oc$pmol_oc_DGTS)[c(1:8)],fitted(linfit_low.DGTS_DGTA.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",c(8:9)]
x = rev(Stds_20161005_oc$pmol_oc_DGTS)[8:9]
linfit_hi.DGTS_DGTA.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_DGTS)[8:9],fitted(linfit_hi.DGTS_DGTA.20161005),col="blue",pch="+")

DGTS_DGTA_std_breakpoint.20161005 = Std_peakareas.20161005[rownames(Std_peakareas.20161005)=="DGTS_DGTA 32:0",8]

### TAG standards from 20161004 ####

# run by Kevin Becker, QE002840-QE002847
# the standard mix includes odd fatty acid TAGs as well, so re-ran the dataset
# in LOBSTAHS with exclude.oddFA = F
# standards were run in duplicate, will take averages

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/TAG_Standards_20161004_pos.RData") # load processed standard data
TAGstd_pos_20161004.raw = getLOBpeaklist(TAG_Standards_20161004_pos) # generate peaklist

# extract standards for each TAG in the mix, and DNPPE
TAGstd_peakareas.20161004 = TAGstd_pos_20161004.raw[
  TAGstd_pos_20161004.raw$compound_name %in% c("TAG 24:0","TAG 27:0","TAG 30:0",
                                               "TAG 33:0","TAG 36:0","TAG 39:0",
                                               "TAG 42:0","TAG 45:0","TAG 48:0",
                                               "TAG 51:0","TAG 54:0","TAG 57:0",
                                               "TAG 48:3","TAG 54:9","TAG 60:15",
                                               "TAG 66:18","DNPPE"),]

# three of these (24:0, 27:0, 57:0) aren't in the current LOBSTAHS DB, not going
# to worry about them

rownames(TAGstd_peakareas.20161004) = TAGstd_peakareas.20161004$compound_name
TAGstd_peakareas.20161004 = TAGstd_peakareas.20161004[,13:20]
TAGstd_peakareas.20161004 = TAGstd_peakareas.20161004[order(colnames(TAGstd_peakareas.20161004))]

# # separate a QC from the standards
# no QC's in these standards, skipping... 
# TAGstd_peakareas.20161004_QC = TAGstd_peakareas.20161004[,c(ncol(TAGstd_peakareas.20161004))]
# TAGstd_peakareas.20161004 = TAGstd_peakareas.20161004[,-c(ncol(TAGstd_peakareas.20161004))]

# average the respective duplicates

TAGstd_peakareas.20161004.means = as.data.frame(matrix(NA,nrow(TAGstd_peakareas.20161004),4))

colnames(TAGstd_peakareas.20161004.means) = c("ng_oc_0.5","ng_oc_4","ng_oc_15","ng_oc_40")
rownames(TAGstd_peakareas.20161004.means) = rownames(TAGstd_peakareas.20161004)

TAGstd_peakareas.20161004.means[,1] = apply(TAGstd_peakareas.20161004[,1:2],1,mean)
TAGstd_peakareas.20161004.means[,2] = apply(TAGstd_peakareas.20161004[,3:4],1,mean)
TAGstd_peakareas.20161004.means[,3] = apply(TAGstd_peakareas.20161004[,5:6],1,mean)
TAGstd_peakareas.20161004.means[,4] = apply(TAGstd_peakareas.20161004[,7:8],1,mean)

# put in more logical order
TAGstd_peakareas.20161004.means =
  TAGstd_peakareas.20161004.means[order(rownames(TAGstd_peakareas.20161004.means)),]

# define quantities on column for standards (in pmol), per Kevin's notes, with
# DNPPE from 3/31/16

# specify molecular weights for each species, in order of species as they appear
# in TAGstd_peakareas.20161004.means after reordering, just above

TAG.mws_20161004 = c(857.51666,554.45464,596.50159,638.54854,680.59549,722.64244,
                     764.68939,806.73634,800.68939,848.78329,890.83024,872.68939,
                     944.68939,1022.73634)

# DNPPE: assumed added at same concentration (ng) as the TAGs

# calculate pmol o.c. for each of the four standard levels, populate data frame

TAGstds_20161004_pmol_oc = as.data.frame(matrix(NA,ncol(TAGstd_peakareas.20161004.means),
                                                nrow(TAGstd_peakareas.20161004.means)))
colnames(TAGstds_20161004_pmol_oc) = apply(expand.grid(c("pmol_oc_"),rownames(TAGstd_peakareas.20161004.means)),1,paste,collapse="")
colnames(TAGstds_20161004_pmol_oc) = gsub(" ","_",colnames(TAGstds_20161004_pmol_oc))
colnames(TAGstds_20161004_pmol_oc) = gsub(":","_",colnames(TAGstds_20161004_pmol_oc))

rownames(TAGstds_20161004_pmol_oc) = colnames(TAGstd_peakareas.20161004.means)

TAGstds_20161004_pmol_oc[1,] = 0.5/TAG.mws_20161004*1000
TAGstds_20161004_pmol_oc[2,] = 4/TAG.mws_20161004*1000
TAGstds_20161004_pmol_oc[3,] = 15/TAG.mws_20161004*1000
TAGstds_20161004_pmol_oc[4,] = 40/TAG.mws_20161004*1000

# now, build standard curves for each TAG, and DNPPE

# preallocate a list object to hold results of the curve fits, set correct element names

TAGstdlist.names = rownames(TAGstd_peakareas.20161004.means)
TAGstdlist.names = gsub(" ","_",TAGstdlist.names)
TAGstdlist.names = gsub(":","_",TAGstdlist.names)

TAGstds_20161004_linfits = vector("list",length(TAGstdlist.names))
names(TAGstds_20161004_linfits) = TAGstdlist.names

# apply linear fit to each of these, then plot, then store

for (i in 1:length(TAGstds_20161004_linfits)) {
  
  # for a few species, want to omit some bad data points in the standards
  
  if (rownames(TAGstd_peakareas.20161004.means)[i] %in% c("TAG 48:0","TAG 51:0")) {
    
    y = TAGstd_peakareas.20161004.means[i,c(1,3,4)]
    x = TAGstds_20161004_pmol_oc[c(1,3,4),i]
    
  } else if (rownames(TAGstd_peakareas.20161004.means)[i]=="TAG 54:0") {
    
    y = TAGstd_peakareas.20161004.means[i,c(1,2,4)]
    x = TAGstds_20161004_pmol_oc[c(1,2,4),i]
    
  } else {
    
    y = TAGstd_peakareas.20161004.means[i,]
    x = TAGstds_20161004_pmol_oc[,i]
    
  }
  
  linfit = lm(as.numeric(y)~x-1) # fit a linear model, force through origin
  plot(TAGstds_20161004_pmol_oc[,i],
       TAGstd_peakareas.20161004.means[i,],
       pch="+",
       ylab = paste0("Peak area, ",rownames(TAGstd_peakareas.20161004.means)[i]),
       xlab = paste0("pmol o.c., ",rownames(TAGstd_peakareas.20161004.means)[i]))
  points(x,fitted(linfit),col="red",pch="+")
  
  TAGstds_20161004_linfits[[i]] = linfit
  
}

# can generate some relative response factors

RRFs = vector("numeric",length(TAGstds_20161004_linfits)) # preallocate
names(RRFs) = rownames(TAGstd_peakareas.20161004.means)

for (i in 1:length(RRFs)) {
  
  RRFs[i] = TAGstds_20161004_linfits[[i]]$coefficients[1]/
    TAGstds_20161004_linfits[[1]]$coefficients[1]
  
}

# calculate some equivalent carbon #s for the TAGs, then plot with RRFs 

ECNs = vector("numeric",nrow(TAGstd_peakareas.20161004.means)-1) # preallocate
names(ECNs) = rownames(TAGstd_peakareas.20161004.means)[2:nrow(TAGstd_peakareas.20161004.means)]

for (i in 1:length(ECNs)) {
  
  # get no. of C, DB
  
  num_C = TAGstd_pos_20161004.raw$FA_total_no_C[TAGstd_pos_20161004.raw$compound_name==names(ECNs)[i]]
  num_DB = TAGstd_pos_20161004.raw$FA_total_no_DB[TAGstd_pos_20161004.raw$compound_name==names(ECNs)[i]]
    
  ECNs[i] = num_C-2*num_DB

}

plot(ECNs,RRFs[2:length(RRFs)])
text(ECNs, RRFs[2:length(RRFs)], labels=names(ECNs), pos = 4, cex = 0.5)

# # maybe a plot of the reciprocals
# 
# plot(ECNs,1/RRFs[2:length(RRFs)])
# text(ECNs, 1/RRFs[2:length(RRFs)], labels=names(ECNs), pos = 4, cex = 0.5)
# 
# x = ECNs
# y = 1/RRFs[2:length(RRFs)]
# 
# # fit a non-linear model
# 
# nls.fit = nls(y~a/x+b*x+c,list(x,y),c(a=700,b=0.5,c=-40)) # fit a linear model, force through origin
# points(x,fitted(nls.fit),col="red",pch="+")

### Marchetti diatom cultures ####

# will use the 20161005 standards for these files

## pull in data ###

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/UNC_Marchetti_diatom_cultures_pos_withoddFA_LOBSet.RData")
Marchetti_diatom_cultures_pos_withoddFA = getLOBpeaklist(Marchetti_diatom_cultures_pos_withoddFA) # generate peaklist

# first, extract only unoxidized IPL (no TAGs, etc) data, plus DNPPE

Marchetti_diatom_cultures_pos.unox_IPL = Marchetti_diatom_cultures_pos_withoddFA[
  (Marchetti_diatom_cultures_pos_withoddFA$lipid_class %in% c("IP_DAG","DNPPE") & 
     Marchetti_diatom_cultures_pos_withoddFA$degree_oxidation==0),]

## ***** experimental trial of use of msn data using series of xcmsRaw objects ***** ##

# now, use fragmentation spectra for basic confirmation of putative LOBSTAHS IDs

# requires several additional files: (1) annotated xcmsSet object for data in positive ion mode, (2) LOBSet object in positive ion mode, and (3) list containing positive-mode xcmsRaw objects for all samples generated using: 
#             xraw <- xcmsRaw("yourfile.mzXML", includeMSn=TRUE)
# the routine below assumes the objects are stored in the list in the order in which they appear in the xsAnnotate and LOBSet objects

# first, define types and values of MS fragmentation experiments for each IP-DAG class
# (i.e., constant neutral loss (CNL) or product ion (PI))

# create empty data frame
Marchetti_diatom.frag_lookup_classes = as.data.frame(matrix(NA,8,2))
rownames(Marchetti_diatom.frag_lookup_classes) = 
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGCC","DGTS_DGTA")
colnames(Marchetti_diatom.frag_lookup_classes) =
  c("Frag_exp_type","mz_value")

# now, populate our data frame with necessary values
# per Popendorf et al., Lipids (2013) 48:185–195

Marchetti_diatom.frag_lookup_classes[,1] =
  c("CNL","CNL","PI","CNL","CNL","CNL","PI","PI")
Marchetti_diatom.frag_lookup_classes[,2] =
  c(189.040224,141.019094,184.073321,197.089937,359.142212,261.051837,104.106690,236.149249)

# pull in the xsAnnotate object & list containing the xcmsRaw objects

load("data/nice/Orbi_MS_data/xsAnnotate_objects/UNC_Marchetti_diatom_cultures_pos_withoddFA_xsAnnotate.RData")
Marchetti_diatoms_xsA_pos = UNC_Marchetti_diatom_cultures_pos_withoddFA_xsAnnotate # rename so easier to work with

load("data/raw/Orbi_MS_data/xcmsRaw_objects/UNC_Marchetti_diatom_cultures_pos_xcmsRaw.RData")
Marchetti_xsR = UNC_Marchetti_diatom_cultures_pos_xcmsRaw
  
# preallocate three matrices for our results
# Marchetti_diatom.detected_pos_ion_fragSpec: to keep track of how many valid ms2 fragmentation spectra were detected for the feature in positive ion mode
# Marchetti_diatom.fragdata_results: number of ms2 spectra in which the PI or CNL criteria were validated

Marchetti_diatom.detected_pos_ion_fragments = as.data.frame(matrix(NA,nrow(Marchetti_diatom_cultures_pos.unox_IPL),ncol=7))
Marchetti_diatom.fragdata_results = as.data.frame(matrix(NA,nrow(Marchetti_diatom_cultures_pos.unox_IPL),ncol=7))

# necessary functions that will allow us to extract the correct ms2 spectra, evaluate transitions, etc

get.ms2Peaklist = function (precursor.index,sample_ID) {
  
  ms2data.start = Marchetti_xsR[[sample_ID]]@msnScanindex[precursor.index]
  ms2data.end = Marchetti_xsR[[sample_ID]]@msnScanindex[precursor.index+1]-1
  
  scandata =
    data.frame(Marchetti_xsR[[sample_ID]]@env$msnMz[ms2data.start:ms2data.end],
               Marchetti_xsR[[sample_ID]]@env$msnIntensity[ms2data.start:ms2data.end])
  
  colnames(scandata) = c("mz","Intensity")
  
  return(scandata)
  
}

get.topN = function(peaklist,N) {
  
  ordered.peaklist = peaklist[order(peaklist$Intensity, decreasing = TRUE),]
  
  topN.peaklist = ordered.peaklist[1:N,]
  
  return(topN.peaklist)
  
}

eval.PIspecies = function(peaklist,species,ppm) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # retrieve product ion m/z
  
  prod.ion = Marchetti_diatom.frag_lookup_classes$mz_value[
    rownames(Marchetti_diatom.frag_lookup_classes)==species]
  
  if (any(abs((prod.ion-peaklist[,1])/prod.ion*1000000)<ppm)) {
    
    # it's a match
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}

eval.CNLspecies = function(peaklist,species,sample_ID,ppm) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # calculate theoretical m/z of the ion that would be produced via the neutral loss
  # throwing in a mean() here in case xcms associated more than one peak with the group in this particular sample
  CNL_product_mz = mean(xcms.peakdata_thisgroup_pos[xcms.peakdata_thisgroup_pos$sample==sample_ID,1])-
    Marchetti_diatom.frag_lookup_classes$mz_value[
      rownames(Marchetti_diatom.frag_lookup_classes)==this.IDclass]
  
  # perform comparison
  
  if (any(abs((CNL_product_mz-peaklist[,1])/CNL_product_mz*1000000)<ppm)) {
    
    # it's a match
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}

sumfrag = function(x) {
  
  if (all(is.na(x))) {
    
    NA
    
  } else {
    
    sum(x, na.rm = T)
    
  }
  
}

# iterate through the IP-DAG IDs by sample, retrieve necessary data from the xsAnnotate and xcmsRaw objects, evaluate, and record the results

for (i in 1:(nrow(Marchetti_diatom_cultures_pos.unox_IPL))) {
  # iterate through each LOBSTAHS ID
  
  # retrieve LOBSTAHS compound ID, lipid species
  
  this.ID = Marchetti_diatom_cultures_pos.unox_IPL$compound_name[i]
  this.IDclass = Marchetti_diatom_cultures_pos.unox_IPL$species[i]
  
  if (this.IDclass %in% rownames(Marchetti_diatom.frag_lookup_classes)) {
    # an escape if the lipid class isn't accounted for in our input parameter table
    
    # retrieve underlying positive-mode xcms group and peak data
    
    xcms.peakIDs_thisgroup_pos =
      Marchetti_diatoms_xsA_pos@xcmsSet@groupidx[[Marchetti_diatom_cultures_pos.unox_IPL$xcms_peakgroup[i]]]
    
    xcms.peakdata_thisgroup_pos =
      as.data.frame(Marchetti_diatoms_xsA_pos@xcmsSet@peaks[xcms.peakIDs_thisgroup_pos,])
    
    # now, iterate through all instances of this putatively identified compound that are present in the xcms peakgroup and in one of the samples of interest (i.e., not a QC) 
    
    for (j in 1:nrow(xcms.peakdata_thisgroup_pos)) {
      
      # retrieve the sample ID
      
      samp_ID = xcms.peakdata_thisgroup_pos$sample[j]
      
      # before beginning, check if the peak is a QC; if so, skip to end
      
      if (samp_ID %in% c(1:7)) {
        
        ### first, pull out fragmentation spectra relevant to this peak, if they exist ###
        
        # retrieve raw (uncorrected) retention times for this peak
        
        RT_min_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmin[j]]
        
        RT_max_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmax[j]]
        
        RT_ctr_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rt[j]]
        
        # retrieve data from ms2 scans which have the correct precursor mz (i.e., the mz of this instance of the LOBSTAHS-ID'd compound (peak) we are currently considering) AND were acquired within the RT window (raw min/max) of the peak
        
        # will use search window of 4 ppm
        ms2_lkup_window.ppm = 4
        
        # first, gather possibly relevant precursor scans based strictly on mass difference
        
        # get an index to these scans
        possible.precursors_pos.ind =
          which(abs(Marchetti_xsR[[samp_ID]]@msnPrecursorMz-xcms.peakdata_thisgroup_pos$mz[j])/
                  xcms.peakdata_thisgroup_pos$mz[j]*1000000 < ms2_lkup_window.ppm)
        
        # whittle this list to make sure the scans *also* fall within the rt window for the parent peak
        # will use a little bit of a buffer to capture any scans falling just outside of the RT range
        
        valid.precursors_pos.ind =
          possible.precursors_pos.ind[Marchetti_xsR[[samp_ID]]@msnRt[possible.precursors_pos.ind]>
                                        (RT_min_pos.raw-10) &
                                        Marchetti_xsR[[samp_ID]]@msnRt[possible.precursors_pos.ind]<
                                        (RT_max_pos.raw+10)]
        
        # to get the actual QE scan numbers corresponding to this index, run the below:
        # Marchetti_xsR[[j]]@msnAcquisitionNum[valid.precursors_pos.ind]
        
        # record the number of valid ms2 spectra in Marchetti_diatom.detected_pos_ion_fragments
        # add to record if there's information already recorded from another isomer; otherwise, simply overwrite the NA placeholder
        
        if (length(valid.precursors_pos.ind)>0) {
          
          if (is.na(Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID])) {
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = length(valid.precursors_pos.ind)
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 
              Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] + length(valid.precursors_pos.ind)
            
          }
          
        } else {
          
          # there is no ms2 data for this parent, at least not how we went about finding it
          
          if (is.na(Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID])) {
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 0
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 
              Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] + 0
            
          }
          
        }
        
        # provide ourselves with an escape here in case there are no valid ms2 scans
        
        if (length(valid.precursors_pos.ind)>0) {
          
          # now, pull out the ms2 data for each of these scans
          # msnScanindex is constructed in such a way that we can recreate the peaklist for
          # a given scan using the following syntax
          
          relevant_ms2data = apply(as.matrix(valid.precursors_pos.ind),1,get.ms2Peaklist,samp_ID)
          
          ### now, can get onto the business of actually examining the spectra for the diagnostic transitions ###
          
          ### scenario 1: class type is diagnosed via presence of product ion ###
          
          if (Marchetti_diatom.frag_lookup_classes$Frag_exp_type[
            rownames(Marchetti_diatom.frag_lookup_classes)==this.IDclass] == "PI") {
            # this class is diagnosed via presence of a product ion

            # apply some logic for product ion scenario: assume that for a product ion-based ID to be "good", we must observe a feature with the mz of the diagnostic ion (+/- some mz tolerance) as one of the top N peaks (by intensity) in at least one of the relevant ms2 scans
            
            # so, extract the top N (right now, 10) most intense features in each scan
            
            top10.features = lapply(relevant_ms2data,get.topN,10)
            
            # evaluate: do the list(s) of the top N most intense fragments contain the diagnostic ion?
            
            PI.eval_result = lapply(top10.features, eval.PIspecies, species = this.IDclass, ppm = 4)
            
            # record result
            
            if (is.na(Marchetti_diatom.fragdata_results[i,samp_ID])) {
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = sum(unlist(PI.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = 
                Marchetti_diatom.fragdata_results[i,samp_ID] + sum(unlist(PI.eval_result))
              
            }
            
          } else if (Marchetti_diatom.frag_lookup_classes$Frag_exp_type[
            rownames(Marchetti_diatom.frag_lookup_classes)==this.IDclass] == "CNL") {
            
            ### scenario 2: class type is diagnosed via constant neutral loss ###
            
            # assume it's a good ID in this case as long as an ion corresponding to the diagnostic CNL is one of the top N (20, for now) peaks (by intensity) in the + mode ms2 spectrum
            
            top20.features = lapply(relevant_ms2data,get.topN,20)
            
            # evaluate & record
            
            CNL.eval_result = lapply(top20.features, eval.CNLspecies, species = this.IDclass, sample_ID = samp_ID, ppm = 4)
            
            if (is.na(Marchetti_diatom.fragdata_results[i,samp_ID])) {
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = sum(unlist(CNL.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = 
                Marchetti_diatom.fragdata_results[i,samp_ID] + sum(unlist(CNL.eval_result))
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
}


# append summary, presence/absence results to Marchetti_diatom_cultures_pos.unox_IPL matrix

Marchetti_diatom.fragdata_result.summary =
  apply(Marchetti_diatom.fragdata_results, 1, sumfrag)
Marchetti_diatom_cultures_pos.unox_IPL$ms2_conf =
  Marchetti_diatom.fragdata_result.summary

Marchetti_diatom.fragdata_presence.summary =
  apply(Marchetti_diatom.detected_pos_ion_fragments, 1, sumfrag)
Marchetti_diatom_cultures_pos.unox_IPL$ms2_present =
  Marchetti_diatom.fragdata_presence.summary

# # convert anything < 1e5 intensity to NA, assuming it's either noise
# # or something so low in concentraton as to be irrelevant from a total lipid
# # perspective
# 
# Marchetti_diatom_cultures_pos.unox_IPL[,14:20] =
#   replace(Marchetti_diatom_cultures_pos.unox_IPL[,14:20],
#           Marchetti_diatom_cultures_pos.unox_IPL[,14:20]<10000,
#           NA)

# get rid of DGCC data since we don't have any standards for these right now

Marchetti_diatom_cultures_pos.unox_IPL.noDGCC = Marchetti_diatom_cultures_pos.unox_IPL[
  Marchetti_diatom_cultures_pos.unox_IPL$species!="DGCC",]

# # assuming only odd-chain FA will be a C15 (based on knowledge of diatom FA biosynthesis),
# # can eliminate some putative odd-chain IDs
# 
# for (i in 1:nrow(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC)) {
#   
#   if (!is.na(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC$FA_total_no_C[i])) {
#     
#     if (!(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC$FA_total_no_C[i] %% 2 == 0)) {
#       
#       # we have an odd-chain species, evaluate
#       
#       if (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC$FA_total_no_C[i]>37) {
#         
#         # i.e., >= 22 + 15 
#         
#         # bulk C number greater than this would be an impossible combination,
#         # given what we know about diatom FA biosynthesis
#         
#         # so, set these to zero
#         
#         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC[i,c(14:20)] = 0
#         
#       }
#       
#     }
#     
#   }
#   
# }

# define classes for which concentrations are to be calculated, the models, and
# the cutoffs in case of split prediction

Marchetti_diatom.conc_classes = as.data.frame(matrix(NA,8,4))
colnames(Marchetti_diatom.conc_classes) =
  c("Lipid_class","Model.low","Model.hi","Cutoff_PA")

Marchetti_diatom.conc_classes[,1] =
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGTS_DGTA","DNPPE")

Marchetti_diatom.conc_classes[,2] =
  c("linfit_low.PG.20161005","linfit_low.PE.20161005","linfit_low.PC.20161005",
    "linfit_low.MGDG.20161005","linfit_low.SQDG.20161005",
    "linfit_low.DGDG.20161005","linfit_low.DGTS_DGTA.20161005",
    "linfit_low.DNPPE.20161005")

Marchetti_diatom.conc_classes[,3] =
  c("linfit_hi.PG.20161005","linfit_hi.PE.20161005","linfit_hi.PC.20161005",
    "linfit_hi.MGDG.20161005","linfit_hi.SQDG.20161005",
    "linfit_hi.DGDG.20161005","linfit_hi.DGTS_DGTA.20161005",
    "linfit_hi.DNPPE.20161005")

Marchetti_diatom.conc_classes[,4] =
c(PG_std_breakpoint.20161005,PE_std_breakpoint.20161005,PC_std_breakpoint.20161005,
MGDG_std_breakpoint.20161005,SQDG_std_breakpoint.20161005,DGDG_std_breakpoint.20161005,
DGTS_DGTA_std_breakpoint.20161005,DNPPE_std_breakpoint.20161005)
# 
# Marchetti_diatom.conc_classes[,2] =
#   c("linfit_low.PG.20161107","linfit_low.PE.20161107","linfit_low.PC.20161107",
#     "linfit_low.MGDG.20161107","linfit_low.SQDG.20161107",
#     "linfit_low.DGDG.20161107","linfit_low.DGTS_DGTA.20161005",
#     "linfit_low.DNPPE.20161107")
# 
# Marchetti_diatom.conc_classes[,3] =
#   c("linfit_hi.PG.20161107","linfit_hi.PE.20161107","linfit_hi.PC.20161107",
#     "linfit_hi.MGDG.20161107","linfit_hi.SQDG.20161107",
#     "linfit_hi.DGDG.20161107","linfit_hi.DGTS_DGTA.20161005",
#     "linfit_hi.DNPPE.20161107")
# 
# Marchetti_diatom.conc_classes[,4] =
#   c(PG_std_breakpoint.20161107,PE_std_breakpoint.20161107,PC_std_breakpoint.20161107,
#     MGDG_std_breakpoint.20161107,SQDG_std_breakpoint.20161107,DGDG_std_breakpoint.20161107,
#     DGTS_DGTA_std_breakpoint.20161005,DNPPE_std_breakpoint.20161107)

# first, calculate pmol o.c.
# preallocate matrix for result
Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc = Marchetti_diatom_cultures_pos.unox_IPL.noDGCC

# calculate for each lipid class, using appropriate standard curve

for (i in 1:nrow(Marchetti_diatom.conc_classes)) {
  
  # first, calculate pmol o.c.
  
  pmol.oc.thisclass =
    apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc[
      Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc$species==
        Marchetti_diatom.conc_classes$Lipid_class[i],14:20],c(1,2),splitpred,
      eval(parse(text = Marchetti_diatom.conc_classes$Model.low[i])),
      eval(parse(text = Marchetti_diatom.conc_classes$Model.hi[i])),
      Marchetti_diatom.conc_classes$Cutoff_PA[i])
  
  # store result as appropriate
  
  Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc[
    Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc$species==
      Marchetti_diatom.conc_classes$Lipid_class[i],14:20] = pmol.oc.thisclass
  
}
  
# now, scale pmol o.c. to pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_BD_Marchetti_diatoms_uL = 20 # amount DNPPE added per sample in uL, per VML B&D protocol

DNPPE_pmol_added_per_samp = DNPPE_mg_mL_BD_extracts_2016*(1/DNPPE_MW)*(10^9)*(1/10^3)*DNPPE_BD_Marchetti_diatoms_uL

Marchetti_diatom_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc[
  Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc$compound_name=="DNPPE",14:20]  # recovery factor

# create final results data frame
Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total = Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.oc
Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20] = sweep(as.matrix(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20]), 2, as.numeric(Marchetti_diatom_DNPPE.samp.RF), "*") # apply RF to samples, calculate total # pmol each species in given sample

# need to simplify the dataset a bit and do some QA

# eliminate features w/calibrated mass < 0 

Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20] =
  replace(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20],
          Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20]<0,
          NA)

# # determine what peak area cutoff we impose that still retains >95% of the ID'd mass in each sample
# 
# apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:18],2,function (x) sum(x[x>=30], na.rm = TRUE)/sum(x, na.rm = TRUE)) # 60 does it for cols 14-18
# 
# apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,19:20],2,function (x) sum(x[x>=4], na.rm = TRUE)/sum(x, na.rm = TRUE)) # 10 does it for cols 19-20
# 
# # apply these constraints
# 
# Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:18] =
#   replace(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:18],
#           Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:18]<30,
#           NA)
# 
# Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,19:20] =
#   replace(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,19:20],
#           Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,19:20]<4,
#           NA)

# now remove elements for which there is no data in any sample (includes elements
# just reduced to NA)
Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total = 
  Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20],1,sum,na.rm = T)>0,]

# can remove DNPPE at this point
Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total = 
  Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[!(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$compound_name=="DNPPE"),]

# now, finally, need to remove some duplicate features still apparently present

# Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total = Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[!(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$match_ID %in% 
#                                                                                                                         c(4910,964,4616,5409,852,285)),]

# some basic data analysis

# bar plots by species of lipid class distribution (molar basis)

# get list of IP DAG classes present
Marchetti_diatoms.IP_DAGclasses =
  unique(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species)
Marchetti_diatoms.IP_DAGclasses = Marchetti_diatoms.IP_DAGclasses[Marchetti_diatoms.IP_DAGclasses!=c("DNPPE")]

# preallocate matrix for results
Marchetti_diatoms.IP_DAGtotals = as.data.frame(matrix(NA,length(Marchetti_diatoms.IP_DAGclasses)+1,7))
rownames(Marchetti_diatoms.IP_DAGtotals) = c(Marchetti_diatoms.IP_DAGclasses,"Total")
colnames(Marchetti_diatoms.IP_DAGtotals) = colnames(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total)[14:20]

# calculate overall and class-specific totals
Marchetti_diatoms.IP_DAGtotals[c("Total"),]=
  apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[,14:20],2,sum,na.rm=T)

for (i in 1:(nrow(Marchetti_diatoms.IP_DAGtotals)-1)) {
  
  Marchetti_diatoms.IP_DAGtotals[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==rownames(Marchetti_diatoms.IP_DAGtotals)[i],14:20],2,sum,na.rm=T)
  
}

# put in rough descending order of abundance of first species listed 

Marchetti_diatoms.IP_DAGtotals = Marchetti_diatoms.IP_DAGtotals[order(Marchetti_diatoms.IP_DAGtotals[,1],decreasing = TRUE),]

# make a plot for thesis appendix A
# some of the samples are actually duplicates, so can eliminate some

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_IP-DAG_dist.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.IP_DAGtotals[2:8,c(1,3,5,7)])),margin=2)
barplot(prop, col=rainbow(length(rownames(prop))), width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401")
        )
legend("topright",inset=c(-0.25,0), fill=rainbow(length(rownames(prop))), density=c(70,60,50,35,30,20,15), legend=rownames(prop))

dev.off()

# an expansion of the 0.95-1 range of the y-axis

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_IP-DAG_dist_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.IP_DAGtotals[2:8,c(1,3,5,7)])),margin=2)
barplot(prop, col=rainbow(length(rownames(prop))), width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1)
)
legend("topright",inset=c(-0.25,0), fill=rainbow(length(rownames(prop))), density=c(70,60,50,35,30,20,15), legend=rownames(prop))

dev.off()

# bar plots by degree of unsaturation (rough categories, molar basis)

# without deliving into the frag spectra or running FAMES, can't tell what the
# distribution of 2x bonds looks like between the two acyl groups
# however, can make a few categories w/certainty assuming max unsaturation will be 
# 6 double bonds (based on known pathways of FA biosynthesis in diatoms)

# define some categories
Marchetti_diatoms.sat_classes = c("Both sn1, sn2 fully saturated","Neither sn1 nor sn2 is more than di-unsaturated",
                                  "Contain other species of middling unsaturation",
                                  "sn1, sn2 both have ≥ 3 double bonds","sn1, sn2 both have/PUFAs ≥ 5 DB")

# need also to define these in terms of numbers that will be used to extract data
Marchetti_diatoms.sat_numbers = c(0,2,9,11)

# preallocate matrix for results
Marchetti_diatoms.sat_totals = as.data.frame(matrix(NA,length(Marchetti_diatoms.sat_classes),7))
rownames(Marchetti_diatoms.sat_totals) = c(Marchetti_diatoms.sat_classes)
colnames(Marchetti_diatoms.sat_totals) = colnames(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total)[14:20]

# calculate saturation class totals

for (i in 1:(nrow(Marchetti_diatoms.sat_totals))) {
  
  if (i==1) {
    
    Marchetti_diatoms.sat_totals[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB==Marchetti_diatoms.sat_numbers[i],14:20],
                                             2,sum,na.rm=T)
    
  } else if (i > 1 & i < 5) {
    
    Marchetti_diatoms.sat_totals[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[(
      Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
        Marchetti_diatoms.sat_numbers[i-1] & 
        Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB<
        Marchetti_diatoms.sat_numbers[i]),14:20],
      2,sum,na.rm=T)
    
  } else if (i==5) {
    
    Marchetti_diatoms.sat_totals[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=Marchetti_diatoms.sat_numbers[i-1],14:20],
                                             2,sum,na.rm=T)
    
  }
  
}

# make a plot
# some of the samples are actually duplicates, so can eliminate some

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# expansion of y-axis

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()


# maybe a plot like this of just PC

# preallocate matrix for results
Marchetti_diatoms.sat_totals.PC = as.data.frame(matrix(NA,length(Marchetti_diatoms.sat_classes),7))
rownames(Marchetti_diatoms.sat_totals.PC) = c(Marchetti_diatoms.sat_classes)
colnames(Marchetti_diatoms.sat_totals.PC) = colnames(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total)[14:20]

# calculate saturation class totals just for PC

for (i in 1:(nrow(Marchetti_diatoms.sat_totals.PC))) {
  
  if (i==1) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB==
         Marchetti_diatoms.sat_numbers[i] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  } else if (i > 1 & i < 5) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
         Marchetti_diatoms.sat_numbers[i-1] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB<
         Marchetti_diatoms.sat_numbers[i] &
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  } else if (i==5) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
         Marchetti_diatoms.sat_numbers[i-1] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  }
  
}

# make a plot
# some of the samples are actually duplicates, so can eliminate some

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_PConly.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals.PC[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# y-axis expansion

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_PConly_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals.PC[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# additional calculations

# % peak area accounted for by DGCC (didn't have a DGCC standard at time of analysis)

Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE = Marchetti_diatom_cultures_pos.unox_IPL[
  Marchetti_diatom_cultures_pos.unox_IPL$species!="DNPPE",]

Marchetti_cultures.peaksums = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE[,14:20],2,sum,na.rm=T)
Marchetti_cultures.DGCC.peaksums = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE[Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE$species=="DGCC",14:20],2,sum,na.rm=T)

Marchetti_cultures.DGCC.peaksums/Marchetti_cultures.peaksums

### Particulate environmental samples from PAL1314, LMG1401 ####

## pull in data ###

# dataset includes putative IDs of species w/odd-chain FA's

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/PAL1314_LMG1401_particulate_enviro_samples_pos_withoddFA.RData")
PAL1314_LMG1401_partic_pos = getLOBpeaklist(PAL1314_LMG1401_particulate_samples_pos_withoddFA) # generate peaklist

# load("data/nice/Orbi_MS_data/LOBSTAHS_processed/PAL1314_LMG1401_particulate_enviro_samples_nooddFA.RData")
# PAL1314_LMG1401_partic_pos = getLOBpeaklist(PAL1314_LMG1401_particulate_enviro_samples_pos) # generate peaklist

# will use the 20161107 standards, with the 20161005 for DGTS 

# extract only unoxidized IPL (no TAGs, etc) data, plus DNPPE

PAL1314_LMG1401_partic_pos.unox_IPL = PAL1314_LMG1401_partic_pos[
  (PAL1314_LMG1401_partic_pos$lipid_class %in% c("IP_DAG","DNPPE") & 
     PAL1314_LMG1401_partic_pos$degree_oxidation==0),]

# # convert anything < 1e5 intensity to NA, assuming it's either noise
# # or something so low in concentraton as to be irrelevant from a total lipid
# # perspective
# 
# PAL1314_LMG1401_partic_pos.unox_IPL[,13:17] =
#   replace(PAL1314_LMG1401_partic_pos.unox_IPL[,13:17],
#           PAL1314_LMG1401_partic_pos.unox_IPL[,13:17]<10000,
#           NA)

# get rid of DGCC data since we don't have any standards for these right now

PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC = PAL1314_LMG1401_partic_pos.unox_IPL[
  PAL1314_LMG1401_partic_pos.unox_IPL$species!="DGCC",]

# define classes for which concentrations are to be calculated, the models, and
# the cutoffs in case of split prediction

PAL1314_LMG1401_partic.conc_classes = as.data.frame(matrix(NA,8,4))
colnames(PAL1314_LMG1401_partic.conc_classes) =
  c("Lipid_class","Model.low","Model.hi","Cutoff_PA")

PAL1314_LMG1401_partic.conc_classes[,1] =
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGTS_DGTA","DNPPE")

PAL1314_LMG1401_partic.conc_classes[,2] =
  c("linfit_low.PG.20161107","linfit_low.PE.20161107","linfit_low.PC.20161107",
    "linfit_low.MGDG.20161107","linfit_low.SQDG.20161107",
    "linfit_low.DGDG.20161107","linfit_low.DGTS_DGTA.20161005",
    "linfit_low.DNPPE.20161107")

PAL1314_LMG1401_partic.conc_classes[,3] =
  c("linfit_hi.PG.20161107","linfit_hi.PE.20161107","linfit_hi.PC.20161107",
    "linfit_hi.MGDG.20161107","linfit_hi.SQDG.20161107",
    "linfit_hi.DGDG.20161107","linfit_hi.DGTS_DGTA.20161005",
    "linfit_hi.DNPPE.20161107")

PAL1314_LMG1401_partic.conc_classes[,4] =
  c(PG_std_breakpoint.20161107,PE_std_breakpoint.20161107,PC_std_breakpoint.20161107,
    MGDG_std_breakpoint.20161107,SQDG_std_breakpoint.20161107,DGDG_std_breakpoint.20161107,
    DGTS_DGTA_std_breakpoint.20161005,DNPPE_std_breakpoint.20161107)

# PAL1314_LMG1401_partic.conc_classes[,2] =
#   c("linfit_low.PG.20161005","linfit_low.PE.20161005","linfit_low.PC.20161005",
#     "linfit_low.MGDG.20161005","linfit_low.SQDG.20161005",
#     "linfit_low.DGDG.20161005","linfit_low.DGTS_DGTA.20161005",
#     "linfit_low.DNPPE.20161005")
# 
# PAL1314_LMG1401_partic.conc_classes[,3] =
#   c("linfit_hi.PG.20161005","linfit_hi.PE.20161005","linfit_hi.PC.20161005",
#     "linfit_hi.MGDG.20161005","linfit_hi.SQDG.20161005",
#     "linfit_hi.DGDG.20161005","linfit_hi.DGTS_DGTA.20161005",
#     "linfit_hi.DNPPE.20161005")
# 
# PAL1314_LMG1401_partic.conc_classes[,4] =
#   c(PG_std_breakpoint.20161005,PE_std_breakpoint.20161005,PC_std_breakpoint.20161005,
#     MGDG_std_breakpoint.20161005,SQDG_std_breakpoint.20161005,DGDG_std_breakpoint.20161005,
#     DGTS_DGTA_std_breakpoint.20161005,DNPPE_std_breakpoint.20161005)

# first, calculate pmol o.c.
# preallocate matrix for result
PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc = PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC

# calculate for each lipid class, using appropriate standard curve

for (i in 1:nrow(PAL1314_LMG1401_partic.conc_classes)) {
  
  # first, calculate pmol o.c.
  
  pmol.oc.thisclass =
    apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc[
      PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc$species==
        PAL1314_LMG1401_partic.conc_classes$Lipid_class[i],13:17],c(1,2),splitpred,
      eval(parse(text = PAL1314_LMG1401_partic.conc_classes$Model.low[i])),
      eval(parse(text = PAL1314_LMG1401_partic.conc_classes$Model.hi[i])),
      PAL1314_LMG1401_partic.conc_classes$Cutoff_PA[i])
  
  # store result as appropriate
  
  PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc[
    PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc$species==
      PAL1314_LMG1401_partic.conc_classes$Lipid_class[i],13:17] = pmol.oc.thisclass
  
}

# now, scale pmol o.c. to pmol per sample using DNPPE (recovery standard added at time of extraction)

DNPPE_BD_Marchetti_diatoms_uL = 20 # amount DNPPE added per sample in uL, per VML B&D protocol
# same amount as for the Marchetti diatom cultures

DNPPE_pmol_added_per_samp = DNPPE_mg_mL_BD_extracts_2016*(1/DNPPE_MW)*(10^9)*(1/10^3)*DNPPE_BD_Marchetti_diatoms_uL

PAL1314_LMG1401_partic.samp.RF = DNPPE_pmol_added_per_samp/PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc[
  PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc$compound_name=="DNPPE",13:17]  # recovery factor

# create final results data frame
PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total = PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.oc
PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17] = sweep(as.matrix(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17]), 2, as.numeric(PAL1314_LMG1401_partic.samp.RF), "*") # apply RF to samples, calculate total # pmol each species in given sample

# need to simplify the dataset a bit and do some QA

# eliminate features w/calibrated mass < 0 

PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17] =
  replace(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17],
          PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17]<0,
          NA)

# # determine what peak area cutoff we impose that still retains >95% of the ID'd mass in each sample
# 
# apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17],2,function (x) sum(x[x>=30], na.rm = TRUE)/sum(x, na.rm = TRUE)) # 30 does it
# 
# # apply these constraints
# 
# PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17] =
#   replace(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17],
#           PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17]<30,
#           NA)

# now remove elements for which there is no data in any sample (includes elements
# just reduced to NA)
PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total = 
  PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17],1,sum,na.rm = T)>0,]

# can remove DNPPE at this point
PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total = 
  PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[!(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$compound_name=="DNPPE"),]

# now, finally, need to remove some duplicate features still apparently present

# Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total = Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[!(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$match_ID %in% 
#                                                                                                                         c(4910,964,4616,5409,852,285)),]

# some basic data analysis

# bar plots by species of lipid class distribution (molar basis)

# get list of IP DAG classes present
PAL1314_LMG1401_partic.IP_DAGclasses =
  unique(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$species)
PAL1314_LMG1401_partic.IP_DAGclasses = PAL1314_LMG1401_partic.IP_DAGclasses[PAL1314_LMG1401_partic.IP_DAGclasses!=c("DNPPE")]

# preallocate matrix for results
PAL1314_LMG1401_partic.IP_DAGtotals = as.data.frame(matrix(NA,length(PAL1314_LMG1401_partic.IP_DAGclasses)+1,5))
rownames(PAL1314_LMG1401_partic.IP_DAGtotals) = c(PAL1314_LMG1401_partic.IP_DAGclasses,"Total")
colnames(PAL1314_LMG1401_partic.IP_DAGtotals) = colnames(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total)[13:17]

# calculate overall and class-specific totals
PAL1314_LMG1401_partic.IP_DAGtotals[c("Total"),]=
  apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[,13:17],2,sum,na.rm=T)

for (i in 1:(nrow(PAL1314_LMG1401_partic.IP_DAGtotals)-1)) {
  
  PAL1314_LMG1401_partic.IP_DAGtotals[i,] = apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$species==rownames(PAL1314_LMG1401_partic.IP_DAGtotals)[i],13:17],2,sum,na.rm=T)
  
}

# put in rough descending order of abundance of first species listed 

PAL1314_LMG1401_partic.IP_DAGtotals = PAL1314_LMG1401_partic.IP_DAGtotals[order(PAL1314_LMG1401_partic.IP_DAGtotals[,1],decreasing = TRUE),]

# make a plot for thesis

# par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
# 
# pdf(file = "PAL1314_LMG1401_partic_IP-DAG_dist.pdf",
#     width = 8, height = 6, pointsize = 12,
#     bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(PAL1314_LMG1401_partic.IP_DAGtotals[2:8,])),margin=2)
barplot(prop, col=rainbow(length(rownames(prop))), width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = colnames(PAL1314_LMG1401_partic.IP_DAGtotals)
)
legend("topright",inset=c(-0.25,0), fill=rainbow(length(rownames(prop))), density=c(70,60,50,35,30,20,15), legend=rownames(prop))

dev.off()

# an expansion of the 0.95-1 range of the y-axis

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_IP-DAG_dist_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.IP_DAGtotals[2:8,c(1,3,5,7)])),margin=2)
barplot(prop, col=rainbow(length(rownames(prop))), width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1)
)
legend("topright",inset=c(-0.25,0), fill=rainbow(length(rownames(prop))), density=c(70,60,50,35,30,20,15), legend=rownames(prop))

dev.off()

# bar plots by degree of unsaturation (rough categories, molar basis)

# without deliving into the frag spectra or running FAMES, can't tell what the
# distribution of 2x bonds looks like between the two acyl groups
# however, can make a few categories w/certainty assuming max unsaturation will be 
# 6 double bonds (based on known pathways of FA biosynthesis in diatoms)

# define some categories
PAL1314_LMG1401_partic.sat_classes = c("Both sn1, sn2 fully saturated","Neither sn1 nor sn2 is more than di-unsaturated",
                                  "Contain other species of middling unsaturation",
                                  "sn1, sn2 both have ≥ 3 double bonds","sn1, sn2 both have/PUFAs ≥ 5 DB")

# need also to define these in terms of numbers that will be used to extract data
PAL1314_LMG1401_partic.sat_numbers = c(0,2,9,11)

# preallocate matrix for results
PAL1314_LMG1401_partic.sat_totals = as.data.frame(matrix(NA,length(PAL1314_LMG1401_partic.sat_classes),5))
rownames(PAL1314_LMG1401_partic.sat_totals) = c(PAL1314_LMG1401_partic.sat_classes)
colnames(PAL1314_LMG1401_partic.sat_totals) = colnames(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total)[13:17]

# calculate saturation class totals

for (i in 1:(nrow(PAL1314_LMG1401_partic.sat_totals))) {
  
  if (i==1) {
    
    PAL1314_LMG1401_partic.sat_totals[i,] = apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB==PAL1314_LMG1401_partic.sat_numbers[i],13:17],
                                             2,sum,na.rm=T)
    
  } else if (i > 1 & i < 5) {
    
    PAL1314_LMG1401_partic.sat_totals[i,] = apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[(
      PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
        PAL1314_LMG1401_partic.sat_numbers[i-1] & 
        PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB<
        PAL1314_LMG1401_partic.sat_numbers[i]),13:17],
      2,sum,na.rm=T)
    
  } else if (i==5) {
    
    PAL1314_LMG1401_partic.sat_totals[i,] = apply(PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total[PAL1314_LMG1401_partic_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=PAL1314_LMG1401_partic.sat_numbers[i-1],13:17],
                                             2,sum,na.rm=T)
    
  }
  
}

# make a plot
# some of the samples are actually duplicates, so can eliminate some

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PAL1314_LMG1401_partic_satur_dist.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(PAL1314_LMG1401_partic.sat_totals[1:5,])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = colnames(PAL1314_LMG1401_partic.sat_totals))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# expansion of y-axis

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()


# maybe a plot like this of just PC

# preallocate matrix for results
Marchetti_diatoms.sat_totals.PC = as.data.frame(matrix(NA,length(Marchetti_diatoms.sat_classes),7))
rownames(Marchetti_diatoms.sat_totals.PC) = c(Marchetti_diatoms.sat_classes)
colnames(Marchetti_diatoms.sat_totals.PC) = colnames(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total)[14:20]

# calculate saturation class totals just for PC

for (i in 1:(nrow(Marchetti_diatoms.sat_totals.PC))) {
  
  if (i==1) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB==
         Marchetti_diatoms.sat_numbers[i] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  } else if (i > 1 & i < 5) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
         Marchetti_diatoms.sat_numbers[i-1] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB<
         Marchetti_diatoms.sat_numbers[i] &
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  } else if (i==5) {
    
    Marchetti_diatoms.sat_totals.PC[i,] = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total[
      (Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$FA_total_no_DB>=
         Marchetti_diatoms.sat_numbers[i-1] & 
         Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total$species==
         "PC"),
      14:20],2,sum,na.rm=T)
    
  }
  
}

# make a plot
# some of the samples are actually duplicates, so can eliminate some

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_PConly.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(8, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals.PC[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# y-axis expansion

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Marchetti_diatom_satur_dist_PConly_expansion.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(16, 4.1, 4.1, 7.1), xpd=TRUE)
prop = prop.table(as.table(as.matrix(Marchetti_diatoms.sat_totals.PC[1:5,c(1,3,5,7)])),margin=2)
barplot(prop, col=cm.colors((length(rownames(prop))+1))[c(1:4,6)], width=2, density=c(70,60,50,35,30,20,15),las=2,
        ylab = "Relative molar abundance",
        xlab = "Species",
        names.arg = c("Actinocyclus	actinochilus UNC 1403","Chaetoceros sp. UNC 1408",
                      "Fragilariopsis cylindrus UNC1301","Thalassiosira antarctica UNC1401"),
        ylim = c(0.95,1))
legend("bottomright",inset=c(-.25,-0.4), fill=cm.colors((length(rownames(prop))+1))[c(1:4,6)], density=c(70,60,50,35,30,20,15), legend=rownames(prop), cex=.5)

dev.off()

# additional calculations

# % peak area accounted for by DGCC (didn't have a DGCC standard at time of analysis)

Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE = Marchetti_diatom_cultures_pos.unox_IPL[
  Marchetti_diatom_cultures_pos.unox_IPL$species!="DNPPE",]

Marchetti_cultures.peaksums = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE[,14:20],2,sum,na.rm=T)
Marchetti_cultures.DGCC.peaksums = apply(Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE[Marchetti_diatom_cultures_pos.unox_IPL.noDNPPE$species=="DGCC",14:20],2,sum,na.rm=T)

Marchetti_cultures.DGCC.peaksums/Marchetti_cultures.peaksums

