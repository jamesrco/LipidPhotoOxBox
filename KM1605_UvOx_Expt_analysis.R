# KM1605_UvOx_Expt_analysis.R

# Analysis of an ultraviolet radiation oxiation experiment conducted with large-volume whole seawater samples during a SCOPE cruise (KM1605) aboard the R/V Kilo Moana

# Created 11/29/16 by J.R.C.

# Follows same general initial path as PAL1314_PAL1516_environmental_samples.R

### workspace preparation ####

# load necessary libraries

library(LOBSTAHS)
library(chemCal)
library(stats)
library(stringr)

# set wd

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")

### load some necessary metadata; define functions ###

# some chemical data

DNPPE_mg_mL_BD_extracts_2016 = 0.051 # concentration of the DNPPE added during B&D extractions
# in Van Mooy Lab in 2016
DNPPE_MW = 875.081 # MW DNPPE, g/mol

### functions for general lookup and computation ###

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

# define function for TAG concentration calculation

TAGconcCalc= function(x,numC,numDB,TAGfitFunction,linDNPPEfit_low,linDNPPEfit_hi,DNPPEfit.cutoff) {
  
  # calculate ECN (equivalent carbon number)
  
  ECN = numC-2*numDB
  
  # create data frame for lookup
  
  ECN.lookup = as.data.frame(ECN)
  colnames(ECN.lookup) = "x"
  
  # retrieve RRF
  
  this.RRF = predict(TAGfitFunction, newdata = ECN.lookup)
  
  # calculate RRF-adjusted peak area
  
  adj.peakarea = this.RRF*x
  
  # calculate pmol o.c. of TAG (by way of DNPPE standard curve, since the TAG RRFs were calculated relative to DNPPE)
  
  TAGpmol.oc = splitpred(adj.peakarea,
            linfit_low = linDNPPEfit_low,
            linfit_hi = linDNPPEfit_hi,
            cutoff = DNPPEfit.cutoff)
  
  return(TAGpmol.oc)
  
}

# define function allow "easy" retrieval of sample metadata based on sample ID in filename

getMetDat = function(fn,metadat.raw,whichdat) {
  
  # fn is sample descriptor (from verbose file name) from which sample ID will be extracted and then matched
  # whichdat is a list of metadata to be retrieved (corresponding to name of column in the metadat.raw table)
  
  match_ind = grepl(sub("^.*_","\\1",fn),metadat.raw$Orbi.sequence.ID) # get index of matching metadata
  
  this.metdat = metadat.raw[match_ind,whichdat] # extract called-for metadata
  
  this.metdat # return extracted data
  
}

### functions for evaluation of fragmentation spectra ###

# necessary functions that will allow us to extract the correct ms2 spectra, evaluate transitions, etc

get.ms2Peaklist = function (precursor.index,sample_ID,xcmsRaw.list) {
  
  ms2data.start = xcmsRaw.list[[sample_ID]]@msnScanindex[precursor.index]
  ms2data.end = xcmsRaw.list[[sample_ID]]@msnScanindex[precursor.index+1]-1
  
  scandata =
    data.frame(xcmsRaw.list[[sample_ID]]@env$msnMz[ms2data.start:ms2data.end],
               xcmsRaw.list[[sample_ID]]@env$msnIntensity[ms2data.start:ms2data.end])
  
  colnames(scandata) = c("mz","Intensity")
  
  return(scandata)
  
}

get.topN = function(peaklist,N) {
  
  ordered.peaklist = peaklist[order(peaklist$Intensity, decreasing = TRUE),]
  
  topN.peaklist = ordered.peaklist[1:N,]
  
  return(topN.peaklist)
  
}

eval.PIspecies = function(peaklist,species,ppm,ms2.lookupClasses) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # retrieve product ion m/z
  
  prod.ion = ms2.lookupClasses$mz_value[
    rownames(ms2.lookupClasses)==species]
  
  if (any(abs((prod.ion-peaklist[,1])/prod.ion*1000000)<ppm)) {
    
    # it's a match
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}

eval.CNLspecies = function(peaklist,species,sample_ID,ppm,ms2.lookupClasses) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # calculate theoretical m/z of the ion that would be produced via the neutral loss
  # throwing in a mean() here in case xcms associated more than one peak with the group in this particular sample
  CNL_product_mz = mean(xcms.peakdata_thisgroup_pos[xcms.peakdata_thisgroup_pos$sample==sample_ID,1])-
    ms2.lookupClasses$mz_value[
      rownames(ms2.lookupClasses)==this.IDclass]
  
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

#### standards ####

### + mode standards  ###

### standards from 20161107  ####

# QE003063-QE003073

# standards from 20161107, for PAL1314 & LMG1401 particulate data
# these didn't include any betaine standards

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/6IPL_Standards_20161107_pos.RData") # load processed standard data
IPLstd_pos_20161107.raw = getLOBpeaklist(x6IPL_Standards_20161107_pos) # generate peaklist

# extract standards for each species, and DNPPE
# **** need to extract more than one species for DGDG, SQDG, MGDG since these were not pure standards of one species, and some of the "secondary" (less abundant) species are still abundant enough that they contribute to the overall mass of the lipid

Std_peakareas.20161107 = IPLstd_pos_20161107.raw[
  IPLstd_pos_20161107.raw$compound_name %in% c("PG 32:0","PE 32:0","PC 32:0",
                                               "MGDG 36:0","MGDG 34:0","SQDG 34:3",
                                               "SQDG 34:2","SQDG 36:6","SQDG 32:0",
                                               "SQDG 32:3","SQDG 36:5","SQDG 34:1",
                                               "SQDG 36:4",
                                               "DGDG 36:4","DGDG 34:2","DGDG 36:3",
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
linfit_low.PG.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 9 standard levels
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
linfit_low.PE.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 9 standard levels
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
linfit_low.PC.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 10 standard levels
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

MGDG.stdareas.20161107 = apply(Std_peakareas.20161107[grep("^MGDG",rownames(Std_peakareas.20161107)),],2,sum,na.rm = T)

y = MGDG.stdareas.20161107[1:9]
x = rev(Stds_20161107_oc$pmol_oc_MGDG)[1:9]
linfit_low.MGDG.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 10 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_MGDG),
     MGDG.stdareas.20161107,
     pch="+",
     ylab = "Peak area, total MGDG",
     xlab = "pmol o.c., total MGDG")
points(rev(Stds_20161107_oc$pmol_oc_MGDG)[1:9],fitted(linfit_low.MGDG.20161107),col="red",pch="+")

# we will need some other fit for levels higher than ~ 155 pmol o.c.

y = MGDG.stdareas.20161107[10:11]
x = rev(Stds_20161107_oc$pmol_oc_MGDG)[10:11]
linfit_hi.MGDG.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_MGDG)[10:11],fitted(linfit_hi.MGDG.20161107),col="blue",pch="+")

MGDG_std_breakpoint.20161107 = MGDG.stdareas.20161107[10]

# SQDG

# curve fitting & diagnostics

SQDG.stdareas.20161107 = apply(Std_peakareas.20161107[grep("^SQDG",rownames(Std_peakareas.20161107)),],2,sum,na.rm = T)

y = SQDG.stdareas.20161107[1:9]
x = rev(Stds_20161107_oc$pmol_oc_SQDG)[1:9]
linfit_low.SQDG.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_SQDG),
     SQDG.stdareas.20161107,
     pch="+",
     ylab = "Peak area, total SQDG",
     xlab = "pmol o.c., total SQDG")
points(rev(Stds_20161107_oc$pmol_oc_SQDG)[1:9],fitted(linfit_low.SQDG.20161107),col="red",pch="+")

# the second-highest standard appears to be messed up, will need to develop a
# different approach for response of highest peak areas

y = SQDG.stdareas.20161107[c(9,11)]
x = rev(Stds_20161107_oc$pmol_oc_SQDG)[c(9,11)]
linfit_hi.SQDG.20161107 = lm(as.numeric(y)~x)
points(rev(Stds_20161107_oc$pmol_oc_SQDG)[c(9,11)],fitted(linfit_hi.SQDG.20161107),col="blue",pch="+")

SQDG_std_breakpoint.20161107 = SQDG.stdareas.20161107[9]

# DGDG

DGDG.stdareas.20161107 = apply(Std_peakareas.20161107[grep("^DGDG",rownames(Std_peakareas.20161107)),],2,sum,na.rm = T)

# curve fitting & diagnostics
# just need one curve here, as long as we omit the second-highest level

y = DGDG.stdareas.20161107[c(1:9,11)]
x = rev(Stds_20161107_oc$pmol_oc_DGDG)[c(1:9,11)]
linfit_low.DGDG.20161107 = lm(as.numeric(y)~x-1) # fit a linear model for the first 9 standard levels
plot(rev(Stds_20161107_oc$pmol_oc_DGDG),
     DGDG.stdareas.20161107,
     pch="+",
     ylab = "Peak area, DGDG 36:4",
     xlab = "pmol o.c., DGDG 36:4")
points(rev(Stds_20161107_oc$pmol_oc_DGDG)[c(1:9,11)],fitted(linfit_low.DGDG.20161107),col="red",pch="+")

# can define the high-range curve to be the same as the low-range curve, in this case

linfit_hi.DGDG.20161107 = linfit_low.DGDG.20161107
DGDG_std_breakpoint.20161107 = DGDG.stdareas.20161107[9]

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
linfit_low.DNPPE.20161107 = lm(as.numeric(y)~x-1) # fit a linear model to first 8 points
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

# **** need to extract more than one species for DGDG, SQDG, MGDG since these were not pure standards of one species, and some of the "secondary" (less abundant) species are still abundant enough that they contribute to the overall mass of the lipid

Std_peakareas.20161005 = IPLstd_pos_20161005.raw[
  IPLstd_pos_20161005.raw$compound_name %in% c("PG 32:0","PE 32:0","PC 32:0",
                                               "MGDG 36:0","MGDG 34:0","SQDG 34:3",
                                               "SQDG 34:2","SQDG 36:6","SQDG 32:0",
                                               "SQDG 32:3","SQDG 36:5","SQDG 34:1",
                                               "SQDG 36:4",
                                               "DGDG 36:4","DGDG 34:2","DGDG 36:3",
                                               "DNPPE",
                                               "DGTS_DGTA 32:0","DGTS_DGTA 34:0"),]

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

MGDG.stdareas.20161005 = apply(Std_peakareas.20161005[grep("^MGDG",rownames(Std_peakareas.20161005)),],2,sum,na.rm = T)

# curve fitting & diagnostics

y = MGDG.stdareas.20161005[1:8]
x = rev(Stds_20161005_oc$pmol_oc_MGDG)[1:8]
linfit_low.MGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 10 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_MGDG),
     MGDG.stdareas.20161005,
     pch="+",
     ylab = "Peak area, total MGDG",
     xlab = "pmol o.c., total MGDG")
points(rev(Stds_20161005_oc$pmol_oc_MGDG)[1:8],fitted(linfit_low.MGDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = MGDG.stdareas.20161005[8:10]
x = rev(Stds_20161005_oc$pmol_oc_MGDG)[8:10]
linfit_hi.MGDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_MGDG)[8:10],fitted(linfit_hi.MGDG.20161005),col="blue",pch="+")

MGDG_std_breakpoint.20161005 = MGDG.stdareas.20161005[8]

# SQDG

SQDG.stdareas.20161005 = apply(Std_peakareas.20161005[grep("^SQDG",rownames(Std_peakareas.20161005)),],2,sum,na.rm = T)

# curve fitting & diagnostics
# something appears to be very weird with the SQDG in these standards
# will generate one curve while ommitting the 8th and 9th points

y = SQDG.stdareas.20161005[1:7]
# x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7,10)]
x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7)]
linfit_low.SQDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 7 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_SQDG),
     SQDG.stdareas.20161005,
     pch="+",
     ylab = "Peak area, total SQDG",
     xlab = "pmol o.c., total SQDG")
points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(1:7)],fitted(linfit_low.SQDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = SQDG.stdareas.20161005[c(7,10)]
x = rev(Stds_20161005_oc$pmol_oc_SQDG)[c(7,10)]
linfit_hi.SQDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_SQDG)[c(7,10)],fitted(linfit_hi.SQDG.20161005),col="blue",pch="+")

SQDG_std_breakpoint.20161005 = SQDG.stdareas.20161005[7]

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

DGDG.stdareas.20161005 = apply(Std_peakareas.20161005[grep("^DGDG",rownames(Std_peakareas.20161005)),],2,sum,na.rm = T)

# curve fitting & diagnostics

y = DGDG.stdareas.20161005[1:8]
x = rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8)]
linfit_low.DGDG.20161005 = lm(as.numeric(y)~x) # fit a linear model for the first 7 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGDG),
     DGDG.stdareas.20161005,
     pch="+",
     ylab = "Peak area, total DGDG",
     xlab = "pmol o.c., total DGDG")
points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(1:8)],fitted(linfit_low.DGDG.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = DGDG.stdareas.20161005[c(8,10)]
x = rev(Stds_20161005_oc$pmol_oc_DGDG)[c(8,10)]
linfit_hi.DGDG.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_DGDG)[c(8,10)],fitted(linfit_hi.DGDG.20161005),col="blue",pch="+")

DGDG_std_breakpoint.20161005 = DGDG.stdareas.20161005[8]

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

DGTS_DGTA.stdareas.20161005 = apply(Std_peakareas.20161005[grep("^DGTS_DGTA",rownames(Std_peakareas.20161005)),],2,sum,na.rm = T)

# curve fitting & diagnostics

y = DGTS_DGTA.stdareas.20161005[1:8]
x = rev(Stds_20161005_oc$pmol_oc_DGTS)[c(1:8)]

linfit_low.DGTS_DGTA.20161005 = lm(as.numeric(y)~x) # fit a model for the first 8 standard levels
plot(rev(Stds_20161005_oc$pmol_oc_DGTS),
     DGTS_DGTA.stdareas.20161005,
     pch="+",
     ylab = "Peak area, DGTS_DGTA 32:0 + 34:0",
     xlab = "pmol o.c., DGTS_DGTA 32:0 + 34:0")
points(rev(Stds_20161005_oc$pmol_oc_DGTS)[c(1:8)],fitted(linfit_low.DGTS_DGTA.20161005),col="red",pch="+")

# we will need some other fit for levels higher than ~ 40 pmol o.c.

y = DGTS_DGTA.stdareas.20161005[8:9]
x = rev(Stds_20161005_oc$pmol_oc_DGTS)[8:9]
linfit_hi.DGTS_DGTA.20161005 = lm(as.numeric(y)~x)
points(rev(Stds_20161005_oc$pmol_oc_DGTS)[8:9],fitted(linfit_hi.DGTS_DGTA.20161005),col="blue",pch="+")

DGTS_DGTA_std_breakpoint.20161005 = DGTS_DGTA.stdareas.20161005[8]

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

# maybe a plot of the reciprocals

plot(ECNs,1/RRFs[2:length(RRFs)])
text(ECNs, 1/RRFs[2:length(RRFs)], labels=names(ECNs), pos = 4, cex = 0.5)

# x = ECNs
# y = 1/RRFs[2:length(RRFs)]

# excluding 30:0, 33:0

x = ECNs[3:length(ECNs)]
y = 1/RRFs[4:length(RRFs)]

# fit a non-linear model of form
# y = a*b^x-c

nls.fit.TAG20161004 = nls(y~a*b^x-c,list(x,y),c(a=0.5,b=1.25,c=30), nls.control(maxiter=5000),
                          lower = c(0.000001,0.5,-1000), algorithm = "port",
                          upper = c(6,3,500))
points(x,fitted(nls.fit.TAG20161004),col="red",pch="+")
text(ECNs[1:2],1/RRFs[2:3],c("EXCLUDED","EXCLUDED"), cex = 0.5, pos = 1, col = "red")

### pull in data ####

# will use the 20161107 standards for these data, with the 20161005 standards for DGTS & DGTA

load("data/nice/Orbi_MS_data/LOBSTAHS_processed/KM1605_UvOX_Expt_pos_withoddFA_LOBSet.RData")
KM1605_UvOX_Expt_pos_withoddFA_particulate = getLOBpeaklist(KM1605_UvOX_Expt_pos_withoddFA_LOBSet) # generate peaklist

# # get rid of the two QC's
# KM1605_UvOX_Expt_pos_withoddFA_particulate = KM1605_UvOX_Expt_pos_withoddFA_particulate[,-c(21:22)]

### ***** experimental trial of use of msn data using series of xcmsRaw objects ***** ####

# now, use fragmentation spectra for basic confirmation of putative LOBSTAHS IDs

# requires several additional files: (1) annotated xcmsSet object for data in positive ion mode, (2) LOBSet object in positive ion mode, and (3) list containing positive-mode xcmsRaw objects for all samples generated using: 
#             xraw <- xcmsRaw("yourfile.mzXML", includeMSn=TRUE)
# the routine below assumes the xcmsRaw objects are stored in the list in the order in which they appear in the xsAnnotate and LOBSet objects

# first, define types and values of MS fragmentation experiments for each IP-DAG class
# (i.e., constant neutral loss (CNL) or product ion (PI))

# create empty data frame
KM1605_UvOX_part.frag_lookup_classes = as.data.frame(matrix(NA,8,2))
rownames(KM1605_UvOX_part.frag_lookup_classes) = 
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGCC","DGTS_DGTA")
colnames(KM1605_UvOX_part.frag_lookup_classes) =
  c("Frag_exp_type","mz_value")

# now, populate our data frame with necessary values
# per Popendorf et al., Lipids (2013) 48:185–195

KM1605_UvOX_part.frag_lookup_classes[,1] =
  c("CNL","CNL","PI","CNL","CNL","CNL","PI","PI")
KM1605_UvOX_part.frag_lookup_classes[,2] =
  c(189.040224,141.019094,184.073321,197.089937,261.051837,359.142212,104.106690,236.149249)

# pull in the xsAnnotate object & list containing the xcmsRaw objects

load("data/nice/Orbi_MS_data/xsAnnotate_objects/KM1605_UvOX_Expt_pos_withoddFA_xsAnnotate.RData")
KM1605_UvOX_p_xsA_pos = KM1605_UvOX_Expt_pos_withoddFA_xsAnnotate # rename so easier to work with

load("data/raw/Orbi_MS_data/xcmsRaw_objects/KM1605_UvOX_Expt_pos_withoddFA_xcmsRaw.RData")
KM1605_UvOX_p_xraw_pos = KM1605_UvOX_Expt_pos_withoddFA_xcmsRaw

# preallocate two matrices for our results
# KM1605_UvOX_p.detected_pos_ion_fragments: to keep track of how many valid ms2 fragmentation spectra were detected for the feature in positive ion mode
# KM1605_UvOX_p.fragdata_results: number of ms2 spectra in which the PI or CNL criteria were validated

KM1605_UvOX_p.detected_pos_ion_fragments = as.data.frame(matrix(NA,nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate),ncol=13))
KM1605_UvOX_p.fragdata_results = as.data.frame(matrix(NA,nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate),ncol=13))

# iterate through the IP-DAG IDs by sample, retrieve necessary data from the xsAnnotate and xcmsRaw objects, evaluate, and record the results

for (i in 1:(nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate))) {
  # iterate through each LOBSTAHS ID
  
  # retrieve LOBSTAHS compound ID, lipid species
  
  this.ID = KM1605_UvOX_Expt_pos_withoddFA_particulate$compound_name[i]
  this.IDclass = KM1605_UvOX_Expt_pos_withoddFA_particulate$species[i]
  
  if (this.IDclass %in% rownames(KM1605_UvOX_part.frag_lookup_classes)) {
    # an escape if the lipid class isn't accounted for in our input parameter table
    
    # retrieve underlying positive-mode xcms group and peak data
    
    xcms.peakIDs_thisgroup_pos =
      KM1605_UvOX_p_xsA_pos@xcmsSet@groupidx[[KM1605_UvOX_Expt_pos_withoddFA_particulate$xcms_peakgroup[i]]]
    
    xcms.peakdata_thisgroup_pos =
      as.data.frame(KM1605_UvOX_p_xsA_pos@xcmsSet@peaks[xcms.peakIDs_thisgroup_pos,])
    
    # now, iterate through all instances of this putatively identified compound that are present in the xcms peakgroup and in one of the samples of interest (i.e., not a QC) 
    
    for (j in 1:nrow(xcms.peakdata_thisgroup_pos)) {
      
      # retrieve the sample ID
      
      samp_ID = xcms.peakdata_thisgroup_pos$sample[j]
      
      # before beginning, check if the peak is a QC; if so, skip to end

      if (!(samp_ID %in% c(6:7))) {
        
        ### first, pull out fragmentation spectra relevant to this peak, if they exist ###
        
        # retrieve raw (uncorrected) retention times for this peak
        
        RT_min_pos.raw = KM1605_UvOX_p_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          KM1605_UvOX_p_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmin[j]]
        
        RT_max_pos.raw = KM1605_UvOX_p_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          KM1605_UvOX_p_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmax[j]]
        
        RT_ctr_pos.raw = KM1605_UvOX_p_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          KM1605_UvOX_p_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rt[j]]
        
        # retrieve data from ms2 scans which have the correct precursor mz (i.e., the mz of this instance of the LOBSTAHS-ID'd compound (peak) we are currently considering) AND were acquired within the RT window (raw min/max) of the peak
        
        # will use search window of 40(!) ppm
        # confirmed this is necessary based on manual inspection of data
        ms2_lkup_window.ppm = 40
        
        # first, gather possibly relevant precursor scans based strictly on mass difference
        
        # get an index to these scans
        possible.precursors_pos.ind =
          which(abs(KM1605_UvOX_p_xraw_pos[[samp_ID]]@msnPrecursorMz-xcms.peakdata_thisgroup_pos$mz[j])/
                  xcms.peakdata_thisgroup_pos$mz[j]*1000000 < ms2_lkup_window.ppm)
        
        # whittle this list to make sure the scans *also* fall within the rt window for the parent peak
        # will use a little bit of a buffer to capture any scans falling just outside of the RT range
        
        valid.precursors_pos.ind =
          possible.precursors_pos.ind[KM1605_UvOX_p_xraw_pos[[samp_ID]]@msnRt[possible.precursors_pos.ind]>
                                        (RT_min_pos.raw-20) &
                                        KM1605_UvOX_p_xraw_pos[[samp_ID]]@msnRt[possible.precursors_pos.ind]<
                                        (RT_max_pos.raw+20)]
        
        # to get the actual QE scan numbers corresponding to this index, run the below:
        # KM1605_UvOX_p_xraw_pos[[j]]@msnAcquisitionNum[valid.precursors_pos.ind]
        
        # record the number of valid ms2 spectra in PAL1314_LMG1401.detected_pos_ion_fragments
        # add to record if there's information already recorded from another isomer; otherwise, simply overwrite the NA placeholder
        
        if (length(valid.precursors_pos.ind)>0) {
          
          if (is.na(KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID])) {
            
            KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] = length(valid.precursors_pos.ind)
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] = 
              KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] + length(valid.precursors_pos.ind)
            
          }
          
        } else {
          
          # there is no ms2 data for this parent, at least not how we went about finding it
          
          if (is.na(KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID])) {
            
            KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] = 0
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] = 
              KM1605_UvOX_p.detected_pos_ion_fragments[i,samp_ID] + 0
            
          }
          
        }
        
        # provide ourselves with an escape here in case there are no valid ms2 scans
        
        if (length(valid.precursors_pos.ind)>0) {
          
          # now, pull out the ms2 data for each of these scans
          # msnScanindex is constructed in such a way that we can recreate the peaklist for
          # a given scan using the following syntax
          
          relevant_ms2data = apply(as.matrix(valid.precursors_pos.ind),1,get.ms2Peaklist,sample_ID = samp_ID,
                                   xcmsRaw.list = KM1605_UvOX_p_xraw_pos)
          
          ### now, can get onto the business of actually examining the spectra for the diagnostic transitions ###
          
          ### scenario 1: class type is diagnosed via presence of product ion ###
          
          if (KM1605_UvOX_part.frag_lookup_classes$Frag_exp_type[
            rownames(KM1605_UvOX_part.frag_lookup_classes)==this.IDclass] == "PI") {
            # this class is diagnosed via presence of a product ion
            
            # apply some logic for product ion scenario: assume that for a product ion-based ID to be "good", we must observe a feature with the mz of the diagnostic ion (+/- some mz tolerance) as one of the top N peaks (by intensity) in at least one of the relevant ms2 scans
            
            # so, extract the top N (right now, 20) most intense features in each scan
            
            top_features.PI = lapply(relevant_ms2data,get.topN,20)
            
            # evaluate: do the list(s) of the top N most intense fragments contain the diagnostic ion?
            
            PI.eval_result = lapply(top_features.PI, eval.PIspecies, species = this.IDclass, ppm = 12,
                                    ms2.lookupClasses = KM1605_UvOX_part.frag_lookup_classes)
            
            # record result
            
            if (is.na(KM1605_UvOX_p.fragdata_results[i,samp_ID])) {
              
              KM1605_UvOX_p.fragdata_results[i,samp_ID] = sum(unlist(PI.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              KM1605_UvOX_p.fragdata_results[i,samp_ID] = 
                KM1605_UvOX_p.fragdata_results[i,samp_ID] + sum(unlist(PI.eval_result))
              
            }
            
          } else if (KM1605_UvOX_part.frag_lookup_classes$Frag_exp_type[
            rownames(KM1605_UvOX_part.frag_lookup_classes)==this.IDclass] == "CNL") {
            
            ### scenario 2: class type is diagnosed via constant neutral loss ###
            
            # assume it's a good ID in this case as long as an ion corresponding to the diagnostic CNL is one of the top N (50, for now) peaks (by intensity) in the + mode ms2 spectrum
            
            top_features.CNL = lapply(relevant_ms2data,get.topN,50)
            
            # evaluate & record
            
            CNL.eval_result = lapply(top_features.CNL, eval.CNLspecies, species = this.IDclass, sample_ID = samp_ID, 
                                     ppm = 12, ms2.lookupClasses = KM1605_UvOX_part.frag_lookup_classes)
            
            if (is.na(KM1605_UvOX_p.fragdata_results[i,samp_ID])) {
              
              KM1605_UvOX_p.fragdata_results[i,samp_ID] = sum(unlist(CNL.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              KM1605_UvOX_p.fragdata_results[i,samp_ID] = 
                KM1605_UvOX_p.fragdata_results[i,samp_ID] + sum(unlist(CNL.eval_result))
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

# append summary, presence/absence results to KM1605_UvOX_Expt_pos_withoddFA_particulate matrix

KM1605_UvOX_p.fragdata_result.summary =
  apply(KM1605_UvOX_p.fragdata_results, 1, sumfrag)
KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_conf =
  KM1605_UvOX_p.fragdata_result.summary

KM1605_UvOX_p.fragdata_presence.summary =
  apply(KM1605_UvOX_p.detected_pos_ion_fragments, 1, sumfrag)
KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_present =
  KM1605_UvOX_p.fragdata_presence.summary

### ***** additional confirmation by presence in negative mode data ***** ####

# as an additional means of confirmation, check to see whether the same IDs were made in negative ion mode

# preallocate a matrix to keep track of results
KM1605_UvOX_p.neg_mode_conf = as.data.frame(matrix(NA,nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate),ncol=1))

# load negative mode LOBSet
load("data/nice/Orbi_MS_data/LOBSTAHS_processed/KM1605_UvOX_Expt_neg_withoddFA_LOBSet.RData")
KM1605_UvOX_Expt_neg_withoddFA_particulate = getLOBpeaklist(KM1605_UvOX_Expt_neg_withoddFA_LOBSet) # generate peaklist

# iterate through the positive mode dataset and compare, then record results

for (i in 1:nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate)) {
  
  # any matches by compound ID
  
  compound_ID.matches = KM1605_UvOX_Expt_neg_withoddFA_particulate[KM1605_UvOX_Expt_pos_withoddFA_particulate$compound_name[i] == KM1605_UvOX_Expt_neg_withoddFA_particulate$compound_name,]
  
  if (any(abs(KM1605_UvOX_Expt_pos_withoddFA_particulate$peakgroup_rt[i]-
              compound_ID.matches$peakgroup_rt)<20)) {
    
    # the negative ion mode match belongs to a peakgroup that has a retention time within 15 seconds of the positive ion mode peakgroup --> call it good
    
    KM1605_UvOX_p.neg_mode_conf[i,1] = 1
    
  } else {
    
    KM1605_UvOX_p.neg_mode_conf[i,1] = 0
    
  }
  
}

# append cross-mode comparison results to the positive mode ID table

KM1605_UvOX_Expt_pos_withoddFA_particulate$neg_mode.conf =
  KM1605_UvOX_p.neg_mode_conf[,1]

# create a field (and populate partially) to assist in curation of the dataset:
# which records will be retained, which will be excised moving forward

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto = NA

# flag for removal certain entire classes which we aren't concerned with at the moment

# (but first, calculate % of total peak area accounted for by TAG)
# .
# .
# .
mean(apply(KM1605_UvOX_Expt_pos_withoddFA_particulate[KM1605_UvOX_Expt_pos_withoddFA_particulate$lipid_class=="TAG",c(16:20,23:28)],2,sum,na.rm=T)/
       apply(KM1605_UvOX_Expt_pos_withoddFA_particulate[,c(16:20,23:28)],2,sum,na.rm=T))

sd(apply(KM1605_UvOX_Expt_pos_withoddFA_particulate[KM1605_UvOX_Expt_pos_withoddFA_particulate$lipid_class=="TAG",c(16:20,23:28)],2,sum,na.rm=T)/
     apply(KM1605_UvOX_Expt_pos_withoddFA_particulate[,c(16:20,23:28)],2,sum,na.rm=T))
# .
# .
# .
# now, actually flag the TAG & some other lipid classes for removal 
KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$lipid_class %in% 
    c("TAG","pigment","PUA","IP_MAG")] = 1

# flag for retention IDs that are *definitely* valid 
KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_conf>0 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$neg_mode.conf>0] = 0

# flag for retention IDs that passed the ms2 criterion but weren't made in the negative mode data
# (essentially, we are deferring to the superior diagnostic power of the fragmentation spectra)
KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_conf>0 &
    (is.na(KM1605_UvOX_Expt_pos_withoddFA_particulate$neg_mode.conf) || 
       KM1605_UvOX_Expt_pos_withoddFA_particulate$neg_mode.conf == 0)] = 0

# flag for removal any putative isomer of an ID we we've already declared valid

for (i in 1:nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate)) {
  
  if (!is.na(KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[i])) {
    
    if (KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[i] == 0) {
      
      KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
        KM1605_UvOX_Expt_pos_withoddFA_particulate$match_ID!=
          KM1605_UvOX_Expt_pos_withoddFA_particulate$match_ID[i] &
          KM1605_UvOX_Expt_pos_withoddFA_particulate$LOBdbase_mz==
          KM1605_UvOX_Expt_pos_withoddFA_particulate$LOBdbase_mz[i] &
          KM1605_UvOX_Expt_pos_withoddFA_particulate$peakgroup_rt==
          KM1605_UvOX_Expt_pos_withoddFA_particulate$peakgroup_rt[i] &
          is.na(KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto)] = 1
    }
    
  }
  
}

# flag for removal any putative isomer of an ID we suspect is valid by additional presence in negative mode data (but which did not necessarily pass the ms2 criterion)

for (i in 1:nrow(KM1605_UvOX_Expt_pos_withoddFA_particulate)) {
  
  if (KM1605_UvOX_Expt_pos_withoddFA_particulate$neg_mode.conf[i] == 1 &
      (is.na(KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_conf[i]) || 
       KM1605_UvOX_Expt_pos_withoddFA_particulate$ms2_conf[i] == 0)) {
    
    KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
      KM1605_UvOX_Expt_pos_withoddFA_particulate$match_ID!=
        KM1605_UvOX_Expt_pos_withoddFA_particulate$match_ID[i] &
        KM1605_UvOX_Expt_pos_withoddFA_particulate$LOBdbase_mz==
        KM1605_UvOX_Expt_pos_withoddFA_particulate$LOBdbase_mz[i] &
        KM1605_UvOX_Expt_pos_withoddFA_particulate$peakgroup_rt==
        KM1605_UvOX_Expt_pos_withoddFA_particulate$peakgroup_rt[i] &
        is.na(KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto)] = 1
    
  }
  
}

# since this is marine environmental data collected at 0.2 um, bacteria are present; therefore any length fatty acid > C13 and < C26 is possible
# however, we will assume that: any C13 or C15 fatty acids must be fully saturated, a C14 can be at most monounsaturated, and a C16 can have at least 4 double bonds
# also, assume that the max # of acyl DB for a single IP-DAG will be 12 (C22:6 + C22:6... C26, if it is present, will be fully saturated)

# flag for removal any putative IDs that do not meet these minimum requirements

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C < 26] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 26 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 0
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 27 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 1
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 28 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 2
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 29 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 1
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 30 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 5
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 31 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 6
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C == 32 &
    KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 8
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_C > 52
  ] = 1

KM1605_UvOX_Expt_pos_withoddFA_particulate$rm.flag_auto[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$FA_total_no_DB > 12
  ] = 1

# if desired, write the results to a csv file

write.csv(KM1605_UvOX_Expt_pos_withoddFA_particulate, file = "KM1605_UvOX_Expt_pos_withoddFA_particulate.csv")

#### further pre-processing ####

# now, extract subset of data (IP-DAG, TAG, oxidized forms), plus DNPPE

KM1605_UvOX_Expt_pos_part.sub = KM1605_UvOX_Expt_pos_withoddFA_particulate[
  KM1605_UvOX_Expt_pos_withoddFA_particulate$lipid_class %in% c("IP_DAG","DNPPE","TAG"),]

KM1605_UvOX_Expt_pos_part.sub.noDGCC = KM1605_UvOX_Expt_pos_part.sub[
  KM1605_UvOX_Expt_pos_part.sub$species!="DGCC",]

# define classes for which concentrations are to be calculated, the models, and
# the cutoffs in case of split prediction

KM1605_UvOX_Expt_pos_part.conc_classes = as.data.frame(matrix(NA,9,4))
colnames(KM1605_UvOX_Expt_pos_part.conc_classes) =
  c("Lipid_class","Model.low","Model.hi","Cutoff_PA")

KM1605_UvOX_Expt_pos_part.conc_classes[,1] =
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGTS_DGTA","DNPPE","TAG")

KM1605_UvOX_Expt_pos_part.conc_classes[,2] =
  c("linfit_low.PG.20161107","linfit_low.PE.20161107","linfit_low.PC.20161107",
    "linfit_low.MGDG.20161107","linfit_low.SQDG.20161107",
    "linfit_low.DGDG.20161107","linfit_low.DGTS_DGTA.20161005",
    "linfit_low.DNPPE.20161107","nls.fit.TAG20161004")

KM1605_UvOX_Expt_pos_part.conc_classes[,3] =
  c("linfit_hi.PG.20161107","linfit_hi.PE.20161107","linfit_hi.PC.20161107",
    "linfit_hi.MGDG.20161107","linfit_hi.SQDG.20161107",
    "linfit_hi.DGDG.20161107","linfit_hi.DGTS_DGTA.20161005",
    "linfit_hi.DNPPE.20161107","nls.fit.TAG20161004")

KM1605_UvOX_Expt_pos_part.conc_classes[,4] =
  c(PG_std_breakpoint.20161107,PE_std_breakpoint.20161107,PC_std_breakpoint.20161107,
    MGDG_std_breakpoint.20161107,SQDG_std_breakpoint.20161107,DGDG_std_breakpoint.20161107,
    DGTS_DGTA_std_breakpoint.20161005,DNPPE_std_breakpoint.20161107,NA)

# first, calculate pmol o.c.
# preallocate matrix for result
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc = KM1605_UvOX_Expt_pos_part.sub.noDGCC

# calculate for each lipid class, using appropriate standard curve

for (i in 1:nrow(KM1605_UvOX_Expt_pos_part.conc_classes)) {
  
  # first, calculate pmol o.c.
  
  if (KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i]!="TAG") {
    
  pmol.oc.thisclass =
    apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc[
      KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc$species==
        KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i],16:28],c(1,2),splitpred,
      eval(parse(text = KM1605_UvOX_Expt_pos_part.conc_classes$Model.low[i])),
      eval(parse(text = KM1605_UvOX_Expt_pos_part.conc_classes$Model.hi[i])),
      KM1605_UvOX_Expt_pos_part.conc_classes$Cutoff_PA[i])
  
  # store result as appropriate
  
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc[
    KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc$species==
      KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i],16:28] = pmol.oc.thisclass
  
  } else if (KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i]=="TAG") {
    
    TAGdata = KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc[
      KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc$species==
        KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i],]
    
    # preallocate
    
    TAGconcs = matrix(NA,nrow(TAGdata),length(16:28))
    
    for (j in 1:nrow(TAGconcs)) {
      
      TAGconcs[j,] =
        sapply(TAGdata[j,16:28],TAGconcCalc,
               numC=TAGdata$FA_total_no_C[j],
               numDB=TAGdata$FA_total_no_DB[j],
               TAGfitFunction = eval(parse(text = KM1605_UvOX_Expt_pos_part.conc_classes$Model.low[i])),
               linDNPPEfit_low = eval(parse(text = KM1605_UvOX_Expt_pos_part.conc_classes$Model.low[
                 KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class=="DNPPE"])),
               linDNPPEfit_hi = eval(parse(text = KM1605_UvOX_Expt_pos_part.conc_classes$Model.hi[
                 KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class=="DNPPE"])),
               DNPPEfit.cutoff = KM1605_UvOX_Expt_pos_part.conc_classes$Cutoff_PA[
                 KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class=="DNPPE"]
               )
      
    }
      
    # store result as appropriate
    
    KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc[
      KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc$species==
        KM1605_UvOX_Expt_pos_part.conc_classes$Lipid_class[i],16:28] = TAGconcs
    
  }
  
}

# now, scale pmol o.c. to pmol per sample (and then volume) using DNPPE recovery standard added at time of extraction, and volumes recorded when samples were sacrificed

DNPPE_BD_KM1605_UvOx_uL = 20 # amount DNPPE added per sample in uL, per VML B&D protocol

DNPPE_pmol_added_per_samp = DNPPE_mg_mL_BD_extracts_2016*(1/DNPPE_MW)*(10^9)*(1/10^3)*DNPPE_BD_KM1605_UvOx_uL

KM1605_UvOx_pos_part_DNPPE.samp.RF = DNPPE_pmol_added_per_samp/KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc[
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc$compound_name=="DNPPE",16:28]  # recovery factor

# avg volume filtered for KM1605 samples (from BVM KM1605 notebook)

KM1605_filter_vol_mL = 1000

# create final results data frame
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L = KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.oc

# apply RF to samples
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28] = sweep(as.matrix(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28]), 2, as.numeric(KM1605_UvOx_pos_part_DNPPE.samp.RF), "*") 

# convert pmol to pmol/L
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28] = sweep(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28], 2, KM1605_filter_vol_mL/1000, "/")  # calculate pmol/mL, using correct volumes

# need to simplify the dataset a bit and do some QA

# eliminate features w/calibrated mass < 0 

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28] =
  replace(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28],
          KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28]<0,
          NA)

# now remove elements for which there is no data in any sample (includes elements
# just reduced to NA)
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L = 
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[,16:28],1,sum,na.rm = T)>0,]

# can remove DNPPE at this point
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L = 
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[!(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L$compound_name=="DNPPE"),]

# another opportunity to write results to file
write.csv(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L,
          file="KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.csv")

# *should* do the manual annotation commented out, but in the interest of time (12/1/16), going to skip for now and do something quick and dirty (for SCOPE Dec 16 meeting poster), like just use everything still present that isn't an automatic eliminate

# # *** at this point, manual curation (based on consideration of all available LOBSTAHS result codes) should be performed on the data exported to KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.csv
# # below, we reimport the results, which are contained in .csv file UNC_Marchetti_diatom_cultures_IP-DAG_pmol_total_man_ann.csv (exported from first worksheet of the .xlsx file UNC_Marchetti_diatom_cultures_IP-DAG_pmol_totals.xlsx after curation was finished)
# # the field rm.flag_manual contains the results of the manual annotation and should be combined with annotations that were already in rm.flag_auto at time of saving the .csv file above
# 
# # import & extract codes from the curation, combine with codes in rm.flag_auto, then apply to the dataset to produce our final list of IDs
# 
# Marchetti_cultures_man_annot = read.csv("data/nice/LOBSTAHS_lipid_identities/UNC_Marchetti_diatom_cultures_IP-DAG_pmol_total_man_ann.csv",
#                                         skip = 0)
# 
# Marchetti_diatom_cultures_pos.pmol_total.fully_annotated = 
#   Marchetti_diatom_cultures_pos.unox_IPL.noDGCC.pmol.total
# Marchetti_diatom_cultures_pos.pmol_total.fully_annotated$rm.flag_auto =
#   Marchetti_cultures_man_annot$rm.flag_auto
# Marchetti_diatom_cultures_pos.pmol_total.fully_annotated$rm.flag_manual =
#   Marchetti_cultures_man_annot$rm.flag_manual
# 
# Marchetti_diatom_cultures_pos.pmol_total.final = 
#   Marchetti_diatom_cultures_pos.pmol_total.fully_annotated[
#     !apply(Marchetti_diatom_cultures_pos.pmol_total.fully_annotated[,c("rm.flag_auto","rm.flag_manual")],1,sum,na.rm = T),]
# 
# # export these final IDs
# 
# write.csv(Marchetti_diatom_cultures_pos.pmol_total.final,
#           file = "UNC_Marchetti_diatom_cultures_IP-DAG_pmol_totals.final.csv")
# 
# # extract, calculate, append total # of C atoms in each ID'd molecule
# 
# Marchetti_diatom_cultures_pos.pmol_total.final$total_no_C =
#   as.numeric(sapply(Marchetti_diatom_cultures_pos.pmol_total.final$elem_formula,
#                     function(x) as.numeric(str_match(as.character(str_match(x, "^C[0-9]*H")), "[0-9]+"))))
# 

# get a subset

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster =
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L[KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L$lipid_class=="TAG" |
                                                (KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L$lipid_class=="IP_DAG" &
                                                KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L$rm.flag_auto %in% c(NA,0))
                                                ,]

# eliminate MGDG that are really oxidized TAG
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster =
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[-c(which(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$species=="MGDG" & 
                                                       KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$peakgroup_rt>(16.5*60) & KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$LOBdbase_mz %in%
         KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$LOBdbase_mz[KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$lipid_class=="TAG"])),]

# a subset without TAGs

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.noTAG = KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$lipid_class!="TAG",]

# some stats

# fraction oxidized

KM1605_UvOX_ox.percent = apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[!(
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$degree_oxidation==0),c(16:20,23:28)],2,sum,na.rm = T)/
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,c(16:20,23:28)],2,sum,na.rm = T)

mean(KM1605_UvOX_ox.percent[1:3])
sd(KM1605_UvOX_ox.percent[1:3])

mean(KM1605_UvOX_ox.percent[4:5])
sd(KM1605_UvOX_ox.percent[4:5])

mean(KM1605_UvOX_ox.percent[6:8])
sd(KM1605_UvOX_ox.percent[6:8])

# one +UVB sample was not filtered correctly; ignore it (per Ben's field notebook)
mean(KM1605_UvOX_ox.percent[c(9:10)])
sd(KM1605_UvOX_ox.percent[c(9,10)])

# barplot

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "KM1605_UvOx_ox_frac_barplot.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey")
treatments <- c("Dark","- UVB","+ UVB")

bardata.mean <- matrix(NA,3,2)

bardata.mean[1,] = c(mean(KM1605_UvOX_ox.percent[4:5]),mean(KM1605_UvOX_ox.percent[1:3]))
  
bardata.mean[2,2] = mean(KM1605_UvOX_ox.percent[6:8])

bardata.mean[3,2] = mean(KM1605_UvOX_ox.percent[c(9:10)])

bardata.se <- matrix(NA,3,2)

bardata.se[1,] = c(sd(KM1605_UvOX_ox.percent[4:5])/sqrt(2),sd(KM1605_UvOX_ox.percent[1:3])/sqrt(3))

bardata.se[2,2] = sd(KM1605_UvOX_ox.percent[6:8])/sqrt(3)

bardata.se[3,2] = sd(KM1605_UvOX_ox.percent[c(9,10)])/sqrt(2)

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,.7), ylab=paste("Mole fraction identified as oxidized lipid"), xlab="Timepoint", names.arg=c("Initial","+ 24 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topleft",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()

# fraction as TAG

KM1605_UvOX_TAG.percent = apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$lipid_class=="TAG",c(16:20,23:28)],2,sum,na.rm = T)/
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,c(16:20,23:28)],2,sum,na.rm = T)

mean(KM1605_UvOX_TAG.percent[1:3])
sd(KM1605_UvOX_TAG.percent[1:3])

mean(KM1605_UvOX_TAG.percent[4:5])
sd(KM1605_UvOX_TAG.percent[4:5])

mean(KM1605_UvOX_TAG.percent[6:8])
sd(KM1605_UvOX_TAG.percent[6:8])

mean(KM1605_UvOX_TAG.percent[c(9:10)])
sd(KM1605_UvOX_TAG.percent[c(9,10)])

# barplot

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "KM1605_UvOx_TAG_frac_barplot.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey")
treatments <- c("Dark","- UVB","+ UVB")

bardata.mean <- matrix(NA,3,2)

bardata.mean[1,] = c(mean(KM1605_UvOX_TAG.percent[4:5]),mean(KM1605_UvOX_TAG.percent[1:3]))

bardata.mean[2,2] = mean(KM1605_UvOX_TAG.percent[6:8])

bardata.mean[3,2] = mean(KM1605_UvOX_TAG.percent[c(9:10)])

bardata.se <- matrix(NA,3,2)

bardata.se[1,] = c(sd(KM1605_UvOX_TAG.percent[4:5])/sqrt(2),sd(KM1605_UvOX_TAG.percent[1:3])/sqrt(3))

bardata.se[2,2] = sd(KM1605_UvOX_TAG.percent[6:8])/sqrt(3)

bardata.se[3,2] = sd(KM1605_UvOX_TAG.percent[c(9,10)])/sqrt(2)

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,.5), ylab=paste("Mole fraction identified as TAG"), xlab="Timepoint", names.arg=c("Initial","+ 24 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topleft",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()

# calculate some treatment & timepoint means

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means= as.data.frame(
  matrix(NA, nrow(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster), ncol = 4))

colnames(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means) =
  c("Initial_pmol_L","Final_dark_control_pmol_L","Final_plus_UVB_pmol_L",
    "Final_minus_UVB_pmol_L")

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means$Initial_pmol_L =
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,19:20],1,mean,na.rm = T)

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means$Final_dark_control_pmol_L =
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,16:18],1,mean,na.rm = T)

# again, eliminating bad +UVB sample
KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means$Final_plus_UVB_pmol_L =
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,c(26,27)],1,mean,na.rm = T)

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means$Final_minus_UVB_pmol_L =
  apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,23:25],1,mean,na.rm = T)

# append these to the main matrix

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster =
cbind(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster,KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster.means)

# some stats

# total lipid mass

KM1605_UvOX_ox.total_mass = apply(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[
  ,c(16:20,23:28)],2,sum,na.rm = T)

mean(KM1605_UvOX_ox.total_mass[1:3])
sd(KM1605_UvOX_ox.total_mass[1:3])

mean(KM1605_UvOX_ox.total_mass[4:5])
sd(KM1605_UvOX_ox.total_mass[4:5])

mean(KM1605_UvOX_ox.total_mass[6:8])
sd(KM1605_UvOX_ox.total_mass[6:8])

# +UVB #2 sample was not filtered correctly; discard it
mean(KM1605_UvOX_ox.total_mass[c(9,10)])
sd(KM1605_UvOX_ox.total_mass[c(9,10)])

# calculate stats & evaluate signficant

  Total_mass.statsub = as.data.frame(KM1605_UvOX_ox.total_mass)
  
  Total_mass.statsub$Treatment = c("Dark_control_final",
                                          "Dark_control_final",
                                          "Dark_control_final",
                                          "Dark_control_initial",
                                          "Dark_control_initial",
                                          "No_UVB_final",
                                          "No_UVB_final",
                                          "No_UVB_final",
                                          "Plus_UVB_final",
                                          "Plus_UVB_final",
                                          "Plus_UVB_final")
  
  colnames(Total_mass.statsub)[1] = c("Total_mass_mol")
  
  # eliminate bad sample
  
  Total_mass.statsub.good = Total_mass.statsub[rownames(Total_mass.statsub)!="Plus_UVB_final_QE002917",]
  
  Totalmass.mod = lm(Total_mass_mol ~ Treatment, data = Total_mass.statsub.good)
  
  print(anova(Totalmass.mod))
  Totalmass.aov = aov(Totalmass.mod)
  tukey = TukeyHSD(Totalmass.aov, conf.level = 0.95)
  print(tukey)
  
  # calculate & display mean differences ± SD for selected treatment pairs
  
  xbar_init = mean(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Dark_control_initial"])
  sd_init = sd(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Dark_control_initial"])
  se_init = sd_init/sqrt(length(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Dark_control_initial"]))
  
  xbar_final = mean(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Plus_UVB_final"])
  sd_final = sd(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Plus_UVB_final"])
  se_final = sd_final/sqrt(length(Total_mass.statsub.good$Total_mass_mol[Total_mass.statsub.good$Treatment=="Plus_UVB_final"]))
  
  delta = xbar_final-xbar_init
  uncert = sqrt(se_final^2 + se_init^2)
  
  cat("+ UVB vs. initial, mean ± uncertainty: ",delta," ± ",uncert,"\n")

# line plot

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "KM1605_UvOx_total_lipid_mass.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

barcolors <- c("black","lightgrey","darkgrey")
treatments <- c("Dark","- UVB","+ UVB")

bardata.mean <- matrix(NA,3,2)

bardata.mean[1,] = c(mean(KM1605_UvOX_ox.total_mass[4:5]),mean(KM1605_UvOX_ox.total_mass[1:3]))

bardata.mean[2,2] = mean(KM1605_UvOX_ox.total_mass[6:8])

bardata.mean[3,2] = mean(KM1605_UvOX_ox.total_mass[c(9:10)])

bardata.se <- matrix(NA,3,2)

bardata.se[1,] = c(sd(KM1605_UvOX_ox.total_mass[4:5])/sqrt(2),sd(KM1605_UvOX_ox.total_mass[1:3])/sqrt(3))

bardata.se[2,2] = sd(KM1605_UvOX_ox.total_mass[6:8])/sqrt(3)

bardata.se[3,2] = sd(KM1605_UvOX_ox.total_mass[c(9,10)])/sqrt(2)

plotTop <- max(bardata.mean+bardata.se,na.rm=T)
barCenters <- barplot(height=bardata.mean, beside=T,space = c(0.1,0.8), col=barcolors, las=1, ylim=c(0,8000), ylab=paste("Total lipids (pmol/L)"), xlab="Timepoint", names.arg=c("Initial","+ 24 h"))
segments(barCenters, bardata.mean-bardata.se, barCenters, bardata.mean+bardata.se, lwd=2) 

legend("topleft",treatments,fill=barcolors,xpd = TRUE) #inset = c(0,-.45))

dev.off()

#### create some heatmaps ####

# necessary libraries

library(gplots)

# extract data, but also a copy of the original DF with all the metadata
KM1605_UvOX.pos.part.pmol.L.HM = KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster[,
 c("Initial_pmol_L","Final_dark_control_pmol_L","Final_plus_UVB_pmol_L",
   "Final_minus_UVB_pmol_L")  
]

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat =
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster

# generate some good colnames and rownames

colnames(KM1605_UvOX.pos.part.pmol.L.HM) = 
  gsub(c("Initial_pmol_L"),c("Dark control, t = 0"),colnames(KM1605_UvOX.pos.part.pmol.L.HM))
colnames(KM1605_UvOX.pos.part.pmol.L.HM) = 
  gsub(c("Final_dark_control_pmol_L"),c("Dark control, t + 24 h"),colnames(KM1605_UvOX.pos.part.pmol.L.HM))
colnames(KM1605_UvOX.pos.part.pmol.L.HM) = 
  gsub(c("Final_plus_UVB_pmol_L"),c("+ UVB, t + 24 h"),colnames(KM1605_UvOX.pos.part.pmol.L.HM))
colnames(KM1605_UvOX.pos.part.pmol.L.HM) = 
  gsub(c("Final_minus_UVB_pmol_L"),c("- UVB, t + 24 h"),colnames(KM1605_UvOX.pos.part.pmol.L.HM))

# append RT data to row names

rownames(KM1605_UvOX.pos.part.pmol.L.HM) =
  apply(cbind(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$compound_name, as.character(round(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.poster$peakgroup_rt/60,2))),1,paste,collapse = ", RT ")

rownames(KM1605_UvOX.pos.part.pmol.L.HM) = apply(cbind(rownames(KM1605_UvOX.pos.part.pmol.L.HM),rep("min.",nrow(KM1605_UvOX.pos.part.pmol.L.HM))),1,paste,collapse = " ")

# change to matrices

KM1605_UvOX.pos.part.pmol.L.HM = as.matrix(KM1605_UvOX.pos.part.pmol.L.HM)

# eliminate any analyte with no data or all 0's

# set all NaN to zero
KM1605_UvOX.pos.part.pmol.L.HM[is.nan(KM1605_UvOX.pos.part.pmol.L.HM)]=0

KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat =
  KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat[c(apply(KM1605_UvOX.pos.part.pmol.L.HM,1,sum)!=0),]

KM1605_UvOX.pos.part.pmol.L.HM = 
  KM1605_UvOX.pos.part.pmol.L.HM[c(apply(KM1605_UvOX.pos.part.pmol.L.HM,1,sum)!=0),]

# fix NA's and small values

KM1605_UvOX.pos.part.pmol.L.HM[KM1605_UvOX.pos.part.pmol.L.HM==0] = 0.000000001

# calculate fold-change relative to initial

KM1605_UvOX.pos.part.pmol.L.HM.foldchange <- sweep(KM1605_UvOX.pos.part.pmol.L.HM, 1, KM1605_UvOX.pos.part.pmol.L.HM[,1], "/")

# # return any initial values (first row) that were originally NA's back to NA's
# 
# screenedpeaks_exptmeans.full_PC.byttp[1,screenedpeaks_exptmeans.full_PC.byttp[1,]==10] <- NA 

# transform into log2 fold change

KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2 <- log2(KM1605_UvOX.pos.part.pmol.L.HM.foldchange)

# # put the NaN's back
# 
# KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[KM1605_UvOX.pos.part.pmol.L.HM==0.000000001] = NaN

# screenedpeaks_exptmeans.full_PC.byttp.rev <- screenedpeaks_exptmeans.full_PC.byttp[dim(screenedpeaks_exptmeans.full_PC.byttp):1,]

# build heatmaps

# # first, combine all features into single matrix
# 
# Exp_13_features.HM.foldchange.log2 = rbind(Exp_13_PC.samp.pmol.mL.norm.HM.foldchange.log2,
#                                            Exp_13_FFA.neg.samp.pmol.mL.mean.HM.foldchange.log2,
#                                            Exp_13_LPC.samp.pmol.mL.neg.mean.HM.foldchange.log2)

# full heatmap

# generate breaks and colors

breaks.neg = seq(-max(abs(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2), na.rm=T),0,length.out=150)
breaks.pos = seq(0,max(abs(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2), na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-150] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

par(oma=c(0,0,0,0)) # set margins

pdf(file = "KM1605_heatmap_full.pdf",
    width = 7, height = 10, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

heatmap.2(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[,2:ncol(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(12,12))

dev.off()

# maybe just most abundant compounds in initial sample

KM1605_UVOx.part.abundant.log2 = KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[
  order(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat$Initial_pmol_L, 
        decreasing = TRUE),]

KM1605_UVOx.part.abundant.log2.top100 = KM1605_UVOx.part.abundant.log2[1:100,]
KM1605_UVOx.part.abundant.log2.top50 = KM1605_UVOx.part.abundant.log2[1:50,]
KM1605_UVOx.part.abundant.log2.top25 = KM1605_UVOx.part.abundant.log2[1:25,]

# generate breaks and colors

breaks.neg = seq(-max(abs(KM1605_UVOx.part.abundant.log2.top25), na.rm=T),0,length.out=150)
breaks.pos = seq(0,max(abs(KM1605_UVOx.part.abundant.log2.top25), na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-150] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

par(oma=c(0,0,0,0)) # set margins

pdf(file = "KM1605_heatmap_mostabundantinitial.pdf",
    width = 7, height = 10, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

heatmap.2(KM1605_UVOx.part.abundant.log2.top25[,2:ncol(KM1605_UVOx.part.abundant.log2.top25)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(14,14))

dev.off()

# maybe just the compounds exhibiting largest changes

# ... based on fold-change

KM1605_UVOx.part.bigdelta.log2 = KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[
  rownames(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2) %in%
  c(rep(rownames(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2),4))[
  order(abs(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2), decreasing = TRUE)[1:25]],]

# generate breaks and colors

breaks.neg = seq(min(KM1605_UVOx.part.bigdelta.log2, na.rm=T),0,length.out=150)
breaks.pos = seq(0,max(KM1605_UVOx.part.bigdelta.log2, na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-150] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

heatmap.2(KM1605_UVOx.part.bigdelta.log2[,2:ncol(KM1605_UVOx.part.bigdelta.log2)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(14,14))

# ... based on largest magnitude change

par(oma=c(0,0,0,0)) # set margins

pdf(file = "KM1605_heatmap.pdf",
        width = 7, height = 10, pointsize = 12,
        bg = "white")

par(mar=c(5,5,1,1))

KM1605_UVOx.part.bigdelta.alt.log2 = KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[
  order(abs(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat$Initial_pmol_L-
              KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat$Final_plus_UVB_pmol_L), 
        decreasing = TRUE)[1:25],]

# generate breaks and colors

breaks.neg = seq(min(KM1605_UVOx.part.bigdelta.alt.log2, na.rm=T),0,length.out=150)
breaks.pos = seq(0,max(KM1605_UVOx.part.bigdelta.alt.log2, na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-150] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

heatmap.2(KM1605_UVOx.part.bigdelta.alt.log2[,2:ncol(KM1605_UVOx.part.bigdelta.alt.log2)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(14,14))

dev.off()

# maybe just those that increased most in the +UVB treatment

par(oma=c(0,0,0,0)) # set margins

pdf(file = "KM1605_heatmap_2.pdf",
    width = 7, height = 10, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

KM1605_UVOx.part.bigdeltaUVB.log2 = KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[
  order(KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[,3],decreasing = TRUE)[1:25],]

# generate breaks and colors

# breaks.neg = seq(min(KM1605_UVOx.part.bigdeltaUVB.log2, na.rm=T),0,length.out=25)
breaks.pos = seq(0,max(KM1605_UVOx.part.bigdeltaUVB.log2, na.rm=T),length.out=150)
breaks.set = c(breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
# breaks.set <- breaks.set[-25] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

heatmap.2(KM1605_UVOx.part.bigdeltaUVB.log2[,2:ncol(KM1605_UVOx.part.bigdeltaUVB.log2)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(14,14))

dev.off()

# maybe just most abundant compounds in +UVB treatment

KM1605_UVOx.part.abundantUVB.log2 = KM1605_UvOX.pos.part.pmol.L.HM.foldchange.log2[
  order(KM1605_UvOX_Expt_pos_part.sub.noDGCC.pmol.L.HM.metdat$Final_plus_UVB_pmol_L, 
        decreasing = TRUE),]

KM1605_UVOx.part.abundantUVB.log2.top100 = KM1605_UVOx.part.abundantUVB.log2[1:100,]
KM1605_UVOx.part.abundantUVB.log2.top50 = KM1605_UVOx.part.abundantUVB.log2[1:50,]

# generate breaks and colors

breaks.neg = seq(min(KM1605_UVOx.part.abundantUVB.log2.top100, na.rm=T),0,length.out=400)
breaks.pos = seq(0,max(KM1605_UVOx.part.abundantUVB.log2.top100, na.rm=T),length.out=150)
breaks.set = c(breaks.neg,breaks.pos)
# breaks.set <- breaks.set[-length(breaks.set)/2]
breaks.set <- breaks.set[-400] # remove the duplicated "0" break point

gradient1 = colorpanel( sum( breaks.set[-1]<=0 ), "darkred", "white" )
gradient2 = colorpanel( sum( breaks.set[-1]>0 ), "white", "darkblue" )
hm.colors = c(gradient1,gradient2)

heatmap.2(KM1605_UVOx.part.abundantUVB.log2.top100[,2:ncol(KM1605_UVOx.part.abundantUVB.log2.top100)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
          density.info=c("none"),margins=c(14,14))

# # maybe just PC 22:6 and derivatives
# 
# # create our subset
# Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2[
#   grepl(paste(c("PC 22:6",
#                 "FFA 22:6",
#                 "LPC 22:6"), collapse = "|"), rownames(Exp_13_features.HM.foldchange.log2)),]
# 
# # also, just want the -HB controls (not the +HB controls)
# Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2.sub226[,-c(4:5)]
# 
# # lastly, get rid of a duplicate 22:6 FFA (ID did not stand up to manual examination of spectra)
# Exp_13_features.HM.foldchange.log2.sub226 = Exp_13_features.HM.foldchange.log2.sub226[-c(10),]
# 
# # create a matrix (yes, manually based on inspection of output from
# # PAL1314_liposome_expts.R - I know, inefficient) of p-values for each block in the
# # subsetted heatmap
# 
# Exp_13_features.HM.foldchange.log2.sub226.pvals = 
#   matrix(nrow = nrow(Exp_13_features.HM.foldchange.log2.sub226),
#          ncol = ncol(Exp_13_features.HM.foldchange.log2.sub226)-1)
# 
# # manually populate based on Tukey HSD test output from PAL1314_liposome_expts.R
# # will use "1" as not significant indicator
# Exp_13_features.HM.foldchange.log2.sub226.pvals[1,] =
#   c(1,1,1,0.05,1,0.01,1,0.01) # PC 22:6, 22:6
# Exp_13_features.HM.foldchange.log2.sub226.pvals[2,] =
#   c(1,1,1,0.01,1,0.01,1,0.05) # PC 44:12 +2O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[3,] =
#   c(1,1,1,0.05,1,0.01,1,1) # PC 44:12 +4O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[4,] =
#   c(rep(1,8)) # PC 44:12 +1O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[5,] =
#   c(rep(1,8)) # PC 44:12 +3O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[6,] =
#   c(rep(1,8)) # PC 44:12 +3O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[7,] =
#   c(rep(1,8)) # FFA 22:6
# Exp_13_features.HM.foldchange.log2.sub226.pvals[8,] =
#   c(1,1,1,1,1,0.0001,1,0.05) # FFA 22:6 +2O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[9,] =
#   c(1,1,1,1,1,0.0001,1,1) # FFA 22:6 +3O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[10,] =
#   c(1,1,1,1,1,0.01,1,1) # FFA 22:6 +10
# Exp_13_features.HM.foldchange.log2.sub226.pvals[11,] =
#   c(1,1,1,1,1,0.001,1,1) # LPC 22:6 +4O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[12,] =
#   c(1,1,1,0.05,1,0.001,1,0.05) # LPC 22:6 +2O
# Exp_13_features.HM.foldchange.log2.sub226.pvals[13,] =
#   c(1,1,1,1,1,0.05,1,1) # LPC 22:6
# Exp_13_features.HM.foldchange.log2.sub226.pvals[14,] =
#   c(1,1,1,1,1,0.01,1,1) # LPC 22:6 +1O
# 
# # create a matrix of actual symbols to be plotted
# 
# Exp_13_features.HM.sigsymbols = Exp_13_features.HM.foldchange.log2.sub226.pvals
# Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==1]=""
# Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.05]="+"
# Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.01]="*"
# Exp_13_features.HM.sigsymbols[Exp_13_features.HM.sigsymbols==0.0001]="***"
# 
# par(oma=c(0,0,0,0)) # set margins
# 
# pdf(file = "Exp13_heatmap_PC22-6plus.pdf",
#     width = 7, height = 9.5, pointsize = 12,
#     bg = "white")
# 
# #par(mar=c(5,5,1,1))
# 
# heatmap.2(Exp_13_features.HM.foldchange.log2.sub226[,2:ncol(Exp_13_features.HM.foldchange.log2.sub226)],breaks=breaks.set,col=hm.colors,scale="none",Colv=TRUE,trace="none",Rowv=TRUE,dendrogram="both",na.color="grey",key=T,
#           density.info=c("none"),margins=c(16,16),cellnote=Exp_13_features.HM.sigsymbols)
# 
# dev.off()
# 
# # now, a small heatmap of just the other unoxidized parent lipids in this same
# # experiment, for comparison purposes
# 
# Exp_13_features.HM.foldchange.log2.subparents = Exp_13_features.HM.foldchange.log2[
#   grepl(paste(c("PC 18:1, 18:1,",
#                 "PC 18:0, 18:0,",
#                 "PC 16:0, 16:0,"), collapse = "|"), rownames(Exp_13_features.HM.foldchange.log2)),]
# 
# # get rid of the +HB controls
# Exp_13_features.HM.foldchange.log2.subparents = Exp_13_features.HM.foldchange.log2.subparents[,-c(4:5)]
# 
# # will need to order the treatments manually to fit with dendrogram clustering
# # of bigger heatmap
# 
# Exp_13_features.HM.foldchange.log2.subparents = 
#   Exp_13_features.HM.foldchange.log2.subparents[,c(2,3,8,6,4,5,7,9)]
# 
# # put species in a more logical order
# 
# Exp_13_features.HM.foldchange.log2.subparents = 
#   Exp_13_features.HM.foldchange.log2.subparents[c(1,3,2),]
# 
# par(oma=c(0,0,0,0)) # set margins
# 
# pdf(file = "Exp13_heatmap_other_parents.pdf",
#     width = 7, height = 10, pointsize = 12,
#     bg = "white")
# 
# #par(mar=c(5,5,1,1))
# 
# heatmap.2(Exp_13_features.HM.foldchange.log2.subparents[,1:ncol(Exp_13_features.HM.foldchange.log2.subparents)],breaks=breaks.set,col=hm.colors,scale="none",Colv=FALSE,trace="none",Rowv=FALSE,dendrogram="none",na.color="grey",key=T,
#           density.info=c("none"),margins=c(16,16))
# 
# dev.off()
# 
# 
# 
