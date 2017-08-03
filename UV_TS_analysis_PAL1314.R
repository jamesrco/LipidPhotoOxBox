# UV_TS_analysis_PAL1314.R
#
# Purpose: Analysis & treatment of UV timeseries data from Palmer 2013-2014
# field season. Uses NOAA ESRL GMD Antarctic UV irradiance data and data from
# Jaz UV-VIS spectrophotometer J.R.C. deployed at Palmer Station during the 
# 2013-2014 field season
#
# Also, calculations of the diffuse attenuation coefficient, Kd, from
# individual JAZ scans
#
# Created: 9/2/16 by James Collins, james.r.collins@aya.yale.edu
# Released under MIT License
#
# For workup and preliminary analysis of the JAZ data, see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem
#
# See Igor Pro file for generation of split-axis timeseries plot
#
# See http://esrl.noaa.gov/gmd/grad/antuv/docs/version0/DescDailyDose.html for
# NOAA ESRL GMD data and field descriptions

# This analysis is based on NOAA's Version 2 (corrected) data product

# need two libraries for the import & parsing of the individual NOAA spectra files
library(stringr)
library(RSEIS)

setwd("/Users/jamesrco/Code/LipidPhotoOxBox")
base.wd = getwd()

# some definitions

W_per_uW = 1/1000000
J_per_kJ = 1/1000
cm2_per_m2 = 10000

##### daily UVA, UVB dosages from NOAA data ##### 

# NOAA data

# set wd to data location
setwd("data/raw/ESRL_GMD_AntUV/ver2") 

# read in NOAA ESRL AntUV data (Version 2, daily dosage data), subset as needed
PAL_AntUV_DD_v2 = read.csv("PAL_all_DailyDose_no_flags_SZA_max_set_to82.csv", 
                  stringsAsFactors = FALSE, skip = 2)

# just the 290-315 nm wavelength dosage data
DD_UVB_incident = as.data.frame(cbind(PAL_AntUV_DD_v2$Date,PAL_AntUV_DD_v2$E290.315))
colnames(DD_UVB_incident) = c("Date_raw_julian","E290.315_kJ_m2")

# the 315-400 nm wavelength dosage data
DD_UVA_incident = as.data.frame(cbind(PAL_AntUV_DD_v2$Date,PAL_AntUV_DD_v2$E315.400))
colnames(DD_UVA_incident) = c("Date_raw_julian","E315.400_kJ_m2")

# convert dates from messed up Excel format
DD_UVB_incident$Date = as.POSIXct('1899-12-30', tz = "GMT")+(DD_UVB_incident$Date_raw_julian*24*60*60)
DD_UVB_incident$Date_simple = strptime(DD_UVB_incident$Date,"%Y-%m-%d", tz = "GMT")

DD_UVA_incident$Date = as.POSIXct('1899-12-30', tz = "GMT")+(DD_UVA_incident$Date_raw_julian*24*60*60)
DD_UVA_incident$Date_simple = strptime(DD_UVA_incident$Date,"%Y-%m-%d", tz = "GMT")

# subset, plot 2013-2014 data
DD_UVB_1314_incident = DD_UVB_incident[DD_UVB_incident$Date>as.POSIXct('2013-06-22') &
                                         DD_UVB_incident$Date<as.POSIXct('2014-06-22'),]

DD_UVA_1314_incident = DD_UVA_incident[DD_UVA_incident$Date>as.POSIXct('2013-06-22') &
                                         DD_UVA_incident$Date<as.POSIXct('2014-06-22'),]

plot(DD_UVB_1314_incident$Date,DD_UVB_1314_incident$E290.315_kJ_m2, "l", lwd = 1)

# save the data files
save(DD_UVB_1314_incident, file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_UVB_dose_incident_NOAA_ESRL_AntUV_DD_PAL1314.RData"))
save(DD_UVA_1314_incident, file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_UVA_dose_incident_NOAA_ESRL_AntUV_DD_PAL1314.RData"))

# export to .csv
write.csv(DD_UVB_1314_incident, file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_UVB_dose_incident_NOAA_ESRL_AntUV_DD_PAL1314.csv"))
write.csv(DD_UVA_1314_incident, file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_UVA_dose_incident_NOAA_ESRL_AntUV_DD_PAL1314.csv"))

# JAZ data

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

setwd(base.wd)
setwd("data/nice/JAZ_UV_VIS")

DD_UVB_1314_subsurf_from_Jaz = read.csv("Daily_int_UVB_dose_0.6m_subsurface_PAL1314_kJ_m2.csv",
         stringsAsFactors = FALSE, header = FALSE)
colnames(DD_UVB_1314_subsurf_from_Jaz) = c("Date_raw","Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2")

# date conversion
DD_UVB_1314_subsurf_from_Jaz$Datetime_GMT = strptime(DD_UVB_1314_subsurf_from_Jaz$Date_raw, "%m/%d/%y", tz = "GMT")

# eliminate some bad data

DD_UVB_1314_subsurf_from_Jaz.QA = DD_UVB_1314_subsurf_from_Jaz[-c(40:44),]

# # match relevant data from the two datasets, calculate average UVB xmiss
# 
# DD_UVB_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 4, nrow = nrow(DD_UVB_1314_subsurf)))
# colnames(DD_UVB_xmiss_PAL1314) = c("Date","DD_UVB_incident_kJ_m2","DD_UVB_subsurf_kJ_m2","Percent_xmiss")
# DD_UVB_xmiss_PAL1314$Date = DD_UVB_1314_subsurf$Date
# DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_kJ_m2 = DD_UVB_1314_subsurf$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2
# 
# for (i in 1:nrow(DD_UVB_xmiss_PAL1314)) {
#   
#   if (any(DD_UVB_1314_incident$Date_simple %in% DD_UVB_xmiss_PAL1314$Date[i])) {
#     # there is matching data in the NOAA dataset
#    
#     DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2[i] =
#       DD_UVB_1314_incident$E290.315_kJ_m2[which(DD_UVB_1314_incident$Date_simple %in% DD_UVB_xmiss_PAL1314$Date[i])]
#     
#   }
#   
# }
# 
# DD_UVB_xmiss_PAL1314$Percent_xmiss = DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_kJ_m2/
#   DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2
# 
# # calculate mean, sd % xmiss
# 
# mean(DD_UVB_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# sd(DD_UVB_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# 
# # save the table
# 
# save(DD_UVB_xmiss_PAL1314, file = "Daily_dose_UVB_PAL1314.RData")
# 
# # match relevant data from the two datasets, calculate average UVA xmiss
# 
# DD_UVA_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 4, nrow = nrow(DD_UVA_1314_subsurf)))
# colnames(DD_UVA_xmiss_PAL1314) = c("Date","DD_UVA_incident_kJ_m2","DD_UVA_subsurf_kJ_m2","Percent_xmiss")
# DD_UVA_xmiss_PAL1314$Date = DD_UVA_1314_subsurf$Date
# DD_UVA_xmiss_PAL1314$DD_UVA_subsurf_kJ_m2 = DD_UVA_1314_subsurf$Daily_UVA_dose_0.6m_subsurface_Palmer_kJ_m2
# 
# for (i in 1:nrow(DD_UVA_xmiss_PAL1314)) {
#   
#   if (any(DD_UVA_1314_incident$Date_simple %in% DD_UVA_xmiss_PAL1314$Date[i])) {
#     # there is matching data in the NOAA dataset
#     
#     DD_UVA_xmiss_PAL1314$DD_UVA_incident_kJ_m2[i] =
#       DD_UVA_1314_incident$E315.400_kJ_m2[which(DD_UVA_1314_incident$Date_simple %in% DD_UVA_xmiss_PAL1314$Date[i])]
#     
#   }
#   
# }
# 
# DD_UVA_xmiss_PAL1314$Percent_xmiss = DD_UVA_xmiss_PAL1314$DD_UVA_subsurf_kJ_m2/
#   DD_UVA_xmiss_PAL1314$DD_UVA_incident_kJ_m2
# 
# # calculate mean, sd % xmiss
# 
# mean(DD_UVA_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# sd(DD_UVA_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# 
# # save the table
# 
# save(DD_UVA_xmiss_PAL1314, file = "Daily_dose_UVA_PAL1314.RData")
# 
# # calculate average UVR (UVA + UVB) xmiss
# 
# DD_UVR_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 4, nrow = nrow(DD_UVA_1314_subsurf)))
# colnames(DD_UVR_xmiss_PAL1314) = c("Date","DD_UVR_incident_kJ_m2","DD_UVR_subsurf_kJ_m2","Percent_xmiss")
# DD_UVR_xmiss_PAL1314$Date = DD_UVA_xmiss_PAL1314$Date
# 
# DD_UVR_xmiss_PAL1314$DD_UVR_incident_kJ_m2 =
#   DD_UVA_xmiss_PAL1314$DD_UVA_incident_kJ_m2+
#   DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2
# 
# DD_UVR_xmiss_PAL1314$DD_UVR_subsurf_kJ_m2 =
#   DD_UVA_xmiss_PAL1314$DD_UVA_subsurf_kJ_m2+
#   DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_kJ_m2
# 
# DD_UVR_xmiss_PAL1314$Percent_xmiss =
#   DD_UVR_xmiss_PAL1314$DD_UVR_subsurf_kJ_m2/
#   DD_UVR_xmiss_PAL1314$DD_UVR_incident_kJ_m2
# 
# # calculate mean, sd % xmiss
# 
# mean(DD_UVR_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# sd(DD_UVR_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# 
# # save the table
# 
# save(DD_UVR_xmiss_PAL1314, file = "Daily_dose_UVR_PAL1314.RData")

##### load, process wavelength-specific ("hi-res") data #####
# 
# NOAA data

# data is in separate files, so have to read them in sequentially

setwd(base.wd)
setwd("data/raw/ESRL_GMD_AntUV/ver2") 

# get list of files
ver2_UV_VIS_specfiles = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.[0-9]{3}")

# preallocate destination matrix
# will be pulling the "corrected" spectral data (col 2) from version 2 data files
# see http://uv.biospherical.com/Version2/Description_Spectra.pdf

PAL1314_NOAA_AntUV_spectra_uW_cm2 = as.data.frame(matrix(data = NA, ncol = 622, nrow = length(ver2_UV_VIS_specfiles)))
colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2) = c("Timestamp_GMT",c(seq(280,340,.2),seq(340.5,400,.5),
                                            seq(401,600,1)))

# cycle through files, parse, collect data

for (i in 1:length(ver2_UV_VIS_specfiles)) {
  
  # first, extract base filename containing necessary metadata
  scanFile = str_extract(ver2_UV_VIS_specfiles[i],"BB[0-9]{6}\\.[0-9]{3}$")
  
  # extract metadata from the filename
  scanJuldate = as.numeric(str_extract(scanFile,"[0-9]{3}$")) # julian date
  scanFile = sub("\\.[0-9]{3}$","",scanFile)
  scanMin = as.numeric(str_extract(scanFile,"[0-9]{2}$"))
  scanFile = sub("[0-9]{2}$","",scanFile)
  scanHour = as.numeric(str_extract(scanFile,"[0-9]{2}$"))
  scanFile = sub("[0-9]{2}$","",scanFile)
  scanYear = as.numeric(paste0(20,str_extract(scanFile,"[0-9]{2}$")))
  scanDay = getmoday(scanJuldate,scanYear)[2]
  scanMonth = getmoday(scanJuldate,scanYear)[1]
  
  # create timestamp
  scanTS = ISOdatetime(scanYear, scanMonth, scanDay, scanHour, scanMin, 0, tz = "GMT")
  scanTS.char = as.character(scanTS)
  
  # cheap hack for case where hour and minute both = 0 (i.e., midnight)
  if (scanHour==0 & scanMin==0) {
    
    scanTS.char = sub("$", " 00:00:00",scanTS.char)
    
  }
    
  # extract data from file
  scanData = read.csv(ver2_UV_VIS_specfiles[i], header = FALSE)
  
  # write timestamp, relevant data to destination matrix
  
  PAL1314_NOAA_AntUV_spectra_uW_cm2[i,1] = scanTS.char
  PAL1314_NOAA_AntUV_spectra_uW_cm2[i,2:622] = scanData[,2]
  
}

# convert date to POSIXlt format

PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT = strptime(PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT, "%Y-%m-%d %T", tz = "GMT")

# save the matrix

save(PAL1314_NOAA_AntUV_spectra_uW_cm2, file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Incident_UV-VIS_spectra_PAL1314_uW_cm2.RData"))

# JAZ data

setwd(base.wd)
setwd("data/nice/JAZ_UV_VIS")

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2 = read.csv("JAZ_UV-VIS_full_spectra_0.6m_subsurface_PAL1314_uW_cm2_QA.csv",
                                     stringsAsFactors = FALSE, header = FALSE)

# read in wavelength metadata for this JAZ instrument

JAZ_wavelengths = read.csv("/Users/jamesrco/Code/Optics_Photochem/JAZ_wavelengths.csv",
                           header = FALSE)

# assign column names

colnames(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2) = c("Date_raw_julian",
                                                             "Inttime_microseconds",
                                                             "Badscans_fullspectrum",
                                                             "Badscans_UVB",
                                                             JAZ_wavelengths$V1)

# create timestamp

# convert dates from messed up Excel format
PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT = as.POSIXct(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Date_raw_julian*60*60*24, tz = "GMT", origin = "0000-01-01")-1*60*60*24

# ##### use data from 11-13 Nov 11 (when JAZ was deployed in open air atop tank) for intercalibration with NOAA radiometer #####
# 
# # load in JAZ data
# 
# setwd(base.wd)
# setwd("data/raw/JAZ_UV_VIS")
# 
# # read in data, from .csv output generated in MATLAB (see the MATLAB script
# # JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)
# 
# PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2 = read.csv("JAZ_UV-VIS_open_air_full_spectra_PAL1314_11_13Nov13_uW_cm2.csv",
#                                                                    stringsAsFactors = FALSE, header = FALSE)
# 
# # assign column names (assumes JAZ_wavelenths has already been loaded above)
# 
# colnames(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2) = c("Date_raw_julian",
#                                                                       "Inttime_microseconds",
#                                                                       "Badscans_fullspectrum",
#                                                                       "Badscans_UVB",
#                                                                       JAZ_wavelengths$V1)
# 
# # create timestamp
# 
# # convert dates from messed up Excel format
# PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT = as.POSIXct(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Date_raw_julian*60*60*24, tz = "GMT", origin = "0000-01-01")-1*60*60*24
# 
# # match elements of NOAA data set with JAZ scans, then take a look at correlation
# 
# # subset NOAA data to range of interest
# 
# PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset = 
#   PAL1314_NOAA_AntUV_spectra_uW_cm2[
#     PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT>=
#       as.POSIXct('2013-11-11 12:00:00', tz = "GMT") & 
#       PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT<=
#       as.POSIXct('2013-11-13 12:00:00', tz = "GMT")
#     ,]
# 
# # first, determine which elements to pull from each JAZ data string
# 
# JAZ.match.index.intercal = vector(length = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)-1), mode = "numeric")
# 
# for (i in 1:length(JAZ.match.index.intercal)) {
#   
#   JAZ.match.index.intercal[i] =
#     which((abs(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[i+1])-
#                  JAZ_wavelengths$V1)
#            ==min(abs(JAZ_wavelengths$V1-as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[i+1])))))
#   
# }
# 
# # preallocate matrices to hold data on surface penetration
# 
# PAL1314_intercal= as.data.frame(matrix(data = NA, ncol = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)-1), nrow = nrow(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)))
# colnames(PAL1314_intercal) = colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)]
# 
# # append timestamps
# 
# PAL1314_intercal$Timestamp_GMT = PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[,1]
# 
# # create matrices for scatterplot
# PAL1314_intercal.NOAA = PAL1314_intercal
# PAL1314_intercal.JAZ = PAL1314_intercal
# 
# # perform matching & calculations
# 
# for (i in 1:nrow(PAL1314_intercal)) {
#   
#   # find matching time in the JAZ data, if it exists
#   # we'll want something +/- 30 seconds of the timestamp in the NOAA data
#   
#   possible.time.match =
#     which((abs(PAL1314_intercal$Timestamp_GMT[i]-
#                  PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT)
#            ==min(abs(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT
#                      -PAL1314_intercal$Timestamp_GMT[i]))))
#   
#   if (as.numeric(abs(PAL1314_intercal$Timestamp_GMT[i]-
#           PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT[possible.time.match]))<=
#       15) { # we have a "matching" reading within 30 seconds --> calculate ratio & record
#     
#     PAL1314_intercal[i,1:(ncol(PAL1314_intercal)-1)] = 
#       
#       PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2[possible.time.match,JAZ.match.index.intercal]/
#       PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[i,2:ncol(PAL1314_intercal)]
#     
#     # also, dump data to our two matrices for scatterplot
#     
#     PAL1314_intercal.NOAA[i,1:(ncol(PAL1314_intercal)-1)] =
#       PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[i,2:ncol(PAL1314_intercal)]
#     PAL1314_intercal.JAZ[i,1:(ncol(PAL1314_intercal)-1)] =
#       PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2[possible.time.match,JAZ.match.index.intercal]
#   }
#   
# }
# 
# # calculate mean values for each wavelength, plot
# 
# PAL1314_NOAA.mean = apply(PAL1314_intercal.NOAA[,1:(ncol(PAL1314_intercal.NOAA)-1)],2,mean,na.rm = TRUE)
# PAL1314_NOAA.sd = apply(PAL1314_intercal.NOAA[,1:(ncol(PAL1314_intercal.NOAA)-1)],2,sd,na.rm = TRUE)
# PAL1314_JAZ.mean = apply(PAL1314_intercal.JAZ[,1:(ncol(PAL1314_intercal.JAZ)-1)],2,mean,na.rm = TRUE)
# PAL1314_JAZ.sd = apply(PAL1314_intercal.JAZ[,1:(ncol(PAL1314_intercal.JAZ)-1)],2,sd,na.rm = TRUE)
# plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#      PAL1314_NOAA.mean,
#      xlim = c(290,600),type="l",col="blue")
# # polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
# #         c(PAL1314_NOAA.mean+PAL1314_NOAA.sd,rev(PAL1314_NOAA.mean-PAL1314_NOAA.sd)),
# #         col="lightskyblue")
# # polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
# #         c(PAL1314_JAZ.mean+PAL1314_JAZ.sd,rev(PAL1314_JAZ.mean-PAL1314_JAZ.sd)),
# #         col="lightsalmon")
# lines(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#       PAL1314_NOAA.mean,
#       col="blue")
# lines(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#       PAL1314_JAZ.mean,col="red")
# 
# # plot intercalibration ratio
# 
# PAL1314_intercal.mean = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,mean,na.rm = TRUE)
# plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#      PAL1314_intercal.mean,
#      ylim = c(-3,3),
#      xlim = c(290,600))
# 
# PAL1314_intercal.mean = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,mean,na.rm = TRUE)
# PAL1314_intercal.sd = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,sd,na.rm = TRUE)
# plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#      PAL1314_intercal.mean, type = "p", pch = 16, cex = 0.3, 
#      ylim = c(-10,10),
#      xlim = c(290,600))
# polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
#   c(PAL1314_intercal.mean+PAL1314_intercal.sd,rev(PAL1314_intercal.mean-PAL1314_intercal.sd)),
#   col=grey(0.95))
# points(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
#      PAL1314_intercal.mean, pch = 16, cex = 0.3)
# 
# NOAA.caldat = unlist(PAL1314_intercal.NOAA[1:(ncol(PAL1314_intercal.NOAA)-1)])
# NOAA.caldat = NOAA.caldat[!is.na(NOAA.caldat)]
# NOAA.caldat = as.matrix(NOAA.caldat)
# 
# JAZ.caldat = unlist(PAL1314_intercal.JAZ[1:(ncol(PAL1314_intercal.JAZ)-1)])
# JAZ.caldat = JAZ.caldat[!is.na(JAZ.caldat)]
# JAZ.caldat = as.matrix(JAZ.caldat)
# 
# # write.csv(NOAA.caldat,file="NOAA.caldat.csv")
# # write.csv(JAZ.caldat,file="JAZ.caldat.csv")
# 
# plot(NOAA.caldat/JAZ.caldat,ylim = c(-3,3))
# 
# fit.lm = lm(JAZ.caldat~NOAA.caldat)
# fit.res = resid(fit.lm)
# plot(NOAA.caldat,fit.res)
# 
# fit.lm = lm(JAZ.caldat~NOAA.caldat)
# summary(fit.lm)
# 
# plot(NOAA.caldat,
#      JAZ.caldat,
#      pch=16,cex=0.1)
# abline(0,1,col="red")
# abline(fit.lm, col="blue")

# ##### match elements of NOAA data set with JAZ scans, then calculate surface penetration ##### 
# 
# # will use the NOAA data as our index since they were sampled at more sparse
# # time & wavelength intervals
# 
# # first, determine which elements to pull from each JAZ data string
# 
# JAZ.match.index = vector(length = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)-1), mode = "numeric")
#   
# for (i in 1:length(JAZ.match.index)) {
#   
#   JAZ.match.index[i] =
#     which((abs(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[i+1])-
#               JAZ_wavelengths$V1)
#         ==min(abs(JAZ_wavelengths$V1-as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[i+1])))))
#   
# }
# 
# # preallocate matrix to hold data on surface penetration
# 
# PAL1314_Surf_pen = as.data.frame(matrix(data = NA, ncol = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)-1), nrow = nrow(PAL1314_NOAA_AntUV_spectra_uW_cm2)))
# colnames(PAL1314_Surf_pen) = colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)]
# 
# # append timestamps
# 
# PAL1314_Surf_pen$Timestamp_GMT = PAL1314_NOAA_AntUV_spectra_uW_cm2[,1]
# 
# # perform matching & calculations
# 
# for (i in 1:nrow(PAL1314_Surf_pen)) {
#   
#   # find matching time in the JAZ data, if it exists
#   # we'll want something +/- 30 seconds of the timestamp in the NOAA data
#   
#   possible.time.match =
#     which((abs(PAL1314_Surf_pen$Timestamp_GMT[i]-
#                  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT)
#            ==min(abs(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT
#                      -PAL1314_Surf_pen$Timestamp_GMT[i]))))
#   
#   if (as.numeric(abs(PAL1314_Surf_pen$Timestamp_GMT[i]-
#            PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT[possible.time.match]))<=
#       30*(1/24/60/60)) { # we have a "matching" reading within 30 seconds --> calculate surface penetration (fraction) & record
#     
#     PAL1314_Surf_pen[i,1:(ncol(PAL1314_Surf_pen)-1)] = 
#         
#       PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[possible.time.match,JAZ.match.index]/
#       PAL1314_NOAA_AntUV_spectra_uW_cm2[i,2:ncol(PAL1314_Surf_pen)]
#       
#   }
#   
# }
# 
# # calculate mean values for each wavelength
# 
# PAL1314_Surf_pen.mean = apply(PAL1314_Surf_pen[,1:(ncol(PAL1314_Surf_pen)-1)],2,mean,na.rm = TRUE)
# plot(as.numeric(colnames(PAL1314_Surf_pen)[1:(ncol(PAL1314_Surf_pen)-1)]),
#      PAL1314_Surf_pen.mean,
#      ylim = c(0,3),
#      xlim = c(290,600))

##### Kd  calc & other workup using static depth UV-VIS profile from 2015-2016 season ##### 

# load in JAZ profile data (made on 15 Dec 2015 at Station B, Arthur Harbor)

setwd(base.wd)
setwd("data/raw/JAZ_UV_VIS")

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2 = read.csv("JAZ_UV-VIS_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.csv", 
                                                          stringsAsFactors = FALSE, header = FALSE)

# assign column names (assumes JAZ_wavelenths has already been loaded above)

colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2) = c("Scan_file_ID",
                                                             "Date_raw_julian",
                                                             "Inttime_microseconds",
                                                             "Badscans_fullspectrum",
                                                             "Badscans_UVB",
                                                             JAZ_wavelengths$V1)

# create timestamp

# convert dates from messed up Excel format
PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2$Timestamp_GMT = as.POSIXct(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2$Date_raw_julian*60*60*24, tz = "GMT", origin = "0000-01-01")-1*60*60*24

# assign depths to the scans (from PAL1516 field notebook, p. 24)

PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2$Depth_m =
  c(8,7.5,7,6.5,6,6,5.5,5,4.5,4,3.5,3,2.5,2,1.5,1,0.5,NaN,0,NaN)

# load in surface PAR data taken concurrently with JAZ profile

setwd(base.wd)
setwd("data/raw/LI_193_PAR")

PAL1516_Stn_B_surface_PAR_20151215_umol_s_m2 = read.csv("PAL1516_Surface_PAR_for_correlation_15Dec15.csv", 
                                                                   stringsAsFactors = FALSE, header = TRUE, skip = 3)

# calculate some surface reference correction factors by which JAZ data can be corrected
# requires PAR data be detrended due to hysteresis

PAR.lm = lm(PAL1516_Stn_B_surface_PAR_20151215_umol_s_m2$PAR.at.surface.from.LI.193SA..umol.s.1.m.2. ~ c(1:17))

JAZ.corrfactors = as.matrix((PAL1516_Stn_B_surface_PAR_20151215_umol_s_m2$PAR.at.surface.from.LI.193SA..umol.s.1.m.2.-fitted(PAR.lm)+PAL1516_Stn_B_surface_PAR_20151215_umol_s_m2$PAR.at.surface.from.LI.193SA..umol.s.1.m.2.[1])/PAL1516_Stn_B_surface_PAR_20151215_umol_s_m2$PAR.at.surface.from.LI.193SA..umol.s.1.m.2.[1])

# surface-corrected JAZ data
# also, eliminate some unwanted scan data

PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected =
  PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2[-c(5,18,20),]

PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[,6:(ncol(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)-2)] = 
  PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[,6:(ncol(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)-2)]*JAZ.corrfactors

# plots

# full spectrum 315-600 nm, with room for interpolated Kd data from 290-315

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

color_ramp.depths = colorRampPalette(c("cadetblue1", "darkblue")) # create a color ramp for shading the points according to lipid class

depth.plot.colors = rev(color_ramp.depths(nrow(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)))

plot(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[338:1141]),
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[1,338:1141],"l",
     col = depth.plot.colors[1], lty = 1, lwd = "1.5",
     ylim = c(0,25), xlim = c(290,600),
     ylab = expression(paste("Irradiance (",mu,"W ",cm^-2," ",nm^-1,")")),
     xlab = "Wavelength (nm)")

text(589,
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[1,1108],
     labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[1]),
     offset = -.25, pos = 3)

# overlay other depths

for (i in c(3,5,7,9,11,13,15,16,17)) {
  
  lines(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[338:1141]),
        PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,338:1141],
        col = depth.plot.colors[i], lty = 1, lwd = "1.5")
  
  if (i %in% c(7,9,11,13,15,16,17)) {
  text(589,
       PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,1108],
       labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[i]),
       offset = -.25, pos = 3)
    
  }
  
}

# # overlay concurrent reading from NOAA radiometer at Palmer, if desired
# # doesnt match up exactly 
# 
# PAL1314_NOAA_AntUV_spectra_uW_cm2.sub =
#   PAL1314_NOAA_AntUV_spectra_uW_cm2[which(abs(PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT-
#                                                 as.POSIXct('2013-12-15 11:11:00', tz = "GMT")) < 10),]
# 
# lines(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)]),
#      PAL1314_NOAA_AntUV_spectra_uW_cm2.sub[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)],
#      col="black")

dev.off()
# 
# # UVB-range inset
# 
# # make plot for inset (290-315 nm)
# 
# par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
# 
# pdf(file = "AH_depth_profile_20151215_inset.pdf",
#     width = 4, height = 4, pointsize = 12,
#     bg = "white")
# 
# plot(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
#      PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,270:1500],"l",
#      col = depth.plot.colors[1], lty = 1, lwd = "1.5",
#      ylim = c(0,5), xlim = c(290,315),
#      ylab = expression(paste("Irradiance (",mu,"W ",cm^-2," ",nm^-1,")")),
#      xlab = "Wavelength (nm)")
# 
# text(315,
#      PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,338],
#      labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[3]),
#      offset = -.5, pos = 3)
# 
# # overlay other depths
# 
# for (i in c(7,11,15,16,17)) {
#   
#   lines(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
#         PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,270:1500],
#         col = depth.plot.colors[i], lty = 1, lwd = "1.5")
#   
#   if (i %in% c(7,9,11,13,15,16,17)) {
#     text(315,
#          PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,338],
#          labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[i]),
#          offset = -.5, pos = 3)
#     
#   }
#   
# }
# 
# dev.off()

# now, calculate Kds for wavelengths 320-700 nm (Jaz data < 320 nm are too noisy)
# 320 nm is 346th element in the Jaz wavelengths vector; 700 nm is 1437th

# preallocate matrix to hold data

PAL1516_AH_Kd_20151215 = as.data.frame(matrix(data = NA, ncol = nrow(JAZ_wavelengths), nrow = 2))
colnames(PAL1516_AH_Kd_20151215) = JAZ_wavelengths$V1

# use 0.5 m (assumed to be the subsurface incident irradiance, E(0))
# and 7 m depth irradiance data (as E(z)) to calculate Kd with z = 7

for (i in 346:1437) {
  
PAL1516_AH_Kd_20151215[1,i] = 
  
  (log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,i+5])-
     log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[16,i+5]))/-7

}

# clean up, transform result 
PAL1516_AH_Kd_20151215_per_meter = as.data.frame(t(PAL1516_AH_Kd_20151215))
colnames(PAL1516_AH_Kd_20151215_per_meter) = c("Kd_per_meter","Kd_per_meter_fitted")

# set any of these calculated Kds with value > 1 = NA
PAL1516_AH_Kd_20151215_per_meter[PAL1516_AH_Kd_20151215_per_meter>=1] = NA

# now, fit an exponential model to estimate Kds for wavelengths 290-320 nm
# will use data from 320-370 nm as basis for fitting curve, then back-extrapolate

Kd_fit_subset = log(PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter[346:482])
Wavelength_fit_subset = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[346:482])

fit.exp = lm(Kd_fit_subset ~ Wavelength_fit_subset)

# get some information
summary(fit.exp)

# predict values for 290-320 nm using model we just created

UVB_pred_subset = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[264:346])
Kd_UVB.pred = exp(predict(fit.exp,list(Wavelength_fit_subset=UVB_pred_subset)))

# create vector of fitted values for range of data we used to fit the curve routine

Kd_fitrange_wavelengths = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[346:482])
Kd_fitrange.fitted = exp(predict(fit.exp,list(Wavelength_fit_subset=Kd_fitrange_wavelengths)))

# quick plot to see if this all makes sense

plot(JAZ_wavelengths$V1,
     PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter,"l",
     col = "black", lty = 1, lwd = "1.5",
     ylim = c(0.05,0.9), xlim = c(290,600),
     ylab = expression(paste(K_d)),
     xlab = "Wavelength (nm)")
lines(UVB_pred_subset,Kd_UVB.pred,lty=3)
lines(Kd_fitrange_wavelengths,Kd_fitrange.fitted,lty=3)

# save corrected irradiance data 

save(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected, 
     file = paste0(base.wd,"/data/nice/JAZ_UV_VIS/JAZ_UV_VIS_corrected_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.RData"))
write.csv(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected, 
          file = paste0(base.wd,"/data/nice/JAZ_UV_VIS/JAZ_UV_VIS_corrected_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.csv"))

# append model-predicted Kd values; save Kds

PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter_fitted[as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)) %in%
                                   Kd_fitrange_wavelengths]=Kd_fitrange.fitted
PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter_fitted[as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)) %in%
                                   UVB_pred_subset]=Kd_UVB.pred

save(PAL1516_AH_Kd_20151215_per_meter, 
     file = paste0(base.wd,"/data/nice/derived_light_calculations/Kd_PAL1516_Arthur_Hbr_Stn_B_20151215.RData"))
write.csv(PAL1516_AH_Kd_20151215_per_meter, 
          file = paste0(base.wd,"/data/nice/derived_light_calculations/Kd_PAL1516_Arthur_Hbr_Stn_B_20151215.csv"))

# a better plot for the manuscript

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215_Kd.pdf",
    width = 6, height = 3, pointsize = 12,
    bg = "white")

# Kds derived directly from observationss
plot(JAZ_wavelengths$V1[483:nrow(JAZ_wavelengths)],
     PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter[483:nrow(PAL1516_AH_Kd_20151215_per_meter)],
     "l",
     col = "black", lty = 1, lwd = "1",
     ylim = c(0.05,0.9), xlim = c(290,600),
     ylab = expression(paste(K_d)),
     xlab = "Wavelength (nm)")

# grey out but retain the range used for curve fitting
lines(JAZ_wavelengths$V1[346:482],
      PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter[346:482],
      col = "lightgrey")
      
# fitted values
lines(UVB_pred_subset,Kd_UVB.pred,lty=3)
lines(Kd_fitrange_wavelengths,Kd_fitrange.fitted,lty=3)

dev.off()

# # just UV wavelengths
# 
# par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
# 
# pdf(file = "AH_depth_profile_20151215_Kd_UV.pdf",
#     width = 4, height = 3, pointsize = 12,
#     bg = "white")
# 
# plot(JAZ_wavelengths$V1,
#      PAL1516_AH_Kd_20151215,"l",
#      col = "black", lty = 1, lwd = "1.5",
#      ylim = c(0.05,0.9), xlim = c(290,315),
#      ylab = expression(K_d),
#      xlab = "Wavelength (nm)")
# 
# dev.off()

##### use newly calculated/estimated Kds to calculate estimated irradiances at depth #####

# first, use Kds to estimate some fluxes at depths from 0-10 m using incident data
# in PAL1314_NOAA_AntUV_spectra_uW_cm2 (1 m intervals)

# define depths, interval for calculation

maxdepth = 10
depthint = 1

depths = seq(from = 0, to = maxdepth, by = 1)

# create a list object to hold the data

PAL1314_est_spectra_at_depth_uW_cm2 =
  vector("list",length(depths))

names(PAL1314_est_spectra_at_depth_uW_cm2) =
  depths

# populate first slot with incident (0 m) data

PAL1314_est_spectra_at_depth_uW_cm2[[1]] =
 as.matrix(PAL1314_NOAA_AntUV_spectra_uW_cm2[,2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)])

# define a function to calculate the estimated irradiance at depth z (in m) and
# wavelength lambda (nm), given incident irradiance incIrrad

# assumes you have a data frame of Napierian Kds in format identical to
# PAL1516_AH_Kd_20151215_per_meter, above

estIrrad = function(incIrrad, lambda, z, Kds) {
  
  # retrieve "best" Kd for this wavelength
  
    # find closest wavelength
  
    Kd_ind = which.min(abs(as.numeric(rownames(Kds))-lambda))
  
    # use Kd derived directly from measurements if it exists; otherwise, use a
    # fitted value
  
    if (!is.na(Kds[Kd_ind,1])) {
      
      Kd = Kds[Kd_ind,1]
      
    } else {
      
      Kd = Kds[Kd_ind,2]
      
    }
    
  # calculate the estimated irradiance at depth z
    
  estIrrad_z = incIrrad*exp(-Kd*z)
  
  # return our estimate irradiance

  return(estIrrad_z)
  
}
    
# make calculations for 1-10 m (or sequence), by wavelength; store in PAL1314_est_spectra_at_depth_uW_cm2

for (i in 2:length(depths)) {
  PAL1314_est_spectra_at_depth_uW_cm2[[i]] = mapply(estIrrad, lambda = as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)]),
              incIrrad = PAL1314_NOAA_AntUV_spectra_uW_cm2[,2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)],
              z = depths[i],
              MoreArgs = list(Kds = PAL1516_AH_Kd_20151215_per_meter))
  
  # change any irradiances < 0 to 0
  
  PAL1314_est_spectra_at_depth_uW_cm2[[i]][PAL1314_est_spectra_at_depth_uW_cm2[[i]]<0] = 0
  
  # set colnames, rownames
  
  colnames(PAL1314_est_spectra_at_depth_uW_cm2[[i]]) = colnames(PAL1314_est_spectra_at_depth_uW_cm2[[1]])
  
}

# also, create a single matrix of estimated irradiances at 0.6 m depth

PAL1314_est_spectra_at_0.6m_uW_cm2 = mapply(estIrrad, lambda = as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)]),
                                                  incIrrad = PAL1314_NOAA_AntUV_spectra_uW_cm2[,2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)],
                                                  z = 0.6,
                                                  MoreArgs = list(Kds = PAL1516_AH_Kd_20151215_per_meter))

# change any irradiances < 0 to 0

PAL1314_est_spectra_at_0.6m_uW_cm2[PAL1314_est_spectra_at_0.6m_uW_cm2<0] = 0

##### calculate daily spectral doses at various depths using our Kd-adjusted irradiances ##### 

# previously used in situ Jaz data to obtain these integrals, but want to use
# NOAA data for sake of consistency since that's where our incident daily
# data come from

# some relevant information for calculation of daily integrals, from NOAA AntUV documentation:

# "Time on daily doses refers to approximate local
# solar noon (01:00 for McMurdo Station; 16:00 for Palmer Station; 12:00 for South Pole
# Station; 17:00 for Ushuaia; 20:00 for San Diego, 21:00 for Barrow, and 15:00 for
# Summit). Integration boundaries for daily dose calculations are the given time +/- 12
# hours. For example: The daily dose for McMurdo assigned to the Date/Time stamp
# 36826.04167 (i.e. 10/27/00 01:00 GMT) is the integral of spectral irradiance between
# 10/26/00 13:00 GMT and 10/27/00 13:00 GMT."

# preallocate objects to hold data
# will make calculations from 15 Oct - 30 Dec 2013
# can use objects holding incident data as starting point

DD_UVB_1314_est_at_0.6m_kJ_m2 = DD_UVB_1314_incident[
  DD_UVB_1314_incident$Date >=
    as.POSIXct('2013-10-15 00:00:00', tz = "GMT") &
    DD_UVB_1314_incident$Date <=
    as.POSIXct('2013-12-30 23:59:59', tz = "GMT")
  ,c(1,3:4)]
DD_UVB_1314_est_at_0.6m_kJ_m2$DD_kJ_m2 = NA

DD_UVA_1314_est_at_0.6m_kJ_m2 = DD_UVA_1314_incident[
  DD_UVA_1314_incident$Date >=
    as.POSIXct('2013-10-15 00:00:00', tz = "GMT") &
    DD_UVA_1314_incident$Date <=
    as.POSIXct('2013-12-30 23:59:59', tz = "GMT")
  ,c(1,3:4)]
DD_UVA_1314_est_at_0.6m_kJ_m2$DD_kJ_m2 = NA

# get vector of wavelengths in NOAA data

NOAA_AntUV_lambdas = 
  as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)])

# now, can calculate some daily integrals at 0.6 m

for (i in 1:nrow(DD_UVB_1314_est_at_0.6m_kJ_m2)) { # iterate through dates
  
  # pull out relevant data for this date, bearing in mind the +/- 12 hr 
  # convention described above

  PAL1314_est_spectra.today =  
  PAL1314_est_spectra_at_0.6m_uW_cm2[
   PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT>
    DD_UVB_1314_est_at_0.6m_kJ_m2$Date[i]-60*60*12 & 
    PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT<=
     DD_UVB_1314_est_at_0.6m_kJ_m2$Date[i]+60*60*12,]
  
  Todays.timestamps = 
    PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT[
      PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT>
        DD_UVB_1314_est_at_0.6m_kJ_m2$Date[i]-60*60*12 & 
        PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT<=
        DD_UVB_1314_est_at_0.6m_kJ_m2$Date[i]+60*60*12]
  
  # create a vector of times in seconds
  
  Timevector.s = Todays.timestamps-(DD_UVB_1314_est_at_0.6m_kJ_m2$Date[i]-60*60*12)
  
  units(Timevector.s) = "secs"
  
  Timevector.s = as.numeric(Timevector.s)
  
  # preallocate vector for spectral integrals
  
  Wavelength_integrals_UVB = 
    vector(mode = "numeric", length(Todays.timestamps))
  
  Wavelength_integrals_UVA = 
    vector(mode = "numeric", length(Todays.timestamps))
  
  for (j in 1:length(Wavelength_integrals_UVB)) { # iterate through tinme, integrate by wavelength
    
    Wavelength_integrals_UVB[j] = caTools::trapz(NOAA_AntUV_lambdas[NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<315],
                     PAL1314_est_spectra.today[j,NOAA_AntUV_lambdas>=290 &
                                                   NOAA_AntUV_lambdas<315])
    
    Wavelength_integrals_UVA[j] = caTools::trapz(NOAA_AntUV_lambdas[NOAA_AntUV_lambdas>=315 & NOAA_AntUV_lambdas<400],
                                                 PAL1314_est_spectra.today[j,NOAA_AntUV_lambdas>=315 &
                                                                             NOAA_AntUV_lambdas<400])
    
  }
  
  # integrate by time; scale; store
  
  DD_UVB_1314_est_at_0.6m_kJ_m2$DD_kJ_m2[i] = 
    
    caTools::trapz(Timevector.s,
                   Wavelength_integrals_UVB*W_per_uW*J_per_kJ*cm2_per_m2)
    
  DD_UVA_1314_est_at_0.6m_kJ_m2$DD_kJ_m2[i] = 
    
    caTools::trapz(Timevector.s,
                   Wavelength_integrals_UVA*W_per_uW*J_per_kJ*cm2_per_m2)
  
}

# collect data and save

# estimates at 0.6 m

DD_1314_est_at_0.6m_kJ_m2 = DD_UVB_1314_est_at_0.6m_kJ_m2
colnames(DD_1314_est_at_0.6m_kJ_m2)[c(2,4)] = c("Date_GMT","DD_UVB_est_at_0.6m_kJ_m2")
DD_1314_est_at_0.6m_kJ_m2$DD_UVA_est_at_0.6m_kJ_m2 = DD_UVA_1314_est_at_0.6m_kJ_m2$DD_kJ_m2

save(DD_1314_est_at_0.6m_kJ_m2, 
     file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_doses_0.6m_from_Kds_and_NOAA_data_PAL1314.RData"))
write.csv(DD_1314_est_at_0.6m_kJ_m2, 
          file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Daily_int_doses_0.6m_from_Kds_and_NOAA_data_PAL1314.csv"))

# estimates for other depths

# append times, wavelengths

PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT = PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT
PAL1314_est_spectra_at_depth_uW_cm2$lambda_nm = NOAA_AntUV_lambdas

save(PAL1314_est_spectra_at_depth_uW_cm2, 
     file = paste0(base.wd,"/data/nice/NOAA_ESRL_GMD_AntUV/Est_UV_VIS_spectra_at_depth_from_Kds_and_NOAA_data_PAL1314_uW_cm2.RData"))

# calculate % transmission to 0.6 m for the various integrals

# UVA

# match relevant data from the two datasets, calculate average UVA xmiss

DD_UVA_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 4, nrow = nrow(DD_UVA_1314_est_at_0.6m_kJ_m2)))
colnames(DD_UVA_xmiss_PAL1314) = c("Datetime_GMT","DD_UVA_incident_kJ_m2","DD_UVA_subsurf_from_Kd_kJ_m2","Percent_xmiss")
DD_UVA_xmiss_PAL1314$Datetime_GMT = DD_UVA_1314_est_at_0.6m_kJ_m2$Date
DD_UVA_xmiss_PAL1314$DD_UVA_subsurf_from_Kd_kJ_m2 = DD_UVA_1314_est_at_0.6m_kJ_m2$DD_kJ_m2

for (i in 1:nrow(DD_UVA_xmiss_PAL1314)) {

  if (any(DD_UVA_1314_incident$Date %in% DD_UVA_xmiss_PAL1314$Datetime_GMT[i])) {
    # there is matching data in the incident NOAA data

    DD_UVA_xmiss_PAL1314$DD_UVA_incident_kJ_m2[i] =
      DD_UVA_1314_incident$E315.400_kJ_m2[which(DD_UVA_1314_incident$Date %in% DD_UVA_xmiss_PAL1314$Datetime_GMT[i])]

  }

}

DD_UVA_xmiss_PAL1314$Percent_xmiss = DD_UVA_xmiss_PAL1314$DD_UVA_subsurf_from_Kd_kJ_m2/
  DD_UVA_xmiss_PAL1314$DD_UVA_incident_kJ_m2

# calculate mean, sd % xmiss

mean(DD_UVA_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
sd(DD_UVA_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
# 
# # save the table
# 
# save(DD_UVA_xmiss_PAL1314, file = "Daily_dose_UVA_PAL1314.RData")

# UVB

# match relevant data from the three datasets, calculate average UVB xmiss

DD_UVB_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 6, nrow = nrow(DD_UVB_1314_est_at_0.6m_kJ_m2)))
colnames(DD_UVB_xmiss_PAL1314) = c("Datetime_GMT","DD_UVB_incident_kJ_m2","DD_UVB_subsurf_from_Kd_kJ_m2","DD_UVB_subsurf_from_Jaz_kJ_m2","Percent_xmiss_NOAA","Percent_xmiss_Jaz")

DD_UVB_xmiss_PAL1314$Datetime_GMT = DD_UVB_1314_est_at_0.6m_kJ_m2$Date
DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Kd_kJ_m2 = DD_UVB_1314_est_at_0.6m_kJ_m2$DD_kJ_m2

for (i in 1:nrow(DD_UVB_xmiss_PAL1314)) {
  
  if (any(DD_UVB_1314_incident$Date %in% DD_UVB_xmiss_PAL1314$Datetime_GMT[i])) {
    # there is matching data in the incident NOAA data
    
    DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2[i] =
      DD_UVB_1314_incident$E290.315_kJ_m2[which(DD_UVB_1314_incident$Date %in% DD_UVB_xmiss_PAL1314$Datetime_GMT[i])]
    
  }
  
}

for (i in 1:nrow(DD_UVB_xmiss_PAL1314)) {
  
  if (any(as.POSIXct(strptime(DD_UVB_1314_subsurf_from_Jaz.QA$Datetime_GMT, format = "%Y-%m-%d"), tz = "GMT") %in% as.POSIXct(strptime(DD_UVB_xmiss_PAL1314$Datetime_GMT, format = "%Y-%m-%d"), tz = "GMT")[i])) {
    # there is matching data in the Jaz data
    
    DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Jaz_kJ_m2[i] =
      DD_UVB_1314_subsurf_from_Jaz.QA$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2[which(as.POSIXct(strptime(DD_UVB_1314_subsurf_from_Jaz.QA$Datetime_GMT, format = "%Y-%m-%d"), tz = "GMT") %in% as.POSIXct(strptime(DD_UVB_xmiss_PAL1314$Datetime_GMT, format = "%Y-%m-%d"), tz = "GMT")[i])]
    
  }
  
}

# remove some NaNs
DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Jaz_kJ_m2[is.nan(DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Jaz_kJ_m2)] = NA

DD_UVB_xmiss_PAL1314$Percent_xmiss_NOAA = DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Kd_kJ_m2/
  DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2

DD_UVB_xmiss_PAL1314$Percent_xmiss_Jaz = DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_from_Jaz_kJ_m2/
  DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2

# calculate mean, sd % xmiss

mean(DD_UVB_xmiss_PAL1314$Percent_xmiss_NOAA, na.rm = TRUE)
sd(DD_UVB_xmiss_PAL1314$Percent_xmiss_NOAA, na.rm = TRUE)

mean(DD_UVB_xmiss_PAL1314$Percent_xmiss_Jaz, na.rm = TRUE)
sd(DD_UVB_xmiss_PAL1314$Percent_xmiss_Jaz, na.rm = TRUE)

#####  other panels for Experiment 13 plot ##### 

# instantaneous/cumulative dosages

# subset to ranges of interest

# JAZ

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT>=
      as.POSIXct('2013-12-14 12:32:00', tz = "GMT") &  # 9:30 local time
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT<=
      as.POSIXct('2013-12-14 20:50:00', tz = "GMT") # 17:50 local time
      ,]

# UVB

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec[,
        c(1:4,
          4+which(JAZ_wavelengths>=290 & JAZ_wavelengths<=315),
          ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec))]


# UVR

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec[,
          c(1:4,
            4+which(JAZ_wavelengths>=290 & JAZ_wavelengths<=395.5),
          ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec))]


# NOAA data

PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec = 
  PAL1314_NOAA_AntUV_spectra_uW_cm2[
    PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT>=
      as.POSIXct('2013-12-14 12:32:00', tz = "GMT") &  # 9:30 local time
      PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT<=
      as.POSIXct('2013-12-14 20:50:00', tz = "GMT") # 17:50 local time
    ,]

PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.UVB = 
  PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec[,which(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec))<=315)]

# instantaneous doses

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB = 
  apply(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB)-1)],1,sum)

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVR = 
  apply(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR)-1)],1,sum)

PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB = 
  apply(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.UVB,1,sum)

# at 0.6 m, from JAZ

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT,
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB,col="red",
  ylim = c(0,500)) 

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVR,col="red",
     ylim = c(0,5000)) 

# superimpose points from NOAA data
points(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec$Timestamp_GMT,
       PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB,col="blue")

# out of curiosity, take a look at average ratio

mean.PAL1314.14Dec.JAZ.instUVB = 
  mean(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB)

mean.PAL1314.14Dec.NOAA.instUVB = 
  mean(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB)

# calculate cumulative dosages

# need library caTools, first have to detach RSEIS since there's a conflict with a function name

detach("package:RSEIS", unload=TRUE)
library(caTools)

# vector of times in seconds, with t = 0 being timepoint at beginning of experiment

Exp13_timeint.s = as.numeric(rev(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT[nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB)]-
                                   PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT))

# UVB

# preallocate vector

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB =
  vector(length = length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB))

for (i in 1:length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB[i] =
    
    caTools::trapz(Exp13_timeint.s[1:i],PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB[1:i])

}

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB.kJ_m2 =
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB*W_per_uW*J_per_kJ*cm2_per_m2

# UVR

# preallocate vector

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR =
  vector(length = length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVR))

for (i in 1:length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR[i] =
    
    caTools::trapz(Exp13_timeint.s[1:i],PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVR[1:i])
  
}

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR.kJ_m2 =
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR*W_per_uW*J_per_kJ*cm2_per_m2

# plot of instaneous and cumulative UVB dosage at 0.6 m

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "14Dec_UVB_timeseries.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,5))

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB,
     col = "black",
     type = "l",
     ylab = expression(paste("Instantaneous UVB irradiance at 0.6 m (",mu,"W ",cm^-2," ",s^-1,")")),
     xlab = "Time (GMT)"
     )

par(new = T)

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB.kJ_m2,
     col = "black",
     type = "l", lty = 2,
     axes = FALSE,
     ylab = "",
     xlab = "")

axis(side = 4)
mtext(side = 4,
      line = 3,
      expression(paste("Cumulative UVB dose received at 0.6 m (kJ ",m^-2,")")))

dev.off()

# plot of instaneous and cumulative UVR dosage at 0.6 m

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "14Dec_UVR_timeseries.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,5))

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVR,
     col = "black",
     type = "l",
     ylab = expression(paste("Instantaneous UVR irradiance at 0.6 m (",mu,"W ",cm^-2," ",s^-1,")")),
     xlab = "Time (GMT)"
)

par(new = T)

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVR.kJ_m2,
     col = "black",
     type = "l", lty = 2,
     axes = FALSE,
     ylab = "",
     xlab = "")

axis(side = 4)
mtext(side = 4,
      line = 3,
      expression(paste("Cumulative UVR dose received at 0.6 m (kJ ",m^-2,")")))

dev.off()

# alternative plots in photon units

# *** requires data frames lambda_nm_JAZ and 
# PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub from
# PAL1314_AQY_calc.R

# subset to just UVR wavelengths

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR =
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub[,
       which(JAZ_wavelengths>=290 & JAZ_wavelengths<=395.5)]

# instantaneous dosages

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR = 
  apply(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR,1,sum)

# cumulative

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR =
  vector(length = length(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR))

for (i in 1:length(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR[i] =
    
    caTools::trapz(Exp13_timeint.s[1:i],PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR[1:i])
  
}

# a quick check of math, comparing to a number from PAL1314_AQY_calc.R which should be same as last value in
# PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR

sum(PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2[which(JAZ_wavelengths>=290 & JAZ_wavelengths<=395.5)])

# also calculating cumulative, adjusted for reduction in transmissivity in quartz vials
# *** requires "FracTrans" and some additional objects from PAL1314_AQY_calc.R

# preallocate matrix, get list of UVR lambdas from JAZ

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR.quartz =
  matrix(NA, nrow = nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR),
         ncol = ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR))

UVR.lambdas =
  JAZ_wavelengths[which(JAZ_wavelengths>=290 & JAZ_wavelengths<=395.5),1]

# adjust values

for (i in 1:ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR.quartz)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR.quartz[,i] =
    PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR[,i]*
    FracTrans$transmittance_quartz_pct[abs(FracTrans$lambda_nm-UVR.lambdas[i])==min(abs(FracTrans$lambda_nm-UVR.lambdas[i]))]
  
}

# instantaneous dosages, adjusted for quartz transmissivity

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR.quartz = 
  apply(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.UVR.quartz,1,sum)

# cumulative, adjusted for quartz transmissivity

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR.quartz =
  vector(length = length(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR.quartz))

for (i in 1:length(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR.quartz)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR.quartz[i] =
    
    caTools::trapz(Exp13_timeint.s[1:i],PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR.quartz[1:i])
  
}

# now, the plot...

# plot of instaneous and cumulative UVR dosage at 0.6 m, in photon units

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "14Dec_UVR_timeseries_umol_photons_m2.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,5))

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub.totalinstUVR,
     col = "black",
     type = "l",
     ylab = expression(paste("Instantaneous UVR dose at 0.6 m (",mu,"mol photons",m^-2," ",s^-1,")")),
     xlab = "Time (GMT)"
)

par(new = T)

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVR$Timestamp_GMT,
     PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2.14Dec.cumUVR,
     col = "black",
     type = "l", lty = 2,
     axes = FALSE,
     ylab = "",
     xlab = "")

axis(side = 4)
mtext(side = 4,
      line = 3,
      expression(paste("Cumulative UVR dose received at 0.6 m (umol photons ",m^-2,")")))

dev.off()


##### pull out JAZ data subset for 2 December 2013 (Exp_03a) #####

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT>=
      as.POSIXct('2013-12-02 11:40:00', tz = "GMT") &  # 8:40 local time
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT<=
      as.POSIXct('2013-12-02 21:28:00', tz = "GMT") # 18:28 local time
    ,]

# ... and create a vector of times in seconds, with t = 0 being timepoint at beginning of experiment

Exp03a_timeint.s = as.numeric(rev(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec$Timestamp_GMT[nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec)]-
                                   PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec$Timestamp_GMT))
