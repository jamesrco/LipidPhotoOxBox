# UV_TS_analysis_PAL1314.R
#
# Purpose: Analysis & treatment of UV timeseries data from Palmer 2013-2014
# field season. Uses NOAA ESRL GMD Antarctic UV irradiance data and data from
# Jaz UV-VIS spectrophotometer J.R.C. deployed at Palmer Station during the 
# 2013-2014 field season
#
# Also, some calculations of the diffuse attenuation coefficient, Kd, from
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

# assuming initial wd is that of the LipidPhotoOxBox repository, get current
# working directory so we can reset it later
initial.wd = getwd()

# # or, this:
# setwd("/Users/jrcollins/Code/LipidPhotoOxBox")

##### daily UVB dosages ##### 

# NOAA data

# set wd to data location
setwd("data/raw/ESRL_GMD_AntUV/ver2") 

# read in NOAA ESRL AntUV data (Version 2, daily dosage data), subset as needed
PAL_AntUV_DD_v2 = read.csv("PAL_all_DailyDose_no_flags_SZA_max_set_to82.csv", 
                  stringsAsFactors = FALSE, skip = 2)

# just the 290-315 nm wavelength dosage data
DD_UVB_incident = as.data.frame(cbind(PAL_AntUV_DD_v2$Date,PAL_AntUV_DD_v2$E290.315))
colnames(DD_UVB_incident) = c("Date_raw_julian","E290.315_kJ_m2")

# convert dates from messed up Excel format
DD_UVB_incident$Date = as.POSIXct('1899-12-30')+(DD_UVB_incident$Date_raw_julian*24*60*60)
DD_UVB_incident$Date_simple = strptime(DD_UVB_incident$Date,"%Y-%m-%d", tz = "GMT")

# subset, plot 2013-2014 data
DD_UVB_1314_incident = DD_UVB_incident[DD_UVB_incident$Date>as.POSIXct('2013-06-22') &
                                         DD_UVB_incident$Date<as.POSIXct('2014-06-22'),]

# plot(DD_UVB_1314_incident$Date,DD_UVB_1314_incident$E290.315_kJ_m2, "l", lwd = 1)
# 
# # save the data file
# save(DD_UVB_1314_incident, file = "NOAA_ESRL_AntUV_DD_UVB_1314.RData")
# 
# # export to .csv
# write.csv(DD_UVB_1314_incident, file = "NOAA_ESRL_AntUV_DD_UVB_1314.csv")

# JAZ data

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

setwd(initial.wd)
setwd("data/nice/JAZ_UV_VIS") 

DD_UVB_1314_subsurf = read.csv("Daily_int_UVB_dose_0.6m_subsurface_PAL1314_kJ_m2.csv", 
         stringsAsFactors = FALSE, header = FALSE)
colnames(DD_UVB_1314_subsurf) = c("Date_raw","Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2")

# date conversion
DD_UVB_1314_subsurf$Date = strptime(DD_UVB_1314_subsurf$Date_raw, "%m/%d/%y", tz = "GMT")

# match relevant data from the two datasets, calculate average UVB xmiss

DD_UVB_xmiss_PAL1314 = as.data.frame(matrix(data = NA, ncol = 4, nrow = nrow(DD_UVB_1314_subsurf)))
colnames(DD_UVB_xmiss_PAL1314) = c("Date","DD_UVB_incident_kJ_m2","DD_UVB_subsurf_kJ_m2","Percent_xmiss")
DD_UVB_xmiss_PAL1314$Date = DD_UVB_1314_subsurf$Date
DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_kJ_m2 = DD_UVB_1314_subsurf$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2

for (i in 1:nrow(DD_UVB_xmiss_PAL1314)) {
  
  if (any(DD_UVB_1314_incident$Date_simple %in% DD_UVB_xmiss_PAL1314$Date[i])) {
    # there is matching data in the NOAA dataset
   
    DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2[i] =
      DD_UVB_1314_incident$E290.315_kJ_m2[which(DD_UVB_1314_incident$Date_simple %in% DD_UVB_xmiss_PAL1314$Date[i])]
    
  }
  
}

DD_UVB_xmiss_PAL1314$Percent_xmiss = DD_UVB_xmiss_PAL1314$DD_UVB_subsurf_kJ_m2/
  DD_UVB_xmiss_PAL1314$DD_UVB_incident_kJ_m2

# calculate mean, sd % xmiss

mean(DD_UVB_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)
sd(DD_UVB_xmiss_PAL1314$Percent_xmiss, na.rm = TRUE)

# save the table

save(DD_UVB_xmiss_PAL1314, file = "Daily_dose_UVB_PAL1314.RData")

##### load, process wavelength-specific ("hi-res") data ##### 

# JAZ data

setwd(initial.wd)
setwd("data/nice/JAZ_UV_VIS") 

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2 = read.csv("JAZ_UV-VIS_full_spectra_0.6m_subsurface_PAL1314_uW_cm2_QA.csv", 
                                     stringsAsFactors = FALSE, header = FALSE)

# read in wavelength metadata for this JAZ instrument

JAZ_wavelengths = read.csv("/Users/jrcollins/Code/Optics_Photochem/JAZ_wavelengths.csv",
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

# NOAA data

# data is in separate files, so have to read them in sequentially

setwd(initial.wd)
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

save(PAL1314_NOAA_AntUV_spectra_uW_cm2, file = "Incident_UV-VIS_spectra_PAL1314_uW_cm2.RData")

##### use data from 11-13 Nov 11 (when JAZ was deployed in open air atop tank) for intercalibration with NOAA radiometer #####

# load in JAZ data

setwd(initial.wd)
setwd("data/raw/JAZ_UV_VIS")

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2 = read.csv("JAZ_UV-VIS_open_air_full_spectra_PAL1314_11_13Nov13_uW_cm2.csv", 
                                                                   stringsAsFactors = FALSE, header = FALSE)

# assign column names (assumes JAZ_wavelenths has already been loaded above)

colnames(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2) = c("Date_raw_julian",
                                                                      "Inttime_microseconds",
                                                                      "Badscans_fullspectrum",
                                                                      "Badscans_UVB",
                                                                      JAZ_wavelengths$V1)

# create timestamp

# convert dates from messed up Excel format
PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT = as.POSIXct(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Date_raw_julian*60*60*24, tz = "GMT", origin = "0000-01-01")-1*60*60*24

# match elements of NOAA data set with JAZ scans, then take a look at correlation

# subset NOAA data to range of interest

PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset = 
  PAL1314_NOAA_AntUV_spectra_uW_cm2[
    PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT>=
      as.POSIXct('2013-11-11 12:00:00', tz = "GMT") & 
      PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT<=
      as.POSIXct('2013-11-13 12:00:00', tz = "GMT")
    ,]

# first, determine which elements to pull from each JAZ data string

JAZ.match.index.intercal = vector(length = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)-1), mode = "numeric")

for (i in 1:length(JAZ.match.index.intercal)) {
  
  JAZ.match.index.intercal[i] =
    which((abs(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[i+1])-
                 JAZ_wavelengths$V1)
           ==min(abs(JAZ_wavelengths$V1-as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[i+1])))))
  
}

# preallocate matrices to hold data on surface penetration

PAL1314_intercal= as.data.frame(matrix(data = NA, ncol = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)-1), nrow = nrow(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)))
colnames(PAL1314_intercal) = colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset)]

# append timestamps

PAL1314_intercal$Timestamp_GMT = PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[,1]

# create matrices for scatterplot
PAL1314_intercal.NOAA = PAL1314_intercal
PAL1314_intercal.JAZ = PAL1314_intercal

# perform matching & calculations

for (i in 1:nrow(PAL1314_intercal)) {
  
  # find matching time in the JAZ data, if it exists
  # we'll want something +/- 30 seconds of the timestamp in the NOAA data
  
  possible.time.match =
    which((abs(PAL1314_intercal$Timestamp_GMT[i]-
                 PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT)
           ==min(abs(PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT
                     -PAL1314_intercal$Timestamp_GMT[i]))))
  
  if (as.numeric(abs(PAL1314_intercal$Timestamp_GMT[i]-
          PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2$Timestamp_GMT[possible.time.match]))<=
      15) { # we have a "matching" reading within 30 seconds --> calculate ratio & record
    
    PAL1314_intercal[i,1:(ncol(PAL1314_intercal)-1)] = 
      
      PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2[possible.time.match,JAZ.match.index.intercal]/
      PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[i,2:ncol(PAL1314_intercal)]
    
    # also, dump data to our two matrices for scatterplot
    
    PAL1314_intercal.NOAA[i,1:(ncol(PAL1314_intercal)-1)] =
      PAL1314_NOAA_AntUV_spectra_uW_cm2.intercal.subset[i,2:ncol(PAL1314_intercal)]
    PAL1314_intercal.JAZ[i,1:(ncol(PAL1314_intercal)-1)] =
      PAL1314_JAZ_open_air_full_spectrum_11to13Nov13_uW_cm2[possible.time.match,JAZ.match.index.intercal]
  }
  
}

# calculate mean values for each wavelength, plot

PAL1314_NOAA.mean = apply(PAL1314_intercal.NOAA[,1:(ncol(PAL1314_intercal.NOAA)-1)],2,mean,na.rm = TRUE)
PAL1314_NOAA.sd = apply(PAL1314_intercal.NOAA[,1:(ncol(PAL1314_intercal.NOAA)-1)],2,sd,na.rm = TRUE)
PAL1314_JAZ.mean = apply(PAL1314_intercal.JAZ[,1:(ncol(PAL1314_intercal.JAZ)-1)],2,mean,na.rm = TRUE)
PAL1314_JAZ.sd = apply(PAL1314_intercal.JAZ[,1:(ncol(PAL1314_intercal.JAZ)-1)],2,sd,na.rm = TRUE)
plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
     PAL1314_NOAA.mean,
     xlim = c(290,600),type="l",col="blue")
# polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
#         c(PAL1314_NOAA.mean+PAL1314_NOAA.sd,rev(PAL1314_NOAA.mean-PAL1314_NOAA.sd)),
#         col="lightskyblue")
# polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
#         c(PAL1314_JAZ.mean+PAL1314_JAZ.sd,rev(PAL1314_JAZ.mean-PAL1314_JAZ.sd)),
#         col="lightsalmon")
lines(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
      PAL1314_NOAA.mean,
      col="blue")
lines(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
      PAL1314_JAZ.mean,col="red")

# plot intercalibration ratio

PAL1314_intercal.mean = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,mean,na.rm = TRUE)
plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
     PAL1314_intercal.mean,
     ylim = c(-3,3),
     xlim = c(290,600))

PAL1314_intercal.mean = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,mean,na.rm = TRUE)
PAL1314_intercal.sd = apply(PAL1314_intercal[,1:(ncol(PAL1314_intercal)-1)],2,sd,na.rm = TRUE)
plot(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
     PAL1314_intercal.mean, type = "p", pch = 16, cex = 0.3, 
     ylim = c(-10,10),
     xlim = c(290,600))
polygon(c(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),rev(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]))),
  c(PAL1314_intercal.mean+PAL1314_intercal.sd,rev(PAL1314_intercal.mean-PAL1314_intercal.sd)),
  col=grey(0.95))
points(as.numeric(colnames(PAL1314_intercal)[1:(ncol(PAL1314_intercal)-1)]),
     PAL1314_intercal.mean, pch = 16, cex = 0.3)

NOAA.caldat = unlist(PAL1314_intercal.NOAA[1:(ncol(PAL1314_intercal.NOAA)-1)])
NOAA.caldat = NOAA.caldat[!is.na(NOAA.caldat)]
NOAA.caldat = as.matrix(NOAA.caldat)

JAZ.caldat = unlist(PAL1314_intercal.JAZ[1:(ncol(PAL1314_intercal.JAZ)-1)])
JAZ.caldat = JAZ.caldat[!is.na(JAZ.caldat)]
JAZ.caldat = as.matrix(JAZ.caldat)

# write.csv(NOAA.caldat,file="NOAA.caldat.csv")
# write.csv(JAZ.caldat,file="JAZ.caldat.csv")

plot(NOAA.caldat/JAZ.caldat,ylim = c(-3,3))

fit.lm = lm(JAZ.caldat~NOAA.caldat)
fit.res = resid(fit.lm)
plot(NOAA.caldat,fit.res)

fit.lm = lm(JAZ.caldat~NOAA.caldat)
summary(fit.lm)

plot(NOAA.caldat,
     JAZ.caldat,
     pch=16,cex=0.1)
abline(0,1,col="red")
abline(fit.lm, col="blue")

##### match elements of NOAA data set with JAZ scans, then calculate surface penetration ##### 

# will use the NOAA data as our index since they were sampled at more sparse
# time & wavelength intervals

# first, determine which elements to pull from each JAZ data string

JAZ.match.index = vector(length = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)-1), mode = "numeric")
  
for (i in 1:length(JAZ.match.index)) {
  
  JAZ.match.index[i] =
    which((abs(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[i+1])-
              JAZ_wavelengths$V1)
        ==min(abs(JAZ_wavelengths$V1-as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[i+1])))))
  
}

# preallocate matrix to hold data on surface penetration

PAL1314_Surf_pen = as.data.frame(matrix(data = NA, ncol = (ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)-1), nrow = nrow(PAL1314_NOAA_AntUV_spectra_uW_cm2)))
colnames(PAL1314_Surf_pen) = colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2)]

# append timestamps

PAL1314_Surf_pen$Timestamp_GMT = PAL1314_NOAA_AntUV_spectra_uW_cm2[,1]

# perform matching & calculations

for (i in 1:nrow(PAL1314_Surf_pen)) {
  
  # find matching time in the JAZ data, if it exists
  # we'll want something +/- 30 seconds of the timestamp in the NOAA data
  
  possible.time.match =
    which((abs(PAL1314_Surf_pen$Timestamp_GMT[i]-
                 PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT)
           ==min(abs(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT
                     -PAL1314_Surf_pen$Timestamp_GMT[i]))))
  
  if (as.numeric(abs(PAL1314_Surf_pen$Timestamp_GMT[i]-
           PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT[possible.time.match]))<=
      30*(1/24/60/60)) { # we have a "matching" reading within 30 seconds --> calculate surface penetration (fraction) & record
    
    PAL1314_Surf_pen[i,1:(ncol(PAL1314_Surf_pen)-1)] = 
        
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[possible.time.match,JAZ.match.index]/
      PAL1314_NOAA_AntUV_spectra_uW_cm2[i,2:ncol(PAL1314_Surf_pen)]
      
  }
  
}

# calculate mean values for each wavelength

PAL1314_Surf_pen.mean = apply(PAL1314_Surf_pen[,1:(ncol(PAL1314_Surf_pen)-1)],2,mean,na.rm = TRUE)
plot(as.numeric(colnames(PAL1314_Surf_pen)[1:(ncol(PAL1314_Surf_pen)-1)]),
     PAL1314_Surf_pen.mean,
     ylim = c(0,3),
     xlim = c(290,600))

##### Kd calc & other workup using static depth UV-VIS profile from 2015-2016 season ##### 

# load in JAZ profile data (made on 15 Dec 2015 at Station B, Arthur Harbor)

setwd(initial.wd)
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

setwd(initial.wd)
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

# full spectrum 290-740 nm

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))

color_ramp.depths = colorRampPalette(c("cadetblue1", "darkblue")) # create a color ramp for shading the points according to lipid class

depth.plot.colors = rev(color_ramp.depths(nrow(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)))

plot(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[1,270:1500],"l",
     col = depth.plot.colors[1], lty = 1, lwd = "1.5",
     ylim = c(0,25), xlim = c(290,720),
     ylab = expression(paste("Irradiance (",mu,"W ",cm^-2," ",nm^-1,")")),
     xlab = "Wavelength (nm)")

text(653,
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[1,1300],
     labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[1]),
     offset = -.25, pos = 3)

# overlay other depths

for (i in c(3,5,7,9,11,13,15,16,17)) {
  
  lines(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
        PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,270:1500],
        col = depth.plot.colors[i], lty = 1, lwd = "1.5")
  
  if (i %in% c(7,9,11,13,15,16,17)) {
  text(653,
       PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,1300],
       labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[i]),
       offset = -.25, pos = 3)
    
  }
  
}

# overlay concurrent reading from NOAA radiometer at Palmer, if desired
# doesnt match up exactly 

PAL1314_NOAA_AntUV_spectra_uW_cm2.sub =
  PAL1314_NOAA_AntUV_spectra_uW_cm2[which(abs(PAL1314_NOAA_AntUV_spectra_uW_cm2$Timestamp_GMT-
                                                as.POSIXct('2013-12-15 11:11:00', tz = "GMT")) < 10),]

lines(as.numeric(colnames(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)]),
     PAL1314_NOAA_AntUV_spectra_uW_cm2.sub[2:ncol(PAL1314_NOAA_AntUV_spectra_uW_cm2.sub)],
     col="black")

dev.off()

# UVB-range inset

# make plot for inset (290-315 nm)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215_inset.pdf",
    width = 4, height = 4, pointsize = 12,
    bg = "white")

plot(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,270:1500],"l",
     col = depth.plot.colors[1], lty = 1, lwd = "1.5",
     ylim = c(0,5), xlim = c(290,315),
     ylab = expression(paste("Irradiance (",mu,"W ",cm^-2," ",nm^-1,")")),
     xlab = "Wavelength (nm)")

text(315,
     PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,338],
     labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[3]),
     offset = -.5, pos = 3)

# overlay other depths

for (i in c(7,11,15,16,17)) {
  
  lines(as.numeric(colnames(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected)[270:1500]),
        PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,270:1500],
        col = depth.plot.colors[i], lty = 1, lwd = "1.5")
  
  if (i %in% c(7,9,11,13,15,16,17)) {
    text(315,
         PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[i,338],
         labels = c(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected$Depth_m[i]),
         offset = -.5, pos = 3)
    
  }
  
}

dev.off()

# now, calculate Kds for various wavelengths

# preallocate matrix to hold data

PAL1516_AH_Kd_20151215 = as.data.frame(matrix(data = NA, ncol = nrow(JAZ_wavelengths), nrow = 1))
colnames(PAL1516_AH_Kd_20151215) = JAZ_wavelengths$V1

# use 0.5 m (assumed to be the subsurface incident irradiance, E(0))
# and 7 m depth irradiance data (as E(z)) to calculate Kd with z = 7

for (i in 1:ncol(PAL1516_AH_Kd_20151215)) {
  
PAL1516_AH_Kd_20151215[1,i] = 
  
  (log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[16,i+5]-
     log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[3,i+5]))/7)

}

# save corrected irradiance data and Kds 

save(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected, file = "JAZ_UV_VIS_corrected_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.RData")

write.csv(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected, file = "JAZ_UV_VIS_corrected_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.csv")

PAL1516_AH_Kd_20151215_per_meter = t(PAL1516_AH_Kd_20151215)

colnames(PAL1516_AH_Kd_20151215_per_meter) = c("Kd_per_meter")

save(PAL1516_AH_Kd_20151215_per_meter, file = "Kd_PAL1516_Arthur_Hbr_Stn_B_20151215.RData")

write.csv(PAL1516_AH_Kd_20151215_per_meter, file = "Kd_PAL1516_Arthur_Hbr_Stn_B_20151215.csv")

# plot

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215_Kd.pdf",
    width = 4, height = 3, pointsize = 12,
    bg = "white")

plot(JAZ_wavelengths$V1,
     PAL1516_AH_Kd_20151215,"l",
     col = "black", lty = 1, lwd = "1.5",
     ylim = c(0.05,0.5), xlim = c(290,720),
     ylab = expression(paste(K_d)),
     xlab = "Wavelength (nm)")

dev.off()

# just UV wavelengths

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "AH_depth_profile_20151215_Kd_UV.pdf",
    width = 4, height = 3, pointsize = 12,
    bg = "white")

plot(JAZ_wavelengths$V1,
     PAL1516_AH_Kd_20151215,"l",
     col = "black", lty = 1, lwd = "1.5",
     ylim = c(0.05,0.5), xlim = c(290,315),
     ylab = expression(K_d),
     xlab = "Wavelength (nm)")

dev.off()

### other panels for Experiment 13 plot ###

# instantaneous/cumulative UVB dosage

# subset to ranges of interest

# JAZ

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT>=
      as.POSIXct('2013-12-14 12:32:00', tz = "GMT") &  # 9:30 local time
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT<=
      as.POSIXct('2013-12-14 20:50:00', tz = "GMT") # 17:50 local time
      ,]

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec[,
        c(1:4,
          4+which(JAZ_wavelengths>=290 & JAZ_wavelengths<=315),
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

# instantaneous UVB doses

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB = 
  apply(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB)-1)],1,sum)

PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB = 
  apply(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.UVB,1,sum)

# at 0.6 m, from JAZ

plot(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT,
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB,col="red",
  ylim = c(0,500)) 

# superimpose points from NOAA data
points(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec$Timestamp_GMT,
       PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB,col="blue")

# out of curiosity, take a look at average ratio

mean.PAL1314.14Dec.JAZ.instUVB = 
  mean(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB)

mean.PAL1314.14Dec.NOAA.instUVB = 
  mean(PAL1314_NOAA_AntUV_spectra_uW_cm2.14Dec.instUVB)

# calculate cumulative dosage

# need library caTools, first have to detach RSEIS since there's a conflict with a function name

detach("package:RSEIS", unload=TRUE)
library(caTools)

# some definitions

W_per_uW = 1/1000000
J_per_kJ = 1/1000
cm2_per_m2 = 10000

# preallocate vector

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB =
  vector(length = length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB))

# vector of times in seconds, with t = 0 being timepoint at beginning of experiment

Exp13_timeint.s = as.numeric(rev(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT[nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB)]-
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.UVB$Timestamp_GMT))

for (i in 1:length(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB[i] =
    
    caTools::trapz(Exp13_timeint.s[1:i],PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.instUVB[1:i])

}

PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB.kJ_m2 =
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.cumUVB*W_per_uW*J_per_kJ*cm2_per_m2

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
