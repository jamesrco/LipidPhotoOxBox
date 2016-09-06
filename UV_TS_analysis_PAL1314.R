# UV_TS_analysis_PAL1314.R
#
# Purpose: Analysis & treatment of UV timeseries data from Palmer 2013-2014
# field season. Uses NOAA ESRL GMD Antarctic UV irradiance data and data from
# Jaz UV-VIS spectrophotometer J.R.C. deployed at Palmer Station during the 
# 2013-2014 field season
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

# assuming initial wd is that of the LipidPhotoOxBox repository, get current
# working directory so we can reset it later
initial.wd = getwd()

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
setwd("data/nice/JAZ_UV-VIS") 

DD_UVB_1314_subsurf = read.csv("Daily_int_UVB_dose_0.6m_subsurface_Palmer_kJ_m2.csv", 
         stringsAsFactors = FALSE)
colnames(DD_UVB_1314_subsurf)[1] = c("Date_raw")

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

##### wavelength-specific ("hi-res") data (for Kd calcs) ##### 

# JAZ data

setwd(initial.wd)
setwd("data/nice/JAZ_UV-VIS") 

# read in data, from .csv output generated in MATLAB (see the MATLAB script
# JAZ_data_read.m in the GitHub repo https://github.com/jamesrco/Optics_Photochem)

UVB_PAL1314_subsurf_hires_uW_cm2 = read.csv("UVB_spectra_0.6m_subsurface_PAL1314_uW_cm2.csv", 
                                     stringsAsFactors = FALSE, skip = 4)

# NOAA data

# need to track down hi-res Palmer spectra from 2013, not currently on NOAA web site
