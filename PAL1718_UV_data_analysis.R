# PAL1718_UV_data_analysis.R

# for PAL1718 data

# We have actinometry this year. So, checking on some older PAL UV data to see
# what the daily fluxes look like in the active wavebands for the actinometers

# Ref for actinometry is Jankowski, Kieber, Mopper, Photochemistry and
# Photobiology, 1999, 70(3): 319–328

##### Calculate some time-integrated doses from the NOAA data #####

# *** requires the object PAL1314_est_spectra_at_depth_uW_cm2 from
# UV_TS_analysis_PAL1314.R ***

# preallocate an object to hold the data, will be equiv to E_n,p,sigma in Eq. 10
# in the now-rejected GCA manuscript

PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident =
  as.data.frame(matrix(NA,length(lipidox.calcdates),length(NOAA_AntUV_lambdas)))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident) = lipidox.calcdates
colnames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident) = NOAA_AntUV_lambdas

for (i in 1:length(lipidox.calcdates)) {
  
  # iterate through dates on which we want an estimate, from lipidox.calcdates
  
  # integrated figures in units of photons/cm2/wavelength
  
  # extract subset of data for this date
  
  PAL1314_est_spectra.today.uW_cm2 =  
    PAL1314_est_spectra_at_depth_uW_cm2[[1]][
      PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT>
        lipidox.calcdates[i]-60*60*12 & 
        PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT<=
        lipidox.calcdates[i]+60*60*12,]
  
  Todays.timestamps = 
    PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT[
      PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT>
        lipidox.calcdates[i]-60*60*12 & 
        PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT<=
        lipidox.calcdates[i]+60*60*12]
  
  # create a vector of times in seconds
  
  Timevector.s = Todays.timestamps-(lipidox.calcdates[i]-60*60*12)
  
  units(Timevector.s) = "secs"
  
  Timevector.s = as.numeric(Timevector.s)
  
  # convert to units of photons/cm2, step by step
  PAL1314_est_spectra.today.W_m2 = 
    PAL1314_est_spectra.today.uW_cm2/100
  
  PAL1314_est_spectra.today.umol_photons_m2_s =
    matrix(NA, ncol = length(NOAA_AntUV_lambdas), nrow = nrow(PAL1314_est_spectra.today.uW_cm2))
  
  for (k in 1:ncol(PAL1314_est_spectra.today.umol_photons_m2_s)) {
    
    PAL1314_est_spectra.today.umol_photons_m2_s[,k] =
      PAL1314_est_spectra.today.W_m2[,k]*NOAA_AntUV_lambdas[k]*0.836*(10^-2)
    
  }
  
  # integrate over time, by wavelength
  
  # first, time integration at each wavelength
  
  # preallocate vector
  PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today = vector(mode = "double",
                                                          length = ncol(PAL1314_est_spectra.today.umol_photons_m2_s))
  
  for (l in 1:length(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today)) {
    
    # requires vector of times in seconds, with t = 0 being beginning of day
    
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today[l] =
      caTools::trapz(Timevector.s[1:length(Timevector.s)],PAL1314_est_spectra.today.umol_photons_m2_s[,l])
    
  }
  
  # store data in appropriate location
  
  PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,] =
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today
  
}

# now, can calculate some integrals for certain wavebands

# nitrite actinometer response bw is 330-380 nm
# nitrate is 311-333

PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm =
  as.data.frame(matrix(NA,nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident),1))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm) =
  rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident)

for (i in 1:nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm)) {
  PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm[i,1] =
    caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=330 & NOAA_AntUV_lambdas<=380)],
                   as.numeric(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,NOAA_AntUV_lambdas>=330 & 
                                                                                  NOAA_AntUV_lambdas<=380]))
}

PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm =
  as.data.frame(matrix(NA,nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident),1))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm) =
  rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident)

for (i in 1:nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm)) {
  PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm[i,1] =
    caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=311 & NOAA_AntUV_lambdas<=333)],
                   as.numeric(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,NOAA_AntUV_lambdas>=311 & 
                                                                                  NOAA_AntUV_lambdas<=333]))
}

# for Jaz data from Dec 14 2013

PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2_330_380m =
  caTools::trapz(JAZ_wavelengths[(JAZ_wavelengths>=330 & JAZ_wavelengths<=380)],
                 as.numeric(PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2[JAZ_wavelengths>=330 & 
                                                                            JAZ_wavelengths<=380]))

PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2_311_333m =
  caTools::trapz(JAZ_wavelengths[(JAZ_wavelengths>=311 & JAZ_wavelengths<=333)],
                 as.numeric(PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2[JAZ_wavelengths>=311 & 
                                                                            JAZ_wavelengths<=333]))

##### A new plot showing percent xmiss of the PC lipids and some materials together #####


# Napierian molar absorption coefficients

absPlotCol = hsv(c(0.1, 0.35, 0.6, 0.85), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash","dotted") # define lty

transPlotCol = c("darkblue","skyblue","darkred","coral1","black") # define colors
transPlotLty = c("solid","dashed","dotted","dotdash","twodash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff_with_DHA_Napierian_plus_materials.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,2,5))
plot(LipidAbsData_init[1:801,c("lambda_nm")],log(LipidAbsData_init[1:801,c("kappa_M_cm_PC22_0")]),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,8), xlim = c(225,500),
     ylab = expression(paste("log ",kappa[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init[1:801,c("lambda_nm")],log(LipidAbsData_init[1:801,c("kappa_M_cm_PC22_1")]),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData_init[1:801,c("lambda_nm")],log(LipidAbsData_Jan17[1:801,c("kappa_M_cm_PC22_6_011387mM")]),
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")
lines(LipidAbsData_init[1:801,c("lambda_nm")],log(LipidAbsData_Jan17[1:801,c("kappa_M_cm_DHA_01032mM")]),
      col = absPlotCol[4], lty = absPlotLty[4], lwd = "1.5")

# append xmiss spectra 

par(new = T)
plot(PctTransData$lambda_nm,PctTransData$transmittance_quartz_pct,"l",
     col = transPlotCol[1], lty = transPlotLty[1], lwd = "1",
     ylim = c(0,100), xlim = c(225,500),
     axes=F, xlab=NA, ylab=NA)
lines(PctTransData$lambda_nm,PctTransData$transmittance_borosilicate_pct,
      col = transPlotCol[2], lty = transPlotLty[2], lwd = "2")
lines(PctTransData$lambda_nm,PctTransData$transmittance_tedlar_PVF_pct,
      col = transPlotCol[3], lty = transPlotLty[3], lwd = "2")
lines(PctTransData$lambda_nm,PctTransData$transmittance_PVF_PET_pct,
      col = transPlotCol[4], lty = transPlotLty[4], lwd = "2")
lines(as.numeric(as.character(PC_bottle[-c(1:8),1])),
      as.numeric(as.character(PC_bottle[-c(1:8),2])),
      col = transPlotCol[5], lty = transPlotLty[5], lwd = "2")

axis(side = 4)
mtext(side = 4, line = 3, "Percent transmittance")

legend(x = 370, y = 35, bty = "n",
       legend = c("Quartz glass vial","Borosilicate glass vial","PVF incubation bag",
                  "PVF bag w/PET screen","PC bottle"),
       col = transPlotCol, lty = transPlotLty, lwd = 2)

legend(x = 370, y = 60, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC","DHA"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

##### Taking a look at data from 2017 C-OPS profiles #####

library(Hmisc)

# read in data, extract profile

C_OPS_StnE_20171116.raw = mdb.get("data/raw/C_OPS/2017/Profile_171116_1216_stnE.mdb")
C_OPS_StnE_20171116.profile.raw = C_OPS_StnE_20171116.raw$Profile_001
C_OPS_StnB_20171116.raw = mdb.get("data/raw/C_OPS/2017/Profile_171116_1500_stnB.mdb")
C_OPS_StnB_20171116.profile.raw = C_OPS_StnB_20171116.raw$Profile_001
COPS_wavelengths_nm = read.csv("data/raw/C_OPS/2017/COPS_wavelengths_nm.csv",
                               skip = 1, header = FALSE, stringsAsFactors = FALSE)

# quick look at data

plot(C_OPS_StnB_20171116.profile.raw$EdZPAR,-C_OPS_StnB_20171116.profile.raw$LuZDepth,'l')
plot(C_OPS_StnE_20171116.profile.raw$EdZPAR,-C_OPS_StnE_20171116.profile.raw$LuZDepth,'l')

# some QA, per 
# Morrow, J.H., C.R. Booth, R.N. Lind, and S.B. Hooker, 2010: "The Compact
# Optical Profiling System (C-OPS)." In: J.H. Morrow, S.B. Hooker, C.R. Booth,
# G. Bernhard, R.N. Lind, and J.W. Brown, Advances in Measuring the Apparent
# Optical Properties (AOPs) of Optically Complex Waters, NASA Tech. Memo. 
# 2010–215856, NASA Goddard Space Flight Center, Greenbelt, Maryland, pgs 42–50.
# Available at: http://www.biospherical.com/images/pdf/cops-215856-chap4w.cover.pdf
# and
# Ian G. Brosnan, Calibration and data processing Biospherical C-Ops system used during the
# Ocean Optics 2015 Class Cruise. NASA Ames Research Center, Moffett Field, CA
# Document Version 1.0, August 4th, 2015 adapted from MURI Hawaii documentation
# by scott.freeman@nasa.gov; May 15, 2014)
# Available at: https://seabass.gsfc.nasa.gov/archive/MAINE/boss/OO2015/OO2015_Sampling/documents/Ocean_Optics_2015_COPS_Report.pdf

# subset data to exclude excessive pitch and roll

C_OPS_StnB_20171116.profile = subset(C_OPS_StnB_20171116.profile.raw,
                                     abs(EdZRoll) < 5 & abs(EdZPitch) < 5 & abs(Ed0Roll) < 5 & abs(Ed0Pitch) < 5)
C_OPS_StnE_20171116.profile = subset(C_OPS_StnE_20171116.profile.raw,
                                     abs(EdZRoll) < 5 & abs(EdZPitch) < 5 & abs(Ed0Roll) < 5 & abs(Ed0Pitch) < 5)

# correct for instrument offset (Ed - 0.05, Lu + 0.25 in most cases)

C_OPS_StnE_20171116.profile$EdZDepth.corr =
  C_OPS_StnE_20171116.profile$LuZDepth - 0.05
C_OPS_StnE_20171116.profile$LuZDepth.corr =
  C_OPS_StnE_20171116.profile$LuZDepth + 0.25
C_OPS_StnB_20171116.profile$EdZDepth.corr =
  C_OPS_StnB_20171116.profile$LuZDepth - 0.05
C_OPS_StnB_20171116.profile$LuZDepth.corr =
  C_OPS_StnB_20171116.profile$LuZDepth + 0.25

# subset data to exclude positive depth values

C_OPS_StnB_20171116.profile = subset(C_OPS_StnB_20171116.profile, EdZDepth.corr > 0)
C_OPS_StnE_20171116.profile = subset(C_OPS_StnE_20171116.profile, EdZDepth.corr > 0)

# subset data to exclude negative instrument readings

if (length(C_OPS_StnB_20171116.profile$EdZDepth.corr)>0) {
  C_OPS_StnB_20171116.profile[,10:65][ C_OPS_StnB_20171116.profile[,10:65] <= 0 ] = NA
}
if (length(C_OPS_StnE_20171116.profile$EdZDepth.corr)>0) {
  C_OPS_StnE_20171116.profile[,10:65][ C_OPS_StnE_20171116.profile[,10:65] <= 0 ] = NA
}

# subset data to exclude upcasting

if (length(C_OPS_StnB_20171116.profile$EdZDepth.corr)>0) {
  C_OPS_StnB_20171116.profile = C_OPS_StnB_20171116.profile[1:which.max(C_OPS_StnB_20171116.profile$EdZDepth.corr),]
  }
if (length(C_OPS_StnE_20171116.profile$EdZDepth.corr)>0) {
  C_OPS_StnE_20171116.profile = C_OPS_StnE_20171116.profile[1:which.max(C_OPS_StnE_20171116.profile$EdZDepth.corr),]
}

# take another look 

plot(C_OPS_StnB_20171116.profile$EdZPAR,-C_OPS_StnB_20171116.profile$EdZDepth.corr,'l',
     xlim = c(0,20))
lines(C_OPS_StnB_20171116.profile$EdZ443,-C_OPS_StnB_20171116.profile$EdZDepth.corr)
lines(C_OPS_StnB_20171116.profile$EdZ412,-C_OPS_StnB_20171116.profile$EdZDepth.corr)
lines(C_OPS_StnB_20171116.profile$EdZ340,-C_OPS_StnB_20171116.profile$EdZDepth.corr)

plot(C_OPS_StnE_20171116.profile$EdZPAR,-C_OPS_StnE_20171116.profile$EdZDepth.corr,'l',
     xlim = c(0,20))
lines(C_OPS_StnE_20171116.profile$EdZ443,-C_OPS_StnE_20171116.profile$EdZDepth.corr)
lines(C_OPS_StnE_20171116.profile$EdZ412,-C_OPS_StnE_20171116.profile$EdZDepth.corr)
lines(C_OPS_StnE_20171116.profile$EdZ340,-C_OPS_StnE_20171116.profile$EdZDepth.corr)

# check if our "0 m" depth cutoffs were sufficient (want to exclude bad data when 
# the instrument was being moved around at the surface)

plot(C_OPS_StnB_20171116.profile$EdZDepth.corr,C_OPS_StnB_20171116.profile$Ed0710)
plot(C_OPS_StnE_20171116.profile$EdZDepth.corr,C_OPS_StnE_20171116.profile$Ed0710)

# calculate some surface reference correction factors by which downwelling data can be corrected

# preallocate
PAL17_StnB_Kd_20171116.corrfactors = as.data.frame(matrix(data = NA, nrow = nrow(C_OPS_StnB_20171116.profile),
                                              ncol = 19)) 
PAL17_StnE_Kd_20171116.corrfactors = as.data.frame(matrix(data = NA, nrow = nrow(C_OPS_StnE_20171116.profile),
                                                          ncol = 19)) 

for (i in 1:nrow(PAL17_StnB_Kd_20171116.corrfactors)) {
  for (j in 1:19) {
    PAL17_StnB_Kd_20171116.corrfactors[i,j] =
      (C_OPS_StnB_20171116.profile[i,j+28])/mean(C_OPS_StnB_20171116.profile[,j+28])
  }
}

for (i in 1:nrow(PAL17_StnE_Kd_20171116.corrfactors)) {
  for (j in 1:19) {
    PAL17_StnE_Kd_20171116.corrfactors[i,j] =
      (C_OPS_StnE_20171116.profile[i,j+28])/mean(C_OPS_StnE_20171116.profile[,j+28])
  }
}

# now, calculate some Kds

# preallocate

PAL17_StnB_Kd_20171116.raw = as.data.frame(matrix(data = NA, nrow = nrow(C_OPS_StnB_20171116.profile),
                                              ncol = 20))
PAL17_StnE_Kd_20171116.raw = as.data.frame(matrix(data = NA, nrow = nrow(C_OPS_StnE_20171116.profile),
                                                  ncol = 20))

# set col names; collect metadata

colnames(PAL17_StnB_Kd_20171116.raw) =
  c("Depth_m",COPS_wavelengths_nm$V1)
PAL17_StnB_Kd_20171116.raw$Depth_m = C_OPS_StnB_20171116.profile$EdZDepth.corr
colnames(PAL17_StnE_Kd_20171116.raw) =
  c("Depth_m",COPS_wavelengths_nm$V1)
PAL17_StnE_Kd_20171116.raw$Depth_m = C_OPS_StnE_20171116.profile$EdZDepth.corr

# calculate Kds for each depth

for (j in 1:19) { # loop by wavelength
  
  # calculate the corrected below-surface irradiance for this wavelwngth
  
  I_0minus = mean(C_OPS_StnB_20171116.profile[C_OPS_StnB_20171116.profile$EdZDepth.corr<0.1,j+53]/
                    PAL17_StnB_Kd_20171116.corrfactors[C_OPS_StnB_20171116.profile$EdZDepth.corr<0.1,j])
  
  for (i in 1:nrow(PAL17_StnB_Kd_20171116.raw)) {
    
    I_Z = C_OPS_StnB_20171116.profile[i,j+53]/
      PAL17_StnB_Kd_20171116.corrfactors[i,j]
    
    PAL17_StnB_Kd_20171116.raw[i,j+1] =
      (log(I_Z)-log(I_0minus))/(0-C_OPS_StnB_20171116.profile$EdZDepth.corr[i])
    
  }
  
}

for (j in 1:19) { # loop by wavelength
  
  # calculate the corrected below-surface irradiance for this wavelwngth
  
  I_0minus = mean(C_OPS_StnE_20171116.profile[C_OPS_StnE_20171116.profile$EdZDepth.corr<0.1,j+53]/
                    PAL17_StnE_Kd_20171116.corrfactors[C_OPS_StnE_20171116.profile$EdZDepth.corr<0.1,j])
  
  for (i in 1:nrow(PAL17_StnE_Kd_20171116.raw)) {
    
    I_Z = C_OPS_StnE_20171116.profile[i,j+53]/
      PAL17_StnE_Kd_20171116.corrfactors[i,j]
    
    PAL17_StnE_Kd_20171116.raw[i,j+1] =
      (log(I_Z)-log(I_0minus))/(0-C_OPS_StnE_20171116.profile$EdZDepth.corr[i])
    
  }
  
}

# mean Kds for each wavelength

PAL17_StnB_Kd_20171116 = as.data.frame(matrix(data = NA, nrow = 19,
                                              ncol = 3))

colnames(PAL17_StnB_Kd_20171116) = c("Wavelength_nm","Kd_per_m.mean","Kd_per_m.sd")
PAL17_StnB_Kd_20171116$Wavelength_nm = COPS_wavelengths_nm$V1

for (i in 1:nrow(PAL17_StnB_Kd_20171116)) {
  
  PAL17_StnB_Kd_20171116$Kd_per_m.mean[i] = mean(PAL17_StnB_Kd_20171116.raw[150:nrow(PAL17_StnB_Kd_20171116.raw),i+1],
                                                 na.rm = T)
  PAL17_StnB_Kd_20171116$Kd_per_m.sd[i] = sd(PAL17_StnB_Kd_20171116.raw[150:nrow(PAL17_StnB_Kd_20171116.raw),i+1],
                                             na.rm = T)
  
}

PAL17_StnE_Kd_20171116 = as.data.frame(matrix(data = NA, nrow = 19,
                                              ncol = 3))

colnames(PAL17_StnE_Kd_20171116) = c("Wavelength_nm","Kd_per_m.mean","Kd_per_m.sd")
PAL17_StnE_Kd_20171116$Wavelength_nm = COPS_wavelengths_nm$V1

for (i in 1:nrow(PAL17_StnE_Kd_20171116)) {
  
  PAL17_StnE_Kd_20171116$Kd_per_m.mean[i] = mean(PAL17_StnE_Kd_20171116.raw[250:nrow(PAL17_StnE_Kd_20171116.raw),i+1],
                                                 na.rm = T)
  PAL17_StnE_Kd_20171116$Kd_per_m.sd[i] = sd(PAL17_StnE_Kd_20171116.raw[250:nrow(PAL17_StnE_Kd_20171116.raw),i+1],
                                             na.rm = T)
  
}


# compare these to the Jaz-derived values from PAL1516

Kd.matches_from_PAL1516 = as.data.frame(matrix(data = NA, nrow = 18, ncol = 3))
colnames(Kd.matches_from_PAL1516) = c("Wavelength_nm","Kd_per_meter","Kd_per_meter_fitted")

for (i in 1:nrow(Kd.matches_from_PAL1516)) {
  
  Kd.matches_from_PAL1516$Wavelength_nm[i] =
    as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))[
      which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-
                      as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[i])))]
  
  Kd.matches_from_PAL1516$Kd_per_meter[i] =
    PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter[
      which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-
                      as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[i])))]
  
  Kd.matches_from_PAL1516$Kd_per_meter_fitted[i] =
    PAL1516_AH_Kd_20151215_per_meter$Kd_per_meter_fitted[
      which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-
                      as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[i])))]
  
}

# take a look at a plot

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

# superimpose Kds from the 2017 C-OPS data
points(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18],PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18])

arrows(x0=as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18]),
       y0=PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18]-PAL17_StnB_Kd_20171116$Kd_per_m.sd[1:18],
       x1=as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18]),
       y1=PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18]+PAL17_StnB_Kd_20171116$Kd_per_m.sd[1:18],
       angle=90,
       code=3,
       length=0.04,
       lwd=0.4)

# points(PAL17_StnE_Kd_20171116$Wavelength_nm[1:18],PAL17_StnE_Kd_20171116$Kd_per_m.mean[1:18], pch = "+")

##### messing around #####

PAL1516_AH_Kd_20151215.new = as.data.frame(matrix(data = NA, ncol = nrow(JAZ_wavelengths), nrow = 14))
colnames(PAL1516_AH_Kd_20151215.new) = JAZ_wavelengths$V1

for (i in 346:1437) {
  
  for (j in 1:14) {
    
    PAL1516_AH_Kd_20151215.new[j,i] = 
      
      (log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[j,i+5])-
         log(PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2.corrected[16,i+5]))/-PAL1516_JAZ_Stn_B_profile_20151215_full_spectrum_uW_cm2$Depth_m[j]
    
  }
  
}


PAL1516_AH_Kd_20151215.new.mean = as.data.frame(matrix(data = NA, nrow = nrow(JAZ_wavelengths),
                                              ncol = 3))

colnames(PAL1516_AH_Kd_20151215.new.mean) = c("Wavelength_nm","Kd_per_m.mean","Kd_per_m.sd")
PAL1516_AH_Kd_20151215.new.mean$Wavelength_nm = JAZ_wavelengths$V1

for (i in 1:nrow(PAL1516_AH_Kd_20151215.new.mean)) {
  
  PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[i] = mean(PAL1516_AH_Kd_20151215.new[,i],
                                                 na.rm = T)
  PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.sd[i] = sd(PAL1516_AH_Kd_20151215.new[,i],
                                                        na.rm = T)
  
}


# set any of these calculated Kds with value > 1 = NA
PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean>=1] = NA

# now, fit an exponential model to estimate Kds for wavelengths 290-320 nm
# will use data from 320-370 nm as basis for fitting curve, then back-extrapolate

Kd_fit_subset.alt = log(PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[346:482])
Wavelength_fit_subset.alt = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[346:482])

fit.exp.alt = lm(Kd_fit_subset.alt ~ Wavelength_fit_subset.alt)

# get some information
summary(fit.exp.alt)

# predict values for 290-320 nm using model we just created

UVB_pred_subset.alt = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[264:346])
Kd_UVB.pred.alt = exp(predict(fit.exp.alt,list(Wavelength_fit_subset.alt=UVB_pred_subset.alt)))

# create vector of fitted values for range of data we used to fit the curve routine

Kd_fitrange_wavelengths.alt = as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter)[346:482])
Kd_fitrange.fitted.alt = exp(predict(fit.exp.alt,list(Wavelength_fit_subset.alt=Kd_fitrange_wavelengths.alt)))


# take a look at a plot

plot(JAZ_wavelengths$V1[483:nrow(JAZ_wavelengths)],
     PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[483:length(PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean)],
     "l",
     col = "black", lty = 1, lwd = "1",
     ylim = c(0.05,0.9), xlim = c(290,600),
     ylab = expression(paste(K_d)),
     xlab = "Wavelength (nm)")

# grey out but retain the range used for curve fitting
lines(JAZ_wavelengths$V1[346:482],
      PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[346:482],
      col = "lightgrey")

# fitted values
lines(UVB_pred_subset.alt,Kd_UVB.pred.alt,lty=3)
lines(Kd_fitrange_wavelengths.alt,Kd_fitrange.fitted.alt,lty=3)

# superimpose Kds from the 2017 C-OPS data
points(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18],PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18])

arrows(x0=as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18]),
       y0=PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18]-PAL17_StnB_Kd_20171116$Kd_per_m.sd[1:18],
       x1=as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[1:18]),
       y1=PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18]+PAL17_StnB_Kd_20171116$Kd_per_m.sd[1:18],
       angle=90,
       code=3,
       length=0.04,
       lwd=0.4)

# compare these to the Jaz-derived values from PAL1516

Kd.matches_from_PAL1516.alt = as.data.frame(matrix(data = NA, nrow = 18, ncol = 3))
colnames(Kd.matches_from_PAL1516.alt) = c("Wavelength_nm","Kd_per_meter","Kd_per_meter_fitted")

for (i in 1:nrow(Kd.matches_from_PAL1516.alt)) {
  
  Kd.matches_from_PAL1516.alt$Wavelength_nm[i] =
    as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))[
      which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-
                      as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[i])))]
  
  Kd.matches_from_PAL1516.alt$Kd_per_meter[i] =
    PAL1516_AH_Kd_20151215.new.mean$Kd_per_m.mean[
      which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-
                      as.numeric(PAL17_StnB_Kd_20171116$Wavelength_nm[i])))]
  
}

plot(Kd.matches_from_PAL1516.alt$Kd_per_meter,
     PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18])

plot(Kd.matches_from_PAL1516$Kd_per_meter,
     PAL17_StnB_Kd_20171116$Kd_per_m.mean[1:18])
