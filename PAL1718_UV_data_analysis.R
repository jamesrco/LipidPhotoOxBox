# PAL1718_UV_data_analysis.R

# for PAL1718 data

# We have actinometry this year. So, checking on some older PAL UV data to see
# what the daily fluxes look like in the active wavebands for the actinometers\\

# Ref for actinometry is Jankowski, Kieber, Mopper, Photochemistry and
# Photobiology, 1999, 70(3): 319–328

# Calculate some time-integrated doses from the NOAA data
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

# A new plot showing percent xmiss of the PC lipids and some materials together


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
