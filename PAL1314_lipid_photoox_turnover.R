# PAL1314_lipid_photoox_turnover.R

# Purpose: Combine in situ UV irradiance data, lipid data from PAL1314_PAL1516_environmental_samples.R, and AQYs calculated in PAL1314_AQY_calc.R to calculate potential lipid turnover rates on WAP during study period

# ****** Assumes all variables created by other scripts in Git repo LipidPhotoOxBox are already
# in user's workspace; these scripts can be found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 11/23/16 by J.R.C.

# first, pull out data (in pmol/L) for the PUFA and very-high-PUFA fractions of the PAL1314 & LMG1401 environmental samples, and the Marchetti diatom samples, analyzed in PAL1314_PAL1516_environmental_samples.R

PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L = PAL1314_LMG1401.sat_totals[c(4,5),]
PUFA_fracs.PAL1314_LMG1401_PConly.pmol_lipids_L = PAL1314_LMG1401.sat_totals.PC[c(4,5),]
PUFA_fracs.PAL1314_LMG1401.pmol_C_L = PAL1314_LMG1401.total_pmol_C.bysatclass[c(4,5),]

PUFA_fracs.Marchetti_diatoms.pmol_total_lipids = Marchetti_diatoms.sat_totals[c(4,5),]
PUFA_fracs.Marchetti_diatoms.pmol_total_C = Marchetti_diatoms.total_pmol_C.bysatclass[c(4,5),]
PUFA_fracs.Marchetti_diatoms_PConly.pmol_total_lipids = Marchetti_diatoms.sat_totals.PC[c(4,5),]

# define some constants

z_surf = 5 # 5 m surface layer depth

# now, calculate daily UVA+UVB and UVB-range photon fluxes from the Palmer deck tank Jaz data (0.6 m depth)

# will obtain integrated final figures in units of umol photons/m2/d

# need package caTools, if not loaded already

detach("package:RSEIS", unload=TRUE)
library(caTools)

# list of dates on which we had a full day of Jaz data
Jazfulldata.dates = DD_UVB_1314_subsurf$Date[!is.nan(DD_UVB_1314_subsurf$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2)]

# preallocate matrix for totals

PAL1314_daily_doses_0.6m_subsurf.umol_photons_m2_d = as.data.frame(matrix(NA,length(Jazfulldata.dates),2))
colnames(PAL1314_daily_doses_0.6m_subsurf.umol_photons_m2_d) = c("UVB_umol_photons_m2_d","UVA_UVB_umol_photons_m2_d")

# calculate the totals for each date

for (i in 1:nrow(PAL1314_daily_doses_0.6m_subsurf.umol_photons_m2_d)) {
  
  # extract Jaz wavelength-specific data for this date
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub = 
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2[
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT>=
        as.POSIXct(paste0(Jazfulldata.dates[i]," 03:00:00"), tz = "GMT") &  # 00:00 local time
        PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2$Timestamp_GMT<=
        as.POSIXct(paste0(Jazfulldata.dates[i]," 02:59:00"), tz = "GMT") + 86400 # 23:59 local time (02:59 GMT the next day)
      ,]
  
  # convert to units of photons/cm2, step by step
  PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.sub = 
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub)-1)]/100
  
  lambda_nm_JAZ = as.numeric(colnames(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub)[5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub)-1)])
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.sub =
    matrix(NA, ncol = length(lambda_nm_JAZ), nrow = nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub))
  
  for (j in 1:ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.sub)) {
    
    PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.sub[,j] =
      PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.sub[,j]*lambda_nm_JAZ[j]*0.836*(10^-2)
    
  }
  
  # integrate over time, by wavelength
  
  # first, time integration at each wavelength
  
  # need a vector of cumulative time elapsed at each timepoint in the subset, in seconds
  
  Daily_cumtimeint.s = as.numeric(rev(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub$Timestamp_GMT[nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub)]-
                                        PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.sub$Timestamp_GMT))
  
  # preallocate vector
  PAL1314.JAZ.sub_umol_photons_m2_d_nm = vector(mode = "double",
                                                length = ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.sub))
  
  for (k in 1:length(PAL1314.JAZ.sub_umol_photons_m2_d_nm)) {
    
    # requires vector of times in seconds, with t = 0 being beginning of day
    
    PAL1314.JAZ.sub_umol_photons_m2_d_nm[k] =
      caTools::trapz(Daily_cumtimeint.s[1:length(Daily_cumtimeint.s)],PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.sub[,k])
    
  }
  
  # record integrated figures for UVB and UVB+UVA range
  
  PAL1314_daily_doses_0.6m_subsurf.umol_photons_m2_d[i,1] =
    sum(PAL1314.JAZ.sub_umol_photons_m2_d_nm[lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315])
  
  PAL1314_daily_doses_0.6m_subsurf.umol_photons_m2_d[i,2] =
    sum(PAL1314.JAZ.sub_umol_photons_m2_d_nm[lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5])
  
}

# now, enter into AQY equation with these data, the PUFA-lipid concentrations, and the AQYs to get some estimates of total # of moles of PUFA-lipids transformed by UVR, and total # moles C (by way of the total # of moles of C bound up in each of lipid species)

# first, for the PAL1314_LMG1401 particulate data

# high-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_d.UVB.hi_PUFA.PAL1314_LMG1401 = matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L))
PUFA_lipids_xformed_pmol_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 = matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L))
C_xformed_pmol_d.UVB.hi_PUFA.PAL1314_LMG1401 = matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L))
C_xformed_pmol_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 = matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L))

for (i in 1:nrow(PUFA_lipids_xformed_pmol_d.UVB.hi_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_d.UVB.hi_PUFA.PAL1314_LMG1401[i] =
    Theta_Exp13_PC_22_6[3]
  
}