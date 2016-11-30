# PAL1314_lipid_photoox_turnover.R

# Purpose: Combine in situ UV irradiance data, lipid data from PAL1314_PAL1516_environmental_samples.R, and AQYs calculated in PAL1314_AQY_calc.R to calculate potential lipid turnover rates on WAP during study period

# ****** Assumes all variables created by other scripts in Git repo LipidPhotoOxBox are already
# in user's workspace; these scripts can be found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 11/23/16 by J.R.C.

# first, pull out data (in pmol/L) for the PUFA and very-high-PUFA fractions of the PAL1314 & LMG1401 environmental samples, and the Marchetti diatom samples, analyzed in PAL1314_PAL1516_environmental_samples.R

PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L = PAL1314_LMG1401.sat_totals[c(4,5),]
PUFA_fracs.PAL1314_LMG1401_PConly.pmol_lipids_L = PAL1314_LMG1401.sat_totals.PC[c(4,5),]
PUFA_fracs.PAL1314_LMG1401.pmol_C_L = PAL1314_LMG1401.total_pmol_C.bysatclass[c(4,5),]

# PUFA_fracs.Marchetti_diatoms.pmol_total_lipids = Marchetti_diatoms.sat_totals[c(4,5),]
# PUFA_fracs.Marchetti_diatoms.pmol_total_C = Marchetti_diatoms.total_pmol_C.bysatclass[c(4,5),]
# PUFA_fracs.Marchetti_diatoms_PConly.pmol_total_lipids = Marchetti_diatoms.sat_totals.PC[c(4,5),]

# obtain list of dates on which we had a full day of Jaz data
# these are the dates for which we will calculate the lipid photooxidation rates

Jazfulldata.dates = DD_UVB_1314_subsurf$Date[!is.nan(DD_UVB_1314_subsurf$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2)]

# define some constants

z_surf_m = 5 # 5 m surface layer depth
m3_per_L = 0.001
cm_per_m = 100
mol_per_pmol = 1/(10^12)
pmol_per_mol = 10^12

# need package caTools, if not loaded already

detach("package:RSEIS", unload=TRUE)
library(caTools)

# calculations
# per Ch 4, Equations 9 and then 6

# preallocate matrices for results

# high-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)


# high-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

# mid-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)


# high-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 = as.data.frame(matrix(NA,length(Jazfulldata.dates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
rownames(C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401) = Jazfulldata.dates
colnames(C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)


for (i in 1:length(Jazfulldata.dates)) {
  
  # iterate through dates on which we had a full day of Jaz data
  
  #### compute integrated photon fluxes for this day ####
  
  # first, extract in situ irradiance data and compute daily integrated photon fluxes (by wavelength) for this date
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
  
  #### calculate Ks ####
  # per Ch 4, Equation 9
  
  # Ks will be in units of mol photons/mol rxn/time (time is implied; it is the 24-hour day)
  
  Ks_thisdate = matrix(NA, ncol = 1, nrow =
                         length(PAL1314.JAZ.sub_umol_photons_m2_d_nm))
  
  for (l in 1:nrow(Ks_thisdate)) { # iterate by wavelength
    
    if (lambda_nm_JAZ[l]<800) { # since don't have any lipid absorbance data above 800 nm
      
      # retrieve necessary values for this wavelength; some manipulations since
      # in some cases, variables were sampled at different intervals from each other
      
      # time-integrated photon flux at this wavelength (over sampling interval) 
      E_mol_photons_m2 = PAL1314.JAZ.sub_umol_photons_m2_d_nm[l]/1000000
      
      # decadic molar absorption coefficient from LipidAbsData$epsilon_M_cm_PC22_6_dil
      # calculated in LipidUV_VISAbsPlots.R
      epsilon_per_M_per_cm = LipidAbsData$epsilon_M_cm_PC22_6_dil[LipidAbsData$lambda_nm==round(lambda_nm_JAZ[l])]
      
      # downwelling attenuation coefficient (Kd) from Arthur Harbor depth profile
      Kd_per_m = PAL1516_AH_Kd_20151215_per_meter[abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-lambda_nm_JAZ[l])==min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-lambda_nm_JAZ[l]))]
      
      # finally, calculate Ks for this wavelength
      # have to make some unit conversions
      
      # this is Ch 4, Equation 9
      
      Ks_thisdate[l] =
        (E_mol_photons_m2*epsilon_per_M_per_cm*cm_per_m*m3_per_L*(1-10^(-Kd_per_m*z_surf_m)))/
        (Kd_per_m*z_surf_m)
      
    }
    
  }
  

  #### integrated Ks ####
  
  # calculate integrated Ks over UVB range and TUVR range
  # term under integral in denominator of Equation 6
  
  Ks_int_UVB = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315)],Ks_thisdate[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315),1])
  
  Ks_int_TUVR = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5)],Ks_thisdate[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5),1])

  #### photooxidation rate calcs ####
  
  # now, can calculate -d[i]/dt in Equation 6 (the photooxidation rates) for each sample in the PAL1314/LMG1401 dataset
  # also, scale up to pmoles of C
  
  # hi-PUFA fraction 
  
  # UVB
  
  for (m in 1:ncol(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401)) {
    
    PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[i,m] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,m]*mol_per_pmol)*
         Ks_int_UVB*Theta_Exp13_PC_22_6[3])*pmol_per_mol
    
    PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[i,m] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,m]*mol_per_pmol)*
         Ks_int_UVB*Theta_Exp13_PC_22_6_sigma[3])*pmol_per_mol
    
    C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[i,m] = 
      PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[i,m]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,m]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,m])
    
    C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[i,m] =
      PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[i,m]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,m]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,m])
    
  }
  
  # TUVR
  
  for (n in 1:ncol(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401)) {
    
    PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[i,n] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,n]*mol_per_pmol)*
         Ks_int_TUVR*Theta_Exp13_PC_22_6[1])*pmol_per_mol
    
    PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[i,n] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,n]*mol_per_pmol)*
         Ks_int_TUVR*Theta_Exp13_PC_22_6_sigma[1])*pmol_per_mol
    
    C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[i,n] = 
      PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[i,n]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,n]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,n])
    
    C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[i,n] =
      PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[i,n]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,n]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,n])
    
  }
  
  # mid-PUFA fraction 
  
  # UVB
  
  for (o in 1:ncol(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401)) {
    
    PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[i,o] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,o]*mol_per_pmol)*
         Ks_int_UVB*Theta_Exp13_PC_22_6[3])*pmol_per_mol
    
    PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[i,o] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,o]*mol_per_pmol)*
         Ks_int_UVB*Theta_Exp13_PC_22_6_sigma[3])*pmol_per_mol
    
    C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[i,o] = 
      PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[i,o]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,o]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,o])
    
    C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[i,o] =
      PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[i,o]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,o]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,o])
    
  }
  
  # TUVR
  
  for (p in 1:ncol(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401)) {
    
    PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[i,p] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,p]*mol_per_pmol)*
         Ks_int_TUVR*Theta_Exp13_PC_22_6[1])*pmol_per_mol
    
    PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[i,p] =
      ((PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,p]*mol_per_pmol)*
         Ks_int_TUVR*Theta_Exp13_PC_22_6_sigma[1])*pmol_per_mol
    
    C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[i,p] = 
      PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[i,p]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,p]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,p])
    
    C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[i,p] =
      PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[i,p]*
      (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,p]/
         PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,p])
    
  }
  
}

#### workup and visualization of results ####

# # high-PUFA fraction, all lipids, UVB
# PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 
# PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 
# C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 
# C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 
# 
# # high-PUFA fraction, all lipids, TUVR
# PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401
# PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 
# C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401
# C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 
# 
# # mid-PUFA fraction, all lipids, UVB
# PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401 
# PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401
# C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401
# C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401
# 
# # mid-PUFA fraction, all lipids, TUVR
# PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401
# PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401
# C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401
# C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401

# convert the numbers in the high-PUFA fraction to units of ug C/m3/d

# high-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 
PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 
C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3
C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3

# high-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401
PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 
C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3
C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3

# mid-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401
PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 
C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3
C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3

# make some objects with these data with continuous dates (no gaps)

Jazfulldata.dates.nogaps = DD_UVB_1314_subsurf$Date

C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

for (i in 1:nrow(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps)) {
  
  if (rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps)[i] %in% 
      rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401)) {
    
    C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401),]
    
  }
  
}

# pull in PAL1314 BP data, for comparison

PAL1314.BP = read.csv("/Users/jrcollins/Code/LipidPhotoOxBox/data/raw/PAL1314_LMG1401_PAL_LTER_data/PAL1314 Bacterial Production (Station).csv", skip = 1)

# convert leucine numbers to units of ug C/m3/d, using most conservative assumptions (ID = 1, conversion factor of 1.5 kg C/mol leu)

PAL1314.BP$BP_ugC_m3_d = PAL1314.BP$Leucine.Incorp...pmol.L.hr. * 36

# make a plot; use most proximate lipid concentration data (temporally)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Lipid_AQY_results_w_BP.pdf",
    width = 6.5, height = 3.5, pointsize = 10,
    bg = "white")

par(mar=c(5,5,1,5))

plot(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(0,1000),
     xlim = c(as.numeric(Jazfulldata.dates.nogaps[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "ug C per m3 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
     )

axis.POSIXct(1, at = seq(Jazfulldata.dates.nogaps[1], as.POSIXct("2014-01-06"), by = "5 days"), format = "%d %b")

polygon(c(Jazfulldata.dates[1:27],
          rev(Jazfulldata.dates[1:27])),
        c(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]+
            C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27],
          rev(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]-
                C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27])),
        border = NA, col = "mistyrose")

polygon(c(Jazfulldata.dates[27:35],
          rev(Jazfulldata.dates[27:35])),
        c(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]+
            C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35],
          rev(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]-
                C_xformed_ug_C_m3_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35])),
        border = NA, col = "mistyrose")

polygon(c(Jazfulldata.dates[1:27],
          rev(Jazfulldata.dates[1:27])),
        c(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]+
            C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27],
          rev(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]-
                C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27])),
        border = NA, col = "azure")

polygon(c(Jazfulldata.dates[27:35],
          rev(Jazfulldata.dates[27:35])),
        c(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]+
            C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35],
          rev(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]-
                C_xformed_ug_C_m3_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35])),
        border = NA, col = "azure")

lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred")

lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], type = "o", pch = 22, bg = "black", cex = 0.4, col = "darkred")
# lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40], lty=2, type = "o", pch = 22, bg = "white")
# lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], lty=2, type = "o", pch = 22, bg = "white")

lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40], lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

# superimpose BP data

# extract relevant data

BP.subset.B = PAL1314.BP[(PAL1314.BP$Depth..m.==0 & PAL1314.BP$Station.Name %in% c("B")),]
BP.subset.E = PAL1314.BP[(PAL1314.BP$Depth..m.==0 & PAL1314.BP$Station.Name %in% c("E")),]
BP.subset.SWI = PAL1314.BP[PAL1314.BP$Station.Name=="SWI",]

BP.subset.SWI.means = as.data.frame(matrix(NA,length(unique(BP.subset.SWI$Date.GMT)),2))
colnames(BP.subset.SWI.means) =  c("Date.GMT","BP_ugC_m3_d")
BP.subset.SWI.means$Date.GMT = unique(BP.subset.SWI$Date.GMT)

for (i in 1:nrow(BP.subset.SWI.means)) {
  
  BP.subset.SWI.means$BP_ugC_m3_d[i] =
    mean(BP.subset.SWI$BP_ugC_m3_d[BP.subset.SWI$Date.GMT==BP.subset.SWI.means$Date.GMT[i]])
  
}

points(as.POSIXct(as.character(BP.subset.E$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.E$BP_ugC_m3_d, pch = 22, bg = "black", cex = 1.2)

points(as.POSIXct(as.character(BP.subset.SWI.means$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.SWI.means$BP_ugC_m3_d, pch = 22, bg = "white", cex = 1.2)

points(as.POSIXct(as.character(BP.subset.B$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.B$BP_ugC_m3_d, pch = 21, bg = "black", cex = 1.2)

dev.off()

# second plot with the higher-value BP samples


par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Lipid_AQY_results_high_BP.pdf",
    width = 6.5, height = 3.5, pointsize = 10,
    bg = "white")

par(mar=c(5,5,1,5))

plot(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(1000,2500),
     xlim = c(as.numeric(Jazfulldata.dates.nogaps[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "ug C per m3 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
)

# superimpose BP data

# extract relevant data

BP.subset.B = PAL1314.BP[(PAL1314.BP$Depth..m.==0 & PAL1314.BP$Station.Name %in% c("B")),]
BP.subset.E = PAL1314.BP[(PAL1314.BP$Depth..m.==0 & PAL1314.BP$Station.Name %in% c("E")),]
BP.subset.SWI = PAL1314.BP[PAL1314.BP$Station.Name=="SWI",]

BP.subset.SWI.means = as.data.frame(matrix(NA,length(unique(BP.subset.SWI$Date.GMT)),2))
colnames(BP.subset.SWI.means) =  c("Date.GMT","BP_ugC_m3_d")
BP.subset.SWI.means$Date.GMT = unique(BP.subset.SWI$Date.GMT)

for (i in 1:nrow(BP.subset.SWI.means)) {
  
  BP.subset.SWI.means$BP_ugC_m3_d[i] =
    mean(BP.subset.SWI$BP_ugC_m3_d[BP.subset.SWI$Date.GMT==BP.subset.SWI.means$Date.GMT[i]])
  
}

points(as.POSIXct(as.character(BP.subset.E$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.E$BP_ugC_m3_d, pch = 22, bg = "black", cex = 1.2)

points(as.POSIXct(as.character(BP.subset.SWI.means$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.SWI.means$BP_ugC_m3_d, pch = 22, bg = "white", cex = 1.2)

points(as.POSIXct(as.character(BP.subset.B$Date.GMT), format = "%m/%d/%y", tz = "GMT"),
       BP.subset.B$BP_ugC_m3_d, pch = 21, bg = "black", cex = 1.2)

dev.off()

# some stats

# mean C:lipid (mol:mol) in the two IP-DAG fractions:

# mid-PUFA fraction
mean(unlist(PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,]/PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,]))
sd(unlist(PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,]/PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,]))

# high-PUFA fraction
mean(unlist(PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,]/PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,]))
sd(unlist(PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,]/PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,]))

# means of daily oxidation rate estimates

apply(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ug_C_m3_d.TUVR.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ug_C_m3_d.TUVR.hi_PUFA.PAL1314_LMG1401,2,sd)

# BP averages, untuil 1/2/14

mean(BP.subset.B$BP_ugC_m3_d[1:5])
sd(BP.subset.B$BP_ugC_m3_d[1:5])

mean(BP.subset.E$BP_ugC_m3_d[1:2])
sd(BP.subset.E$BP_ugC_m3_d[1:2])

mean(BP.subset.SWI.means$BP_ugC_m3_d[1:5])
sd(BP.subset.SWI.means$BP_ugC_m3_d[1:5])

