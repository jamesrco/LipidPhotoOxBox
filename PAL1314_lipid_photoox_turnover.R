# PAL1314_lipid_photoox_turnover.R

# Purpose: Combine in situ UV irradiance data, lipid data from PAL1314_PAL1516_environmental_samples.R, and AQYs calculated in PAL1314_AQY_calc.R to calculate potential lipid turnover rates on WAP during study period

# ****** Assumes all variables created by other scripts in Git repo LipidPhotoOxBox are already
# in user's workspace; these scripts can be found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 11/23/16 by J.R.C.

#### pull out some necessary lipid data; define constants ####

# first, pull out data (in pmol/L) for the PUFA and very-high-PUFA fractions of the PAL1314 & LMG1401 environmental samples, and the Marchetti diatom samples, analyzed in PAL1314_PAL1516_environmental_samples.R

PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L = PAL1314_LMG1401.sat_totals[c(4,5),]
PUFA_fracs.PAL1314_LMG1401_PConly.pmol_lipids_L = PAL1314_LMG1401.sat_totals.PC[c(4,5),]
PUFA_fracs.PAL1314_LMG1401.pmol_C_L = PAL1314_LMG1401.total_pmol_C.bysatclass[c(4,5),]

# PUFA_fracs.Marchetti_diatoms.pmol_total_lipids = Marchetti_diatoms.sat_totals[c(4,5),]
# PUFA_fracs.Marchetti_diatoms.pmol_total_C = Marchetti_diatoms.total_pmol_C.bysatclass[c(4,5),]
# PUFA_fracs.Marchetti_diatoms_PConly.pmol_total_lipids = Marchetti_diatoms.sat_totals.PC[c(4,5),]

# # obtain list of dates on which we had a full day of Jaz data
# # these are the dates for which we will calculate the lipid photooxidation rates
# 
# Jazfulldata.dates = DD_UVB_1314_subsurf$Date[!is.nan(DD_UVB_1314_subsurf$Daily_UVB_dose_0.6m_subsurface_Palmer_kJ_m2)]

# list of dates for which we have subsurface UVB integrals

lipidox.calcdates = DD_UVB_1314_est_at_0.6m_kJ_m2$Date[!is.nan(DD_UVB_1314_est_at_0.6m_kJ_m2$DD_kJ_m2)]

# define some constants

m3_per_L = 0.001
cm_per_m = 100
mol_per_pmol = 1/(10^12)
pmol_per_mol = 10^12

# need package caTools, if not loaded already

detach("package:RSEIS", unload=TRUE)
library(caTools)

#### preallocate results objects for volumetric rate calculations ####

# will make calculations for each of 9 depths, 1-10 m (for which we calculated photon fluxes in UV_TS_analysis_PAL1314.R)

# high-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

# high-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

# high-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

# mid-PUFA fraction, all lipids, UVB
PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

# mid-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

# high-PUFA fraction, all lipids, TUVR
PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
  rownames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
  
}

C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401 =
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}

C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401 = 
  vector("list",length(depths)-1)

for (i in 1:length(C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401)) {
  
  C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)))
  rownames(C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = lipidox.calcdates
  colnames(C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i]]) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_C_L)
  
}


##### calculation of total radiation received at each wavelength, each day, each depth #####

# produces the E_n,p,sigma in Eq. 10 in the manuscript

# preallocate an object to hold the data

PAL1314_E_n_p_sigma_umol_photons_m2_d_nm =
  vector("list",length(depths)-1)

names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm) =
  c(depths[2:length(depths)])

for (i in 1:length(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)) {
  
  PAL1314_E_n_p_sigma_umol_photons_m2_d_nm[[i]] =
    as.data.frame(matrix(NA,length(lipidox.calcdates),length(NOAA_AntUV_lambdas)))
  rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm[[i]]) = lipidox.calcdates
  colnames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm[[i]]) = NOAA_AntUV_lambdas
  
}

for (i in 1:length(lipidox.calcdates)) {
  
  # iterate through dates on which we want an estimate, from lipidox.calcdates
  
  for (j in 2:length(depths)) {
    
    # iterate through depths, 1-10 m
    
    # integrated figures in units of photons/cm2/wavelength
    
    # extract subset of data for this date and depth
    
    PAL1314_est_spectra.today.thisdepth_uW_cm2 =  
      PAL1314_est_spectra_at_depth_uW_cm2[[j]][
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
    PAL1314_est_spectra.today.thisdepth_W_m2 = 
      PAL1314_est_spectra.today.thisdepth_uW_cm2/100
    
    PAL1314_est_spectra.today.thisdepth_umol_photons_m2_s =
      matrix(NA, ncol = length(NOAA_AntUV_lambdas), nrow = nrow(PAL1314_est_spectra.today.thisdepth_uW_cm2))
    
    for (k in 1:ncol(PAL1314_est_spectra.today.thisdepth_umol_photons_m2_s)) {
      
      PAL1314_est_spectra.today.thisdepth_umol_photons_m2_s[,k] =
        PAL1314_est_spectra.today.thisdepth_W_m2[,k]*NOAA_AntUV_lambdas[k]*0.836*(10^-2)
      
    }
    
    # integrate over time, by wavelength
    
    # first, time integration at each wavelength
    
    # preallocate vector
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today.thisdepth = vector(mode = "double",
                                                                      length = ncol(PAL1314_est_spectra.today.thisdepth_umol_photons_m2_s))
    
    for (l in 1:length(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today.thisdepth)) {
      
      # requires vector of times in seconds, with t = 0 being beginning of day
      
      PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today.thisdepth[l] =
        caTools::trapz(Timevector.s[1:length(Timevector.s)],PAL1314_est_spectra.today.thisdepth_umol_photons_m2_s[,l])
      
    }
    
    # store data in appropriate location
    
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm[[which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[j])]][i,] =
      PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today.thisdepth
    
  }
  
}

#### calculate values of F_lipid (Eq. 11 in manuscript) ####

# preallocate some objects to hold results

F_lipid_lambda =
  array(data = NA, dim = c(nrow(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           length(NOAA_AntUV_lambdas)))

dimnames(F_lipid_lambda) =
  list(rownames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       NOAA_AntUV_lambdas)

# make calculations

for (i in 1:length(NOAA_AntUV_lambdas)) {
  
  if (NOAA_AntUV_lambdas[i]<500) { # since don't have any lipid absorbance data above 500 nm
    
    # retrieve appropriate value of Kd 
    
    # retrieve "best" Kd for this wavelength
    
    # find closest wavelength
    
    Kd_ind = which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-NOAA_AntUV_lambdas[i]))
    
    # use Kd derived directly from measurements if it exists; otherwise, use a
    # fitted value
    
    if (!is.na(PAL1516_AH_Kd_20151215_per_meter[Kd_ind,1])) {
      
      Kd.per_meter = PAL1516_AH_Kd_20151215_per_meter[Kd_ind,1]
      
    } else {
      
      Kd.per_meter = PAL1516_AH_Kd_20151215_per_meter[Kd_ind,2]
      
    }
    
    for (j in 1:nrow(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) {
      
      for (k in 1:ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) {
        
        F_lipid_lambda[j,k,i] =
          
          ((LipidAbsData$kappa_M_cm_PC22_6[LipidAbsData$lambda_nm==round(NOAA_AntUV_lambdas[i])]*
              PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[j,k])*cm_per_m*mol_per_pmol)/
          Kd.per_meter
        
      }
      
    }
    
  }
  
}

#### Calculate integrands in Eq. 10 in manuscript ####

# integrands: units of mol photons/wavelength/volume/time (time is implied; t = 1 day)
# once integrated (below), will be in units of mol photons/volume/time

# preallocate objects for the two saturation fractions (high and mid PUFA)

Integrands_mol_photons_m3_nm.hi_PUFA =
  array(data = NA, dim = c(length(lipidox.calcdates),
                           length(NOAA_AntUV_lambdas),
                           ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           length(depths)-1))

dimnames(Integrands_mol_photons_m3_nm.hi_PUFA) =
  list(as.character(lipidox.calcdates),
       NOAA_AntUV_lambdas,
       colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       depths[2:length(depths)])

Integrands_mol_photons_m3_nm.mid_PUFA =
  array(data = NA, dim = c(length(lipidox.calcdates),
                           length(NOAA_AntUV_lambdas),
                           ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           length(depths)-1))

dimnames(Integrands_mol_photons_m3_nm.mid_PUFA) =
  list(as.character(lipidox.calcdates),
       NOAA_AntUV_lambdas,
       colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       depths[2:length(depths)])

for (i in 1:length(NOAA_AntUV_lambdas)) { # iterate by wavelength
  
  if (NOAA_AntUV_lambdas[i]<500) { # since don't have any lipid absorbance data above 500 nm
    
    # retrieve appropriate value of Kd 
    
    # retrieve "best" Kd for this wavelength
    
    # find closest wavelength
    
    Kd_ind = which.min(abs(as.numeric(rownames(PAL1516_AH_Kd_20151215_per_meter))-NOAA_AntUV_lambdas[i]))
    
    # use Kd derived directly from measurements if it exists; otherwise, use a
    # fitted value
    
    if (!is.na(PAL1516_AH_Kd_20151215_per_meter[Kd_ind,1])) {
      
      Kd.per_meter = PAL1516_AH_Kd_20151215_per_meter[Kd_ind,1]
      
    } else {
      
      Kd.per_meter = PAL1516_AH_Kd_20151215_per_meter[Kd_ind,2]
      
    }
    
    for (j in 2:length(depths)) { # iterate by depth
      
      for (k in 1:length(lipidox.calcdates)) { # iterate by date
        
        E_n_p_sigma_photons_m2 = PAL1314_E_n_p_sigma_umol_photons_m2_d_nm[[which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[j])]][k,i]/1000000
        
        for (l in 1:ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) { # iterate by sample
          
          Integrands_mol_photons_m3_nm.hi_PUFA[k,i,l,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[j])] =
            ((E_n_p_sigma_photons_m2*(1-exp(-Kd.per_meter*depths[j])))/depths[j])*F_lipid_lambda[2,l,i]
          
          Integrands_mol_photons_m3_nm.mid_PUFA[k,i,l,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[j])] =
            ((E_n_p_sigma_photons_m2*(1-exp(-Kd.per_meter*depths[j])))/depths[j])*F_lipid_lambda[1,l,i]
          
        }
        
      }
      
    }
    
  }
  
}
    
#### Now, integrate term under integral in Eq. 10 ####

# integrals: units of mol photons/volume/time

# preallocate object to hold calcs for the six scenarios (saturation fractions x wavelength bands)

Integral_mol_photons_m3_d =
  array(data = NA, dim = c(length(lipidox.calcdates),
                           nrow(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
                           3,
                           length(depths)-1))

dimnames(Integral_mol_photons_m3_d) =
  list(as.character(lipidox.calcdates),
       rownames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L),
       c("UVB","UVA","TUVR"),
       depths[2:length(depths)])

for (i in 2:length(depths)) {
  
  for (j in 1:length(lipidox.calcdates)) {
    
      for (l in 1:ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) {
        
        # UVB, hi PUFA fraction

        Integral_mol_photons_m3_d[j,2,l,1,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=315)],
                         Integrands_mol_photons_m3_nm.hi_PUFA[j,
                                                              which(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=315),
                                                              l,
                                                              which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
                         )
        
        # UVB, mid PUFA fraction
        
        Integral_mol_photons_m3_d[j,1,l,1,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=315)],
                         Integrands_mol_photons_m3_nm.mid_PUFA[j,
                                                              which(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=315),
                                                              l,
                                                              which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
          )
        
        # UVA, hi PUFA fraction
        
        Integral_mol_photons_m3_d[j,2,l,2,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>315 & NOAA_AntUV_lambdas<=395.5)],
                         Integrands_mol_photons_m3_nm.hi_PUFA[j,
                                                              which(NOAA_AntUV_lambdas>315 & NOAA_AntUV_lambdas<=395.5),
                                                              l,
                                                              which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
          )
        
        # UVA, mid PUFA fraction
        
        Integral_mol_photons_m3_d[j,1,l,2,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>315 & NOAA_AntUV_lambdas<=395.5)],
                         Integrands_mol_photons_m3_nm.mid_PUFA[j,
                                                               which(NOAA_AntUV_lambdas>315 & NOAA_AntUV_lambdas<=395.5),
                                                               l,
                                                               which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
          )
        
        # TUVR, hi PUFA fraction
        
        Integral_mol_photons_m3_d[j,2,l,3,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=395.5)],
                         Integrands_mol_photons_m3_nm.hi_PUFA[j,
                                                              which(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=395.5),
                                                              l,
                                                              which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
          )
        
        # TUVR, mid PUFA fraction
        
        Integral_mol_photons_m3_d[j,1,l,3,which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])] = 
          caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=395.5)],
                         Integrands_mol_photons_m3_nm.mid_PUFA[j,
                                                               which(NOAA_AntUV_lambdas>=290 & NOAA_AntUV_lambdas<=395.5),
                                                               l,
                                                               which(names(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm)==depths[i])]
          )
        
      }
    
  }
  
}
        
#### Make volumetric lipid oxidation rate calculations for each depth ####

# Eq. 10 in manuscript 

for (i in 2:length(depths)) {
  
  for (j in 1:length(lipidox.calcdates)) {
    
    for (k in 1:ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) {
      
      # UVB, hi PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,1,i-1]*Theta_Exp13_PC_22_6[3]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,1,i-1]*Theta_Exp13_PC_22_6_sigma[3]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      C_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      # UVB, mid PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,1,i-1]*Theta_Exp13_PC_22_6[3]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,1,i-1]*Theta_Exp13_PC_22_6_sigma[3]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
      C_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
      # UVA, hi PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,2,i-1]*Theta_Exp13_PC_22_6[2]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,2,i-1]*Theta_Exp13_PC_22_6_sigma[2]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      # UVA, mid PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,2,i-1]*Theta_Exp13_PC_22_6[2]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,2,i-1]*Theta_Exp13_PC_22_6_sigma[2]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
      C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
      # TUVR, hi PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,3,i-1]*Theta_Exp13_PC_22_6[1]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,2,k,3,i-1]*Theta_Exp13_PC_22_6_sigma[1]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      C_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[2,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[2,k])
      
      # TUVR, mid PUFA fraction
      
      PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,3,i-1]*Theta_Exp13_PC_22_6[1]*
        pmol_per_mol*m3_per_L
      
      PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        Integral_mol_photons_m3_d[j,1,k,3,i-1]*Theta_Exp13_PC_22_6_sigma[1]*
        pmol_per_mol*m3_per_L
      
      C_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
      C_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k]*
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
    }
    
  }
  
}

#### perform integration over depth of mixed layer ####

# *** from here, just carrying forward results using UVA-band AQY ***

# preallocate objects

# high-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401 =
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401) = lipidox.calcdates
colnames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 =
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401) = lipidox.calcdates
colnames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401 = 
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401) = lipidox.calcdates
colnames(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 =
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401) = lipidox.calcdates
colnames(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

# mid-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401 =
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401) = lipidox.calcdates
colnames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 =
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401) = lipidox.calcdates
colnames(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401 = 
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401) = lipidox.calcdates
colnames(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)
C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 = 
  as.data.frame(matrix(NA,length(lipidox.calcdates),ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)))
rownames(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401) = lipidox.calcdates
colnames(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401) = colnames(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)

for (i in 1:length(lipidox.calcdates)) {
  
  # define depth of ML to be used; see manuscript Section 3.1, also spreadsheet
  # "MLD Calculations" in LipidPhotoOxBox/data/raw/PAL1314_LMG1401_PAL_LTER_data/Selected CTD - PAL-LTER 2013-2014.xlsx
  
  if (lipidox.calcdates[i]<as.POSIXct('2013-12-18 00:00:00', tz = "GMT")) {
    
    MLD.AH = 7 # MLD at Station B on 12/12/13
    
  } else if (lipidox.calcdates[i]>=as.POSIXct('2013-12-18 00:00:00', tz = "GMT") &
             lipidox.calcdates[i]<as.POSIXct('2013-12-25 00:00:00', tz = "GMT")
             ) {
    
    MLD.AH = 5 # MLD at Station B on 12/14/13
    
  } else if (lipidox.calcdates[i]>=as.POSIXct('2013-12-25 00:00:00', tz = "GMT")) {
    
    MLD.AH = 10 # MLD at Station B on 12/27/13
    
  }
  
  for (j in 1:ncol(PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L)) {
    
    PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401[i,j] = 
      caTools::trapz(seq(1:MLD.AH), 
                     unlist(lapply(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))*
                       L_per_m3)
    
    PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH), 
                     unlist(lapply(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))*
                       L_per_m3)
    
    C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH),
                     unlist(lapply(C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))* 12.01 * (1/10^6) * 10^3)
    
    C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH),
                     unlist(lapply(C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))* 12.01 * (1/10^6) * 10^3)
    
    PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401[i,j] = 
      caTools::trapz(seq(1:MLD.AH), 
                     unlist(lapply(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))*
                       L_per_m3)
    
    PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH), 
                     unlist(lapply(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))*
                       L_per_m3)
    
    C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH),
                     unlist(lapply(C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))* 12.01 * (1/10^6) * 10^3)
    
    C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[i,j] =
      caTools::trapz(seq(1:MLD.AH),
                     unlist(lapply(C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[1:MLD.AH], function(x) x[i,j]))* 12.01 * (1/10^6) * 10^3)
    
  }
  
}

#### pull in, work up PAL1314 BP data, for comparison ####

PAL1314.BP = read.csv("/Users/jamesrco/Code/LipidPhotoOxBox/data/raw/PAL1314_LMG1401_PAL_LTER_data/PAL1314 Bacterial Production (Station).csv", skip = 1, stringsAsFactors = FALSE)

# change missing values to NA
PAL1314.BP$Leucine.Incorp...pmol.L.hr.[PAL1314.BP$Leucine.Incorp...pmol.L.hr.==-999.00] = NA

# convert leucine numbers to units of ug C/m3/d, using most conservative assumptions (ID = 1, conversion factor of 1.5 kg C/mol leu)

PAL1314.BP$BP_ugC_m3_d = PAL1314.BP$Leucine.Incorp...pmol.L.hr. * 36

# append a vector of dates in POSIXct format

PAL1314.BP$Date.POSIXct = as.POSIXct(PAL1314.BP$Date.GMT, format = "%m/%d/%y", tz = "GMT")
  
# calculate some depth-integrated BP rates for the mixed layer

# load in MLD data

PAL1314.MLD.Depths = read.csv("/Users/jamesrco/Code/LipidPhotoOxBox/data/raw/PAL1314_LMG1401_PAL_LTER_data/PAL1314 - MLD Depths.csv", skip = 3, stringsAsFactors = FALSE)

# append a vector of dates in POSIXct format

PAL1314.MLD.Depths$Date.POSIXct = as.POSIXct(PAL1314.MLD.Depths$Date, format = "%m/%d/%y", tz = "GMT")

# obtain list of dates <= 1/2/14 for which we want to calculate depth-integrated BP

PAL1314.BP.depthint.dates =
  unique(PAL1314.BP$Date.POSIXct)[unique(PAL1314.BP$Date.POSIXct)<=
                                    as.POSIXct('2014-01-02', tz = "GMT")]

# preallocate object for depth int BP calcs

PAL1314.BP.depthint.ugC_m2_d =
  as.data.frame(matrix(NA,length(PAL1314.BP.depthint.dates),4))
PAL1314.BP.depthint.ugC_m2_d[,1] = PAL1314.BP.depthint.dates
rownames(PAL1314.BP.depthint.ugC_m2_d) = PAL1314.BP.depthint.dates
colnames(PAL1314.BP.depthint.ugC_m2_d) = c("Date","SWI","B","E")

# perform calculations

for (i in 1:nrow(PAL1314.BP.depthint)) {
  
  for (j in 2:ncol(PAL1314.BP.depthint)) {
    
    # determine whether we have data for this station on this date
    
    BPdat.subset.ugC_m3_d = 
      PAL1314.BP[PAL1314.BP$Date.POSIXct==PAL1314.BP.depthint.ugC_m2_d$Date[i] &
                   PAL1314.BP$Station.Name==colnames(PAL1314.BP.depthint.ugC_m2_d)[j],
                 c("Depth..m.","BP_ugC_m3_d")]
    
    if (nrow(BPdat.subset.ugC_m3_d)>0) { # we have data, proceed
      
      # eliminate any entries for value = NA
      
      BPdat.subset.ugC_m3_d = BPdat.subset.ugC_m3_d[!(is.na(BPdat.subset.ugC_m3_d$BP_ugC_m3_d)),]
      
      # need to average the data if from SWI since replicates were taken
      
      if (colnames(PAL1314.BP.depthint.ugC_m2_d)[j]=="SWI") {
        
        SWI.BP.mean = mean(BPdat.subset.ugC_m3_d$BP_ugC_m3_d)
        BPdat.subset.ugC_m3_d$Depth..m. = BPdat.subset.ugC_m3_d$Depth..m.[1]
        BPdat.subset.ugC_m3_d$BP_ugC_m3_d = SWI.BP.mean
        BPdat.subset.ugC_m3_d = BPdat.subset.ugC_m3_d[1,]
        
      }
      
      # obtain most relevant MLD; using Station B values for SWI
      
      # pull out possible values
      
      if (colnames(PAL1314.BP.depthint.ugC_m2_d)[j]=="SWI") {
        
        MLDs.thisstation = PAL1314.MLD.Depths[PAL1314.MLD.Depths$Station.Name=="B",]
        
      } else {
        
        MLDs.thisstation = PAL1314.MLD.Depths[PAL1314.MLD.Depths$Station.Name==
                                                colnames(PAL1314.BP.depthint.ugC_m2_d)[j],]
        
      }
      
      # select closest value temporally
      
      relevant.MLD.m = MLDs.thisstation$MLD.depth.m[
        which.min(abs(PAL1314.BP.depthint.ugC_m2_d$Date[i]-MLDs.thisstation$Date.POSIXct))]
      
      # calculate depth-integrated BP
      
      # assume uniform values at all depths in ML for SWI case
      
      if (colnames(PAL1314.BP.depthint.ugC_m2_d)[j]=="SWI") {
        
        PAL1314.BP.depthint.ugC_m2_d[i,j] =
          caTools::trapz(c(0,relevant.MLD.m),rep(BPdat.subset.ugC_m3_d$BP_ugC_m3_d,2))
        
      } else {
        
        if (relevant.MLD.m %in% BPdat.subset.ugC_m3_d$Depth..m.) {
          
          BPdat.subset.ugC_m3_d.ML_only = 
            BPdat.subset.ugC_m3_d[BPdat.subset.ugC_m3_d$Depth..m.<=relevant.MLD.m,]
          
        } else {
          
          # need to estimate BP by linear interpolation at depth of ML
          
          # obtain BP observations for the two depths bracketing the MLD
          
          deeper.MLDs.ind = which((BPdat.subset.ugC_m3_d$Depth..m.-relevant.MLD.m)>0)
          next.deepest.from.target = min(BPdat.subset.ugC_m3_d$Depth..m.[deeper.MLDs.ind])
          
          shallower.MLDs.ind = which((BPdat.subset.ugC_m3_d$Depth..m.-relevant.MLD.m)<0)
          next.shallowest.from.target = max(BPdat.subset.ugC_m3_d$Depth..m.[shallower.MLDs.ind])
          
          BP.bracket = c(BPdat.subset.ugC_m3_d$BP_ugC_m3_d[BPdat.subset.ugC_m3_d$Depth..m.==next.shallowest.from.target],
                         BPdat.subset.ugC_m3_d$BP_ugC_m3_d[BPdat.subset.ugC_m3_d$Depth..m.==next.deepest.from.target])
          
          depth.bracket = c(next.shallowest.from.target,next.deepest.from.target)
          
          # predict BP at MLD by linear interpolation
          
          linfit.BP = lm(BP.bracket~depth.bracket)
          
          BP.est.at.MLD_ugC_m3_d = predict.lm(linfit.BP, data.frame(depth.bracket = relevant.MLD.m))
          
          # construct subset of needed values
          
          BPdat.subset.ugC_m3_d.ML_only =
            rbind(BPdat.subset.ugC_m3_d[BPdat.subset.ugC_m3_d$Depth..m.<relevant.MLD.m,],
                  as.numeric(data.frame(relevant.MLD.m,BP.est.at.MLD_ugC_m3_d)))
          
        }
        
        # now, calculate depth-int BP; store
        
        # reorder data frame
        
        BPdat.subset.ugC_m3_d.ML_only = 
          BPdat.subset.ugC_m3_d.ML_only[order(BPdat.subset.ugC_m3_d.ML_only$Depth..m.),]
        
        PAL1314.BP.depthint.ugC_m2_d[i,j] =
          caTools::trapz(BPdat.subset.ugC_m3_d.ML_only$Depth..m.,
                         BPdat.subset.ugC_m3_d.ML_only$BP_ugC_m3_d)
        
      }
      
    }
    
  }
  
}

# make a plot; use most proximate lipid concentration data (temporally)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Lipid_AQY_results_w_BP.pdf",
    width = 6.5, height = 3.5, pointsize = 10,
    bg = "white")

par(mar=c(5,5,1,5))

plot(lipidox.calcdates,
     C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000,
     "o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(0,10),
     xlim = c(as.numeric(lipidox.calcdates[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "mg C per m2 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
     )

axis.POSIXct(1, at = seq(lipidox.calcdates[1], as.POSIXct("2014-01-06"), by = "5 days"), format = "%d %b")

polygon(c(lipidox.calcdates,
          rev(lipidox.calcdates)),
        c(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000+
            C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000,
          rev(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000-
                C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000)),
        border = NA, col = "mistyrose")

polygon(c(lipidox.calcdates,
          rev(lipidox.calcdates)),
        c(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000+
            C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000,
          rev(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000-
                C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000)),
        border = NA, col = "azure")

lines(lipidox.calcdates,C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000,"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred")

lines(lipidox.calcdates,C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000, lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

# superimpose BP data

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$E)],
       PAL1314.BP.depthint.ugC_m2_d$E[!is.na(PAL1314.BP.depthint.ugC_m2_d$E)]/1000,
       pch = 22, bg = "black", cex = 1.2)

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$SWI)],
       PAL1314.BP.depthint.ugC_m2_d$SWI[!is.na(PAL1314.BP.depthint.ugC_m2_d$SWI)]/1000,
       pch = 22, bg = "white", cex = 1.2)

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$B)],
       PAL1314.BP.depthint.ugC_m2_d$B[!is.na(PAL1314.BP.depthint.ugC_m2_d$B)]/1000,
       pch = 21, bg = "black", cex = 1.2)

dev.off()

# second plot with the higher-value BP samples

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Lipid_AQY_results_high_BP.pdf",
    width = 6.5, height = 3.5, pointsize = 10,
    bg = "white")

par(mar=c(5,5,1,5))

plot(lipidox.calcdates,
     C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000,
     "o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(10,45),
     xlim = c(as.numeric(lipidox.calcdates[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "mg C per m2 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
)

# superimpose BP data

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$E)],
       PAL1314.BP.depthint.ugC_m2_d$E[!is.na(PAL1314.BP.depthint.ugC_m2_d$E)]/1000,
       pch = 22, bg = "black", cex = 1.2)

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$SWI)],
       PAL1314.BP.depthint.ugC_m2_d$SWI[!is.na(PAL1314.BP.depthint.ugC_m2_d$SWI)]/1000,
       pch = 22, bg = "white", cex = 1.2)

points(PAL1314.BP.depthint.ugC_m2_d$Date[!is.na(PAL1314.BP.depthint.ugC_m2_d$B)],
       PAL1314.BP.depthint.ugC_m2_d$B[!is.na(PAL1314.BP.depthint.ugC_m2_d$B)]/1000,
       pch = 21, bg = "black", cex = 1.2)

dev.off()

# last plot, with alternate y-axis

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Lipid_AQY_results_w_BP_lipid_units.pdf",
    width = 6.5, height = 3.5, pointsize = 10,
    bg = "white")

par(mar=c(5,5,1,5))

plot(lipidox.calcdates,
     PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000,
     "o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(0,12),
     yaxp = c(0,12,24),
     xlim = c(as.numeric(lipidox.calcdates[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "umol IP-DAG per m2 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
)

axis.POSIXct(1, at = seq(lipidox.calcdates[1], as.POSIXct("2014-01-06"), by = "5 days"), format = "%d %b")

polygon(c(lipidox.calcdates,
          rev(lipidox.calcdates)),
        c(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000+
            PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000,
          rev(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000-
                PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000)),
        border = NA, col = "mistyrose")

polygon(c(lipidox.calcdates,
          rev(lipidox.calcdates)),
        c(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000+
            PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000,
          rev(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000-
                PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000)),
        border = NA, col = "azure")

lines(lipidox.calcdates,PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000,"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred")

lines(lipidox.calcdates,PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120/1000000, lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

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

apply(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ML_ug_C_m2_d.UVA.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_ML_pmol_m2_d.UVA.hi_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ML_ug_C_m2_d.UVA.hi_PUFA.PAL1314_LMG1401,2,sd)

# BP averages, untuil 1/2/14

mean(PAL1314.BP.depthint.ugC_m2_d$B, na.rm = T)
sd(PAL1314.BP.depthint.ugC_m2_d$B, na.rm = T)

mean(PAL1314.BP.depthint.ugC_m2_d$E, na.rm = T)
sd(PAL1314.BP.depthint.ugC_m2_d$E, na.rm = T)

mean(PAL1314.BP.depthint.ugC_m2_d$SWI, na.rm = T)
sd(PAL1314.BP.depthint.ugC_m2_d$SWI, na.rm = T)

