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

#### preallocate results objects ####

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
        PUFA_lipids_xformed_pmol_L_d.UVB.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
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
        PUFA_lipids_xformed_pmol_L_d.UVB.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
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
        PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
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
        PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
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
        PUFA_lipids_xformed_pmol_L_d.TUVR.hi_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
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
        PUFA_lipids_xformed_pmol_L_d.TUVR.mid_PUFA.sigma.PAL1314_LMG1401[[i-1]][j,k] =
        (PUFA_fracs.PAL1314_LMG1401.pmol_C_L[1,k]/
           PUFA_fracs.PAL1314_LMG1401.pmol_lipids_L[1,k])
      
    }
    
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

# high-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401 
PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 
C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3
C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3

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

# mid-PUFA fraction, all lipids, UVA
PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401
PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 
C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3
C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 =
  C_xformed_pmol_L_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401 * 12.01 * (1/10^6) * 10^3

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

C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401.nogaps) =
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

C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps) =
  Jazfulldata.dates.nogaps

C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401.nogaps = as.data.frame(matrix(NA,length(Jazfulldata.dates.nogaps),ncol(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401)))
colnames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401.nogaps) =
  colnames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401)
rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401.nogaps) =
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
    
    C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401),]
    
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
    
    C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401),]
    
    C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401.nogaps[i,] =
      C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401[
        rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401.nogaps)[i]==
          rownames(C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401),]
    
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

plot(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
     ylim = c(0,1000),
     xlim = c(as.numeric(Jazfulldata.dates.nogaps[1]),as.numeric(as.POSIXct("2014-01-05"))),
     ylab = "ug C per m3 per day",
     xlab = "Date (2013-2014)",
     xaxt = "n"
     )

axis.POSIXct(1, at = seq(Jazfulldata.dates.nogaps[1], as.POSIXct("2014-01-06"), by = "5 days"), format = "%d %b")

polygon(c(Jazfulldata.dates[1:27],
          rev(Jazfulldata.dates[1:27])),
        c(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]+
            C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27],
          rev(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]-
                C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27])),
        border = NA, col = "mistyrose")

polygon(c(Jazfulldata.dates[27:35],
          rev(Jazfulldata.dates[27:35])),
        c(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]+
            C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35],
          rev(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]-
                C_xformed_ug_C_m3_d.UVA.hi_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35])),
        border = NA, col = "mistyrose")

polygon(c(Jazfulldata.dates[1:27],
          rev(Jazfulldata.dates[1:27])),
        c(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]+
            C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27],
          rev(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27]-
                C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:27])),
        border = NA, col = "azure")

polygon(c(Jazfulldata.dates[27:35],
          rev(Jazfulldata.dates[27:35])),
        c(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]+
            C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35],
          rev(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401$LMG1401_004_QE003116[27:35]-
                C_xformed_ug_C_m3_d.UVA.mid_PUFA.sigma.PAL1314_LMG1401$LMG1401_004_QE003116[27:35])),
        border = NA, col = "azure")

lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred")

lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], type = "o", pch = 22, bg = "black", cex = 0.4, col = "darkred")
# lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40], lty=2, type = "o", pch = 22, bg = "white")
# lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.UVB.hi_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], lty=2, type = "o", pch = 22, bg = "white")

lines(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40], lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

lines(Jazfulldata.dates.nogaps[41:63],C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401.nogaps$LMG1401_004_QE003116[41:63], lty=3, type = "o", pch = 21, bg = "cyan", cex = 0.4, col = "cyan")

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

plot(Jazfulldata.dates.nogaps[1:40],C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401.nogaps$PAL1314_Stn_E_3m_2Jan14_0.2um_QE003120[1:40],"o", pch = 22, bg = "darkred", cex = 0.4, col = "darkred",
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

apply(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_pmol_L_d.UVA.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ug_C_m3_d.UVA.mid_PUFA.PAL1314_LMG1401,2,sd)

apply(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(PUFA_lipids_xformed_pmol_L_d.UVA.hi_PUFA.PAL1314_LMG1401,2,sd)

apply(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401,2,mean)
apply(C_xformed_ug_C_m3_d.UVA.hi_PUFA.PAL1314_LMG1401,2,sd)

# BP averages, untuil 1/2/14

mean(BP.subset.B$BP_ugC_m3_d[1:5])
sd(BP.subset.B$BP_ugC_m3_d[1:5])

mean(BP.subset.E$BP_ugC_m3_d[1:2])
sd(BP.subset.E$BP_ugC_m3_d[1:2])

mean(BP.subset.SWI.means$BP_ugC_m3_d[1:5])
sd(BP.subset.SWI.means$BP_ugC_m3_d[1:5])

